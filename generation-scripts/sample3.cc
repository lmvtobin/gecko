/**
* * * * * * * * * * * * * * * * * * * * *
*                           o---------x *
*     G e c K o         x-o |  #  w w | *
*                     o-|-x |       w | *
*       Grid          x-|---|---o     | *
*       Kirby     x-----o   |   |   # | *
*                 |     x---|---|-----o *
*             x---o     |   |   | x-o   *
*             |     o---|---x o-|-|-x   *
*         o-x |     |   |     | x-o     *
*       x-|-|-o x---|---o     |         *
*       o-x |   |   |     o---x         *
* o---------|---|---x     |             *
* |         |   |   o-----x             *
* |     o---|---|---|-x                 *
* |     |   |   | x-|-o        By       *
* |   x-|---o   | o-x         Lucy      *
* | o-|-x       |             Tobin     *
* x-|-o         |                       *
*   x-----------o                       *
* * * * * * * * * * * * * * * * * * * * *
**/

/**
 * Command line utility which generates a sample of random 3-manifolds, simplifies the
 * resulting triangulations as much as possible, and performs connected sum decompositions.
 * 
 * Connected sum components are guaranteed to be prime, but not necessarily minimal triangulations.
 * 
 * Results will be printed to console, and optionally saved to
 * ../data/sample3/nN_S_I.txt or ../data/sample3/knots_nN_S_I.txt
 * where
 * - N is the grid size
 * - S is the sample size
 * - I = 0,1,2,... is a counter which increments if a file with this N and S already exists
 * 
 * Save format is one line per simplified isomorphism signature / combinatorial type, of the form:
 * sig1 sig2 ... sigK : num_grids : grid_1 grid_2 ... grid_(num_grids)
 * where
 * - sig1,... are isomorphism signatures for triangulations of the prime summands
 * - num_grids is the number of grids whose triangulations were simplified to this
 * - grid_1,... are grids as pairs of permutations, in the form e.g. [0,1,2][1,2,0]
 * 
 * Note the lines are not guaranteed to represent distinct 3-manifolds.
 **/

#include <string>
#include <vector>
#include <queue>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <ranges>

#include <thread>
#include <mutex>
#include <condition_variable>
#include <semaphore>

#include <boost/process.hpp>
#include <boost/asio.hpp>

#include <census/census.h>
#include <triangulation/isosigtype.h>
#include <snappea/snappeatriangulation.h>
#include <subcomplex/standardtri.h>

#include "../gecko/geckocore.h"
#include "../mcmc/mcmc3.h"

std::queue<std::pair<std::string,int>> countQueue;
std::queue<std::string> sigQueue;
std::mutex queue_mtx;
std::condition_variable queue_cv;
std::counting_semaphore<20> queue_sem {20};
bool producerDone = false;

std::mutex output_mtx;

std::mutex count_mtx;
std::map<std::string,int> finalCounts;
std::map<std::string,int> unidentifiedCounts;
std::map<std::vector<std::string>,int> decomposedCounts;
std::map<std::string,int> couldntDecomposeCounts;
std::map<std::string,std::vector<gecko::Grid>> gridMap;

std::map<std::string,int> heightChecked;
std::map<std::string,std::string> hardComponentsMap;

int timeout_mins = 1440;


void sampleWorker(int n, int repeats, bool knotOnly) {

    {
        std::unique_lock<std::mutex> lck(output_mtx);
        std::clog << "worker started\n";
    }

    std::vector<std::pair<regina::Triangulation<3>,gecko::Grid>> samples = gecko::sample3Manifolds(n,repeats,knotOnly);

    std::map<std::string,int> counts;
    std::map<std::string,std::vector<gecko::Grid>> grids;

    for (auto [trig,grid] : samples) {
        for (int i=0; i<10; ++i) {
            trig.simplify();
        }
        std::string sig = trig.isoSig<regina::IsoSigDegrees<3,1>>();
        ++counts[sig];
        if (grids.contains(sig)) {
            grids[sig].push_back(grid);
        } else {
            grids[sig] = {grid};
        }
    }

    {
        std::scoped_lock<std::mutex> lck(count_mtx);
        for (auto [sig,count] : counts) {
            unidentifiedCounts[sig]+=count;
            if (gridMap.contains(sig)) {
                gridMap[sig].insert(gridMap[sig].end(),grids[sig].begin(),grids[sig].end());
            } else {
                gridMap[sig] = grids[sig];
            }
        }
    }

    {
        std::unique_lock<std::mutex> lck(output_mtx);
        std::clog << "worker done\n";
    }

}

void decomposeConsumer(int workerID) {

    {
        std::unique_lock<std::mutex> lck(output_mtx);
        std::clog << "worker " << workerID << " started\n";
    }

    auto workerStart = std::chrono::high_resolution_clock::now();

    std::pair<std::string,int> pair;

    while (true) {

        {
            // wait for countQueue to be nonempty and get the lock
            // once producer has loaded all grids, continue until countQueue is empty
            std::unique_lock<std::mutex> lck(queue_mtx);
            queue_cv.wait(lck, []() { return !countQueue.empty() or producerDone; });

            if (countQueue.empty()) {
                // all sigs generated and queue empty: done!
                break;
            } else {
                // get grid
                pair = countQueue.front();
                countQueue.pop();
                queue_sem.release();
            }
        }

        {
            std::unique_lock<std::mutex> lck(output_mtx);
            std::cout << "worker " << workerID << " is decomposing " << pair.first << "\n";
        }

        regina::Triangulation<3> trig(pair.first);
        trig.simplify();

        // try {

        // std::vector<regina::Triangulation<3>> summands = summandWrapper(trig);
        boost::process::ipstream decomp_stream;
        boost::process::child decomposerProc("./summands-wrapper", trig.isoSig<regina::IsoSigDegrees<3,1>>(), boost::process::std_out > decomp_stream);
        std::vector<regina::Triangulation<3>> summands;
        std::string line;
        // bool complete = decomposerChild.wait_for(std::chrono::minutes(timeout_mins));
        // decomposerProc.wait();
        // boost::asio::execute(decomposerProc,(boost::asio::cancel_after(std::chrono::minutes(timeout_mins)),boost::asio::cancellation_type::terminal));

        bool complete = false;
        int seconds_waited = 0;
        int check_interval = 1;
        // this is an incedibly painful way of implementing a timeout and I hate it
        // boost's wait_for function is broken in v1.74
        // a sensible person would update boost and use asio::cancel_after
        // but I don't want to update ubuntu to update boost, so here we are
        for (int i=0; i<10; ++i) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            if (!decomposerProc.running()) {
                complete = true;
                break;
            }
        }
        seconds_waited += 1;
        if (!complete) {
            while (true) {
                std::this_thread::sleep_for(std::chrono::seconds(check_interval));
                seconds_waited += check_interval;
                if (!decomposerProc.running()) {
                    complete = true;
                    break;
                }
                if (seconds_waited >= 60*timeout_mins) {
                    decomposerProc.terminate();
                    decomposerProc.wait();
                    complete = false;
                    break;
                }
            }
        }

        if (!complete | (complete && decomposerProc.exit_code())) {
            {
                std::scoped_lock<std::mutex> lck(count_mtx);
                couldntDecomposeCounts[pair.first] = pair.second;
            }

            {
                std::unique_lock<std::mutex> lck(output_mtx);
                std::cout << "--> worker " << workerID << ": " << (complete ? "encountered an error" : "timed out" ) << "\n";
            }   
            continue;
        }

        while (decomp_stream && std::getline(decomp_stream,line) && !line.empty()) {
            summands.push_back(regina::Triangulation<3>(line));
        }

        std::vector<std::string> decomp;
        if (summands.empty()) {
            // trig is S3
            decomp.push_back("bkaagj");
        } else {
            for (regina::Triangulation<3> summand : summands) {
                summand.simplify();
                decomp.push_back(summand.isoSig<regina::IsoSigDegrees<3,1>>());
            }
        }
        std::sort(decomp.begin(),decomp.end());

        std::string fullname;
        for (std::string sig : decomp) {
            fullname = fullname + sig + " ";
        }
        fullname.erase(fullname.size()-1,1); // erase last " "

        {
            std::scoped_lock<std::mutex> lck(count_mtx);
            decomposedCounts[decomp]+=pair.second;
            auto nh = gridMap.extract(pair.first);
            if (gridMap.contains(fullname)) {
                gridMap[fullname].insert(gridMap[fullname].end(),nh.mapped().begin(),nh.mapped().end());
            } else {
                gridMap[fullname] = nh.mapped();
            }
        }

        {
            std::unique_lock<std::mutex> lck(output_mtx);
            std::cout << "--> worker " << workerID << ": ";
            for (std::string s : decomp) {
                std::cout << s << " ";
            }
            std::cout << "\n";
        }

    }

    {
        std::unique_lock<std::mutex> lck(output_mtx);
        std::clog << "worker " << workerID << " done | ";
        auto workerNow = std::chrono::high_resolution_clock::now();
        std::cout << "active for " << std::chrono::duration_cast<std::chrono::seconds>(workerNow-workerStart).count() << " s\n\n";
    }

}

void saveToFile(std::map<std::string,int> counts, std::string baseFilename)
{
    std::filesystem::create_directories("../data/sample3"); // create storage/output directory (if it doesn't exist)
    std::string filename = baseFilename;
    int filenum = 0;
    while (true) {
        filename = baseFilename + "_" + std::to_string(filenum);
        if (!std::filesystem::exists("../data/sample3/"+filename+".txt")) {
            break;
        } else {
            ++filenum;
        }
    }
    std::ofstream writer;
    writer.open("../data/sample3/"+filename+".txt");
    for (auto pair : counts)
    {
        writer << pair.first << " : " << pair.second << " : ";
        for (gecko::Grid grid : gridMap[pair.first]) {
            gecko::printVector(grid.x,writer,false);
            gecko::printVector(grid.o,writer,false);
            writer << " ";
        }
        writer << "\n";
    }
    writer.close();

    std::cout << "Saved to ../data/sample3/" << filename << ".txt\n";
}

int secondsSince(std::chrono::time_point<std::chrono::high_resolution_clock> start) {
    return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-start).count();
}

int main(int argc, char* argv[]) {
    
    if (argc < 2) {
        std::cerr << "Usage:\n ./sample3 gridsize [samples=1] [cores=1] [timeout_mins=1440] [-save] [-knot]\n";
        std::cerr << "Example usage for grid size 6 with default values and not saving:\n";
        std::cerr << "./sample3 6\n";
        std::cerr << "Example usage for grid size 6, 50 samples on 10 cores, timing out after 15 minutes of trying to decompose, saving to default location, knots only:\n";
        std::cerr << "./sample3 6 50 10 15 -save -knot\n";
        std::cerr << "Optional positional arguments must be in order, flags (-save, -knot) must appear after including all optional arguments\n";
        return 1;
    }

    auto start = std::chrono::high_resolution_clock::now();

    int n = std::stoi(argv[1]);
    int repeats = 1;
    int cores = 1;
    bool hard = false;
    bool saveResults = false;
    bool knotOnly = false;
    if (argc > 2) {
        repeats = std::stoi(argv[2]);
    }
    if (argc > 3) {
        cores = std::stoi(argv[3]);
    }
    if (argc > 4) {
        timeout_mins = std::stoi(argv[4]);
    }
    for (int i=5; i<argc; ++i) {
        if (std::string(argv[i])=="-save") {
            saveResults = true;
        } else if (std::string(argv[i])=="-knot") {
            knotOnly = true;
        }
    }


    std::cout << "--- (1.1) collect iso sigs ---\n\n";

    std::vector<int> repeatsPerWorker(cores);

    bool done = false;
    int counter = 0;
    while (!done) {
        for (int i=0; i<cores; ++i) {
            if (counter==repeats) {
                done = true;
                break;
            } else {
                ++repeatsPerWorker[i];
                ++counter;
            }
        }
    }

    std::cout << "Sample size per worker: ";
    gecko::printVector(repeatsPerWorker,std::cout);

    std::vector<std::thread> workers;
    for (int i=0; i<cores; ++i) {
        workers.emplace_back(std::thread(sampleWorker,n,repeatsPerWorker[i],knotOnly));
    }

    for (auto& worker : workers) {
        worker.join();
    }
    workers.clear();

    std::cout << "\n";

    std::cout << "elapsed: " << secondsSince(start) << " s\n\n\n";


    std::cout << "--- (1.2) census lookup --- \n\n";

    std::map<std::string,int> newUnidentifiedCounts;

    for (auto pair : unidentifiedCounts) {
        std::string sig = pair.first;
        regina::Triangulation<3> trig(sig);
        std::list<regina::CensusHit> lookup = regina::Census::lookup(trig);
        auto recog = regina::StandardTriangulation::recognise(trig);
        if (!recog && lookup.empty()) {
            newUnidentifiedCounts[sig]+=pair.second;
        } else {
            // our triangulation is in a census, so we assume it is minimal (or near-minimal)
            // and stop trying to simplify it
            // I believe Regina's censuses don't contain any connected sums, but just in case
            // we confirm it is irreducible here (and we don't need to decompose it)
            if (  (!lookup.empty() && lookup.front().name().starts_with("S2 x S1 : #")) |
                  (recog && recog->manifold() && recog->manifold()->name()=="S2 x S1") |
                  trig.isIrreducible() ) {
                finalCounts[sig]+=pair.second;
            } else {
                newUnidentifiedCounts[sig]+=pair.second;
            }
        }
    }

    unidentifiedCounts = newUnidentifiedCounts;
    newUnidentifiedCounts.clear();

    std::cout << "Identified near-minimal triangulations:\n\n";
    if (!finalCounts.empty()) {
        for (auto pair : finalCounts) {
            std::cout << pair.first << ": " << pair.second << "\n";
        }
        std::cout << "\n";
    }

    int checkTotal = 0;
    for (auto pair : finalCounts) {
        checkTotal += pair.second;
    }
    std::cout << "Remaining: " << repeats-checkTotal << " / " << repeats << "\n\n";

    std::cout << "elapsed: " << secondsSince(start) << " s\n\n\n";


    std::cout << "--- (2.1) MCMC simplification, first pass ---\n\n";

    int countSimplified = 0;
    std::map<std::string,int> toSimplify = unidentifiedCounts;
    std::map<std::string,int> newToSimplify;

    for (auto pair : toSimplify) {
        std::string sig = pair.first;

        std::vector<std::future<std::string>> workers;
        for (int i=0; i<cores; ++i) {
            workers.emplace_back(std::async(mcmc::simplify,sig,100,10,false));
        }
        std::string simplifiedSig = workers[0].get();
        for (int i=1; i<cores; ++i) {
            std::string newsig = workers[i].get();
            if (!newsig.empty() && newsig.size() < simplifiedSig.size()) {
                simplifiedSig = newsig;
            }
            if (newsig.empty()) {std::cout << "!";}
        }
        workers.clear();

        if (!simplifiedSig.empty() && simplifiedSig.size() < sig.size()) {
            ++countSimplified;
            // if we've previously simplified something into this sig and then
            // also simplify it again here, apply the extra simplification to
            // the previous too
            // otherwise, sig ends up in newUnidentifiedCounts but ~not~ in gridMap
            // which causes errors down the line
            if (newUnidentifiedCounts.contains(sig)) {
                auto nh = newUnidentifiedCounts.extract(sig);
                newUnidentifiedCounts[simplifiedSig] += nh.mapped();
            }
            newUnidentifiedCounts[simplifiedSig]+=pair.second;
            hardComponentsMap[sig]=simplifiedSig;
            auto nh = gridMap.extract(sig);
            if (gridMap.contains(simplifiedSig)) {
                gridMap[simplifiedSig].insert(gridMap[simplifiedSig].end(),nh.mapped().begin(),nh.mapped().end());
            } else {
                gridMap[simplifiedSig] = nh.mapped();
            }
            std::cout << "S" << std::flush;
        } else {
            newToSimplify[sig]+=pair.second;
            std::cout << "F" << std::flush;
        }
    }

    std::cout << "\n\n";

    toSimplify = newToSimplify;
    newToSimplify.clear();

    std::cout << "Simplified " << countSimplified << " triangulations\n\n";

    std::cout << "elapsed: " << secondsSince(start) << " s\n\n\n";


    std::cout << "--- (2.2) MCMC simplification, second pass ---\n\n";

    countSimplified = 0;

    for (auto pair : toSimplify) {
        std::string sig = pair.first;
        std::vector<std::future<std::string>> workers;
        for (int i=0; i<cores; ++i) {
            workers.emplace_back(std::async(mcmc::simplify,sig,1000,10,false));
        }
        std::string simplifiedSig = workers[0].get();
        for (int i=1; i<cores; ++i) {
            std::string newsig = workers[i].get();
            if (!newsig.empty() &&newsig.size() < simplifiedSig.size()) {
                simplifiedSig = newsig;
            }
            if (newsig.empty()) {std::cout << "!";}
        }
        workers.clear();

        if (!simplifiedSig.empty() && simplifiedSig.size() < sig.size()) {
            ++countSimplified;
            if (newUnidentifiedCounts.contains(sig)) {
                auto nh = newUnidentifiedCounts.extract(sig);
                newUnidentifiedCounts[simplifiedSig] += nh.mapped();
            }
            newUnidentifiedCounts[simplifiedSig]+=pair.second;
            hardComponentsMap[sig]=simplifiedSig;
            auto nh = gridMap.extract(sig);
            if (gridMap.contains(simplifiedSig)) {
                gridMap[simplifiedSig].insert(gridMap[simplifiedSig].end(),nh.mapped().begin(),nh.mapped().end());
            } else {
                gridMap[simplifiedSig] = nh.mapped();
            }
            std::cout << "S" << std::flush;
        } else {
            newToSimplify[sig]+=pair.second;
            std::cout << "F" << std::flush;
        }
    }

    std::cout << "\n\n";

    toSimplify = newToSimplify;
    newToSimplify.clear();

    std::cout << "Simplified " << countSimplified << " triangulations\n\n";

    std::cout << "elapsed: " << secondsSince(start) << " s\n\n\n";


    std::cout << "--- (2.3) MCMC simplification, third pass ---\n\n";

    countSimplified = 0;

    for (auto pair : toSimplify) {
        std::string sig = pair.first;
        std::vector<std::future<std::string>> workers;
        for (int i=0; i<cores; ++i) {
            workers.emplace_back(std::async(mcmc::simplify,sig,10000,100,false));
        }
        std::string simplifiedSig = workers[0].get();
        for (int i=1; i<cores; ++i) {
            std::string newsig = workers[i].get();
            if (!newsig.empty() &&newsig.size() < simplifiedSig.size()) {
                simplifiedSig = newsig;
            }
            if (newsig.empty()) {std::cout << "!";}
        }
        workers.clear();

        if (!simplifiedSig.empty() && simplifiedSig.size() < sig.size()) {
            ++countSimplified;
            if (newUnidentifiedCounts.contains(sig)) {
                auto nh = newUnidentifiedCounts.extract(sig);
                newUnidentifiedCounts[simplifiedSig] += nh.mapped();
            }
            newUnidentifiedCounts[simplifiedSig]+=pair.second;
            hardComponentsMap[sig]=simplifiedSig;
            auto nh = gridMap.extract(sig);
            if (gridMap.contains(simplifiedSig)) {
                gridMap[simplifiedSig].insert(gridMap[simplifiedSig].end(),nh.mapped().begin(),nh.mapped().end());
            } else {
                gridMap[simplifiedSig] = nh.mapped();
            }
            std::cout << "S" << std::flush;
        } else {
            newUnidentifiedCounts[sig]+=pair.second;
            std::cout << "F" << std::flush;
        }
    }

    std::cout << "\n\n";

    unidentifiedCounts = newUnidentifiedCounts;
    newUnidentifiedCounts.clear();

    std::cout << "Simplified " << countSimplified << " triangulations\n\n";

    std::cout << "elapsed: " << secondsSince(start) << " s\n\n\n";


    std::cout << "--- (2.4) census lookup --- \n\n";

    for (auto pair : unidentifiedCounts) {
        std::string sig = pair.first;
        regina::Triangulation<3> trig(sig);
        std::list<regina::CensusHit> lookup = regina::Census::lookup(trig);
        auto recog = regina::StandardTriangulation::recognise(trig);
        if (lookup.empty() && !recog) {
            newUnidentifiedCounts[sig]+=pair.second;
        } else {
            // our triangulation is in a census, so we assume it is minimal (or near-minimal)
            // and stop trying to simplify it
            // I believe Regina's censuses don't contain any connected sums, but just in case
            // we confirm it is irreducible here (and we don't need to decompose it)
            if (  (!lookup.empty() && lookup.front().name().starts_with("S2 x S1 : #")) |
                  (recog && recog->manifold() && recog->manifold()->name()=="S2 x S1") |
                  trig.isIrreducible() ) {
                finalCounts[sig]+=pair.second;
            } else {
                newUnidentifiedCounts[sig]+=pair.second;
            }
        }
    }

    unidentifiedCounts = newUnidentifiedCounts;
    newUnidentifiedCounts.clear();

    std::cout << "Identified near-minimal triangulations:\n\n";
    if (!finalCounts.empty()) {
        for (auto pair : finalCounts) {
            std::cout << pair.first << ": " << pair.second << "\n";
        }
        std::cout << "\n";
    }

    checkTotal = 0;
    for (auto pair : finalCounts) {
        checkTotal += pair.second;
    }
    std::cout << "Remaining: " << repeats-checkTotal << " / " << repeats << "\n\n";

    std::cout << "elapsed: " << secondsSince(start) << " s\n\n\n";


    std::cout << "--- (3.1) connected sums ---\n\n";

    producerDone = false;

    std::vector<std::thread> consumers;
    for (int i=0; i<cores; ++i) {
        consumers.emplace_back(std::thread(decomposeConsumer,i+1));
    }

    int i = 0;
    int total = unidentifiedCounts.size();

    // iterate in reverse order: maps sort by key size/lexicographic, so this (hopefully)
    // puts the long / time consuming ones earlier and lets threads which finish fast 
    // deal with the easy ones
    for (auto pair : unidentifiedCounts | std::views::reverse) {
        queue_sem.acquire();
        std::unique_lock<std::mutex> lck(queue_mtx);
        countQueue.push(pair);
        queue_cv.notify_all();
        ++i;
        if (i % 50 == 0) {
            std::scoped_lock<std::mutex> lck(output_mtx);
            std::cout << "\n***** Loaded " << i << "/" << total << " sigs to queue *****\n";
            std::cout << "elapsed: " << secondsSince(start) << " s\n\n";
        } 
    }

    {
        std::scoped_lock<std::mutex> lck(output_mtx);
        std::cout << "\n***** Loaded " << i << "/" << total << " sigs to queue *****\n\n";
        std::cout << "elapsed: " << secondsSince(start) << " s\n\n";
    } 

    {
        std::unique_lock<std::mutex> lck(queue_mtx);
        producerDone = true;
        queue_cv.notify_all();
    }

    for (auto& consumer : consumers) {
        consumer.join();
    }
    consumers.clear();

    std::cout << "Failed to decompose:\n\n";
    if (!couldntDecomposeCounts.empty()) {
        for (auto pair : couldntDecomposeCounts) {
            std::cout << pair.first << ": " << pair.second << "\n";
        }
        std::cout << "\n";
    }

    std::cout << "elapsed: " << secondsSince(start) << " s\n\n\n";


    std::cout << "--- (3.2) check triangulations we failed to decompose ---\n\n";

    for (auto pair : couldntDecomposeCounts) {

        std::cout << pair.first << "\n";
        regina::Triangulation<3> trig(pair.first);
        bool verifiedPrime = false;

        verifiedPrime = trig.isSphere();
        if (verifiedPrime) {
            std::cout << "\tis a sphere!" << " [elapsed: " << secondsSince(start) << "s]\n";
        } else {
            std::cout << "\tis not a sphere..." << " [elapsed: " << secondsSince(start) << "s]\n";
            verifiedPrime = trig.hasStrictAngleStructure();
            if (verifiedPrime) {
                std::cout << "\tis prime! (hyperbolic - strict angle structure)" << " [elapsed: " << secondsSince(start) << "s]\n";
            } else {
                std::cout << "\tdoesn't have a strict angle structure..." << " [elapsed: " << secondsSince(start) << "s]\n";;
                try {
                    verifiedPrime = (regina::SnapPeaTriangulation(trig).solutionType() == regina::SnapPeaTriangulation::Solution::Geometric);
                } catch (regina::SnapPeaFatalError &e) {} 
                  catch (regina::SnapPeaMemoryFull &e) {}
                if (verifiedPrime) {
                    std::cout << "\tis prime! (hyperbolic - SnapPy found geometric solution)" << " [elapsed: " << secondsSince(start) << "s]\n";
                } else {
                    std::cout << "\tSnapPy couldn't find a geometric solution..." << " [elapsed: " << secondsSince(start) << "s]\n";
                    regina::GroupPresentation gp = trig.group();
                    gp.simplify();
                    verifiedPrime = gp.identifyAbelian();
                    if (verifiedPrime) {
                        std::cout << "\tis prime! (fundamental group is abelian)" << " [elapsed: " << secondsSince(start) << "s]\n";
                    } else {
                        std::cout << "\tdoesn't have an (identifiably) abelian fundamental group..." << " [elapsed: " << secondsSince(start) << "s]\n";
                    }
                }
            }
        }

        // TODO - more/better check for hyperbolicity
        // we expect the problematic cases to be large prime manifolds, which are likely to be hyperbolic,
        // and verifying hyperbolicity *should* be much faster than running summands() / isIrreducible()
        // spawn a new process to have proper use of SnapPy / sage
        // randomise a few times and see if verify_hyperbolicity() gives true

        // TODO - outputs of recogniseGroup() that are definitively not free products

        if (verifiedPrime) {
            std::vector<std::string> decomp = {pair.first};
            decomposedCounts[decomp]+=pair.second;
        } else {
            std::cout << "    is not obviously prime, trying harder to simplify...\n";
            std::vector<std::future<std::string>> workers;
            for (int i=0; i<cores; ++i) {
                workers.emplace_back(std::async(mcmc::simplify,pair.first,100000,100,false));
            }
            std::string simplifiedSig = workers[0].get();
            for (int i=1; i<cores; ++i) {
                std::string newsig = workers[i].get();
                if (!newsig.empty() &&newsig.size() < simplifiedSig.size()) {
                    simplifiedSig = newsig;
                }
                if (newsig.empty()) {std::cout << "!";}
            }
            workers.clear();
            if (simplifiedSig.size() < pair.first.size()) {
                std::cout << " success!" << " [elapsed: " << secondsSince(start) << "s]\n";
                std::cout << "--> " << simplifiedSig << "\n";
                trig = regina::Triangulation<3>(simplifiedSig);
            } else {
                std::cout << "failed" << " [elapsed: " << secondsSince(start) << "s]\n";
            }

            std::cout << "trying to decompose again...\n";
            std::vector<regina::Triangulation<3>> summands = trig.summands();
            std::vector<std::string> decomp;
            if (summands.empty()) {
                // trig is S3
                decomp.push_back("bkaagj");
            } else {
                for (regina::Triangulation<3> summand : summands) {
                    summand.simplify();
                    decomp.push_back(summand.isoSig<regina::IsoSigDegrees<3,1>>());
                }
            }
            std::sort(decomp.begin(),decomp.end());

            std::string fullname;
            for (std::string sig : decomp) {
                fullname = fullname + sig + " ";
            }
            fullname.erase(fullname.size()-1,1); // erase last " "

            decomposedCounts[decomp]+=pair.second;
            auto nh = gridMap.extract(pair.first);
            if (gridMap.contains(fullname)) {
                gridMap[fullname].insert(gridMap[fullname].end(),nh.mapped().begin(),nh.mapped().end());
            } else {
                gridMap[fullname] = nh.mapped();
            }

            std::cout << "--> ";
            for (std::string s : decomp) {
                std::cout << s << " ";
            }
            std::cout << "" << " [elapsed: " << secondsSince(start) << "s]\n";
        }
    }
    std::cout << "\n";

    std::cout << "elapsed: " << secondsSince(start) << " s\n\n\n";


    std::cout << "--- (3.3) census lookup ---\n\n";

    std::map<std::vector<std::string>,int> newDecomposedCounts;
    std::set<std::string> hardComponents; 

    for (auto pair : decomposedCounts) {
        std::string fullname;
        bool identified = true;
        for (std::string sig : pair.first) {
            regina::Triangulation<3> trig(sig);
            std::list<regina::CensusHit> lookup = regina::Census::lookup(trig);
            auto recog = regina::StandardTriangulation::recognise(trig);
            if (lookup.empty() && !recog) {
                hardComponents.insert(sig);
                identified = false;
            } else {
                // note: this time we definitely don't need to check irreducibility,
                // since we've already decomposed into irreducible components
                fullname = fullname + sig + " ";
            }
        }
        if (identified) {
            fullname.erase(fullname.size()-1,1); // erase last " "
            finalCounts[fullname]+=pair.second;
        } else {
            newDecomposedCounts[pair.first]=pair.second;
        }
    }

    decomposedCounts = newDecomposedCounts;
    newDecomposedCounts.clear();

    std::cout << "Identified near-minimal triangulations / decompositions:\n\n";
    if (!finalCounts.empty()) {
        for (auto pair : finalCounts) {
            std::cout << pair.first << ": " << pair.second << "\n";
        }
        std::cout << "\n";
    }

    checkTotal = 0;
    for (auto pair : finalCounts) {
        checkTotal += pair.second;
    }
    std::cout << "Remaining: " << repeats-checkTotal << " / " << repeats << "\n\n";

    std::cout << "elapsed: " << secondsSince(start) << " s\n\n\n";
    

    // merge unidentified decompositions at last step into final counts
    for (auto pair : decomposedCounts) {
        std::string fullname;
        for (std::string sig : pair.first) {
            fullname = fullname + sig + " ";
        }
        fullname.erase(fullname.size()-1,1); // erase last " "
        finalCounts[fullname]+=pair.second;
    }

    std::cout << "Final counts:\n\n";
    for (auto pair : finalCounts) {
        std::cout << pair.first << ": " << pair.second << "\n";
    }
    std::cout << "\n";

    // sanity check
    checkTotal = 0;
    for (auto pair : finalCounts) {
        if (pair.second != gridMap[pair.first].size()) {
            std::cout << "Error: entry should have count of " << pair.second << " but has " << gridMap[pair.first].size() << " grids\n";
            std::cout << "(" << pair.first << ")\n";
        }
        checkTotal += pair.second;
    }
    if (checkTotal != repeats) {
        std::cout << "Error: total count (" << checkTotal << ") does not match the original number of repeats given (" << repeats << ")\n";
    }

    std::cout << "elapsed: " << secondsSince(start) << " s\n\n";

    if (saveResults) {
        std::string filename;
        if (knotOnly) {
            filename = "knots_n" + std::to_string(n) + "_" + std::to_string(repeats);
        } else {
            filename = "n" + std::to_string(n) + "_" + std::to_string(repeats);
        }
        saveToFile(finalCounts,filename);
    }

    return 0;
}