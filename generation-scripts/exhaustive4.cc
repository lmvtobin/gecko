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
 * Command line utility which exhaustively enumerates all 4-manifolds resulting from grids
 * of a given size (without 3-handles), calculates the intersection form for closed 4-manifolds
 * and simplifies the resulting ideal triangulations for non-closed 4-manifolds.
 * 
 * Results will be printed to console, and optionally saved to
 * ../data/exh4/nN_closedSC.txt and ../data/exh4/nN_ideal.txt
 * where N is the grid size
 * 
 * For closed (and simply connected) 4-manifolds, save format is one line per intersection form:
 * form.odd() form.rank() form.signature() : num_grids : grid_1 grid_2 ... grid_(num_grids)
 * the three values specifying the intesection form are all integers, to be interpreted as a
 * bool, size_t, and long respectively (i.e. a formTup)
 * 
 * For non-closed 4-manifolds, save format is one line per simplified isomorphism signature:
 * sig : num_grids : grid_1 grid_2 ... grid_(num_grids)
 * 
 * Note:
 * - in the closed case, lines represent distinct ~topological~ (not smooth) 4-manifolds
 * - in the non-closed case, lines represent smooth 4-manifolds, but are not guaranteed
 *   to all be distinct
 **/

#include <string>
#include <vector>
#include <queue>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <filesystem>

#include <thread>
#include <mutex>
#include <condition_variable>
#include <semaphore>

#include <triangulation/dim4.h>
#include <census/census.h>
#include <triangulation/isosigtype.h>
#include <algebra/intersectionform.h>

#include "../gecko/geckocore.h"

typedef std::tuple<bool,size_t,long> formTup;
typedef std::pair<gecko::Grid,int> gridWithIndex;

std::mutex queue_mtx;
std::condition_variable queue_cv;
std::queue<gecko::Grid> gridQueue;
bool producerDone = false;

std::mutex form_mtx;
std::map<formTup,int> formCounts;
std::map<formTup,std::vector<gridWithIndex>> formGridMap;

// note: without 3-handles, this map (of closed but not simply connected manifolds) was not necessary,
// as it should always be empty
// nonetheless it could be a good sanity check
// (if M is closed with no 3-handles, turn it upside down to get a decomp with no 1-handles
// now we are taking 0+2-handles, which is SC, and adding 3-handles, which can't produce more pi_1)

std::mutex ideal_mtx;
std::map<std::string,int> idealCounts;
std::map<std::string,std::vector<gridWithIndex>> gridMap;

bool closedOnly = false;

void gridConsumer() {

    std::clog << "worker started\n";

    while (true) {
        gecko::Grid grid;

        // get a grid from the queue

        {
            // wait for gridQueue to be nonempty and get the lock
            // once producer has loaded all grids, continue until gridQueue is empty
            std::unique_lock<std::mutex> lck(queue_mtx);
            queue_cv.wait(lck, []() { return !gridQueue.empty() or producerDone; });

            if (gridQueue.empty()) {
                // all grids generated and queue empty: done!
                break;
            } else {
                // get grid
                grid = gridQueue.front();
                gridQueue.pop();
            }
        }

        // process grid

        std::vector<regina::Triangulation<4>> trigs = gecko::buildAll4Manifolds(grid,false,true); // orient and don't simplify
        std::cout << "|" << std::flush;
        std::vector<std::pair<formTup,int>> forms;
        // std::vector<std::pair<std::string,int>> closedNSCSigs;
        std::vector<std::pair<std::string,int>> idealSigs;

        for (int i=0; i<trigs.size(); ++i) {
            regina::Triangulation<4> trig = trigs[i];
            // previously just did the obvious thing and called trig.isClosed()
            // now being extra careful in case regina ever misidentifies something with
            // hard to simplify spherical vertex links as ideal
            bool closed = true;
            for (auto v : trig.vertices()) {
                regina::Triangulation<3> link = v->buildLink();
                link.simplify();
                if (!link.isSphere()) {
                    closed = false;
                    break;
                }
            }
            if (closed) {
                regina::IntersectionForm form = trig.intersectionForm();
                formTup ft = {form.odd(),form.rank(),form.signature()};
                forms.push_back({ft,i});
            } else if (!closedOnly) {
                trig.simplify();
                std::string sig = trig.isoSig<regina::IsoSigDegrees<4,2>>();
                idealSigs.push_back({sig,i});
            }

        }

        {
            std::scoped_lock<std::mutex> lck(form_mtx);
            for (auto [ft,i] : forms) {
                ++formCounts[ft];
                gridWithIndex gwi = {grid,i};
                if (formGridMap.contains(ft)) {
                    formGridMap[ft].push_back(gwi);
                } else {
                    formGridMap[ft] = {gwi};
                }
            }
        }

        {
            std::scoped_lock<std::mutex> lck(ideal_mtx);
            for (auto [sig,i] : idealSigs) {
                ++idealCounts[sig];
                gridWithIndex gwi = {grid,i};
                if (gridMap.contains(sig)) {
                    gridMap[sig].push_back(gwi);
                } else {
                    gridMap[sig] = {gwi};
                }
            }
        }

    }

    std::clog << "\nworker done\n";

}

void saveToFile(std::map<std::string,int> counts, std::string filename)
{
    std::filesystem::create_directories("../data/exh4_data"); // create storage/output directory (if it doesn't exist)
    std::ofstream writer;
    writer.open("../data/exh4_data/"+filename+".txt");
    for (auto pair : counts)
    {
        writer << pair.first << " : " << pair.second << " : ";
        for (auto [grid,i] : gridMap[pair.first]) {
            gecko::printVector(grid.x,writer,false);
            gecko::printVector(grid.o,writer,false);
            writer << i;
            writer << " ";
        }
        writer << "\n";
    }
    writer.close();

    std::cout << "Saved to ../data/exh4_data/" << filename << ".txt\n";
}

void saveToFile(std::map<formTup,int> counts, std::string filename)
{
    std::filesystem::create_directories("../exh4_data"); // create storage/output directory (if it doesn't exist)
    std::ofstream writer;
    writer.open("../exh4_data/"+filename+".txt");
    for (auto pair : counts)
    {
        writer << std::get<0>(pair.first) << " " << std::get<1>(pair.first) << " " << std::get<2>(pair.first) << " : " << pair.second << " : ";
        for (auto [grid,i] : formGridMap[pair.first]) {
            gecko::printVector(grid.x,writer,false);
            gecko::printVector(grid.o,writer,false);
            writer << i;
            writer << " ";
        }
        writer << "\n";
    }
    writer.close();

    std::cout << "Saved to ../exh4_data/" << filename << ".txt\n";
}

int main(int argc, char* argv[]) {
    
    if (argc < 2) {
        std::cerr << "Usage:\n ./exhaustive4 gridsize [cores=1] [-save] [-closed]\n";
        std::cerr << "Example usage for grid size 6 with default values and not saving:\n";
        std::cerr << "./exhaustive4 6\n";
        std::cerr << "Example usage for grid size 6, on 10 cores, saving to default location, closed manifolds only:\n";
        std::cerr << "./exhaustive4 6 10 -save -closed\n";
        std::cerr << "Optional positional arguments must be in order, flags (-save, closed) must appear after including all optional arguments\n";
        return 1;
    }

    auto start = std::chrono::high_resolution_clock::now();

    int n = std::stoi(argv[1]);
    int cores = 1;
    bool saveResults = false;
    if (argc > 2) {
        cores = std::stoi(argv[2]);
    }
    for (int i=3; i<argc; ++i) {
        if (std::string(argv[i])=="-save") {
            saveResults = true;
        } else if (std::string(argv[i])=="-closed") {
            closedOnly = true;
        }
    }

    if (closedOnly) {
        std::cout << "Collecting forms...\n\n";
    } else {
        std::cout << "Collecting forms (SC) and iso sigs (non-SC)...\n\n";
    }

    std::vector<std::thread> consumers;
    for (int i=0; i<cores; ++i) {
        consumers.emplace_back(std::thread(gridConsumer));
    }

    gecko::Grid grid = gecko::firstGrid(n);
    do {
        std::unique_lock<std::mutex> lck(queue_mtx);
        gridQueue.push(grid);
        queue_cv.notify_all();
    } while (gecko::nextGrid(&grid));

    {
        producerDone = true;
        queue_cv.notify_all();
    }

    for (auto& consumer : consumers) {
        consumer.join();
    }
    consumers.clear();

    auto now = std::chrono::high_resolution_clock::now();
    std::cout << "elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(now-start).count() << " s\n\n";

    std::cout << "Closed (and simply connected):\n\n";
    if (!formCounts.empty()) {
        for (auto pair : formCounts) {
            std::cout << (std::get<0>(pair.first)? "odd" : "even") << " rank" << std::get<1>(pair.first) << " sig" << std::get<2>(pair.first) << " : " << pair.second << "\n";
        }
        std::cout << "\n";
    }

    if (!closedOnly) {
        std::cout << "Ideal:\n\n";
        if (!idealCounts.empty()) {
            for (auto pair : idealCounts) {
                std::cout << pair.first << ": " << pair.second << "\n";
            }
            std::cout << "\n";
        }
    }
    now = std::chrono::high_resolution_clock::now();
    std::cout << "elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(now-start).count() << " s\n\n";

    if (saveResults) {
        saveToFile(formCounts,"n"+std::to_string(n)+"_closedSC");
        if (!closedOnly) {
            saveToFile(idealCounts,"n"+std::to_string(n)+"_ideal");
        }
    }
    return 0;

}