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
 * Command line utility which generates a sample of random 3-manifolds, and uses fast
 * heuristics to identify hyperbolic / non-hyperbolic manifolds and compute their
 * hyperbolic volumes
 * 
 * Results will be printed to console, and optionally saved to
 * ../sample3-data/volumes_nN_S_I.txt or ../sample3-data/volumes_knot_nN_S_I.txt
 * where
 * - N is the grid size
 * - S is the sample size
 * - I = 0,1,2,... is a counter which increments if a file with this N and S already exists 
 * 
 * Save format has a preamble with information about the distribution, followed first 
 * by samples identified to be likely hyperbolic, then samples identified to be likely
 * non-hyerbolic, with headers "\nlikely hyperbolic:\n" and "\nlikely non-hyperbolic:\n"
 * 
 * For each sample, there is a line
 * volume grid
 * where volume is a hex string for a floating point, as in std::hexfloat or Python's float.fromhex()
 * e.g. 
 * 0x1.cb17f914c9057p+3 [13,4,5,14,8,17,1,11,9,19,12,7,16,3,10,15,18,20,6,0,2][16,5,13,17,3,2,6,1,18,14,0,11,8,19,20,10,4,9,12,15,7]
 * 
 * Note as yet this is NOT parallelised: so far, the typical use case has been to sample at a 
 * variety of different grid sizes, so it makes more sense to just run multiple instances 
 * simultaneously
 **/

#include <fstream>
#include <filesystem>
#include <chrono>
#include <algorithm>

#include <vector>

#include <triangulation/isosigtype.h>
#include <snappea/snappeatriangulation.h>

#include "../gecko/geckocore.h"

int secondsSince(std::chrono::time_point<std::chrono::high_resolution_clock> start) {
    return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-start).count();
}

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage:\n ./sample3-volume gridsize [samples=1] [-save] [-knot] [-simplify]\n";
        std::cerr << "Example usage for grid size 6 with default values and not saving:\n./sample3-volume 6\n";
        std::cerr << "Example usage for grid size 6, 50 samples, saving to default location, only knots, simplifying triangulations:\n./sample3-volume 6 50 -save -knot -simplify\n";
        std::cerr << "Optional positional arguments must be in order, flags (-save, -knot, -simplify) must appear after including all optional arguments\n";
        return 1;
    }

    auto start = std::chrono::high_resolution_clock::now();

    int n = std::stoi(argv[1]);
    int repeats = 1;
    bool saveResults = false;
    bool knotOnly = false;
    bool simplify = false;
    if (argc > 2) {
        repeats = std::stoi(argv[2]);
    }
    for (int i=3; i<argc; ++i) {
        if (std::string(argv[i])=="-save") {
            saveResults = true;
        } else if (std::string(argv[i])=="-knot") {
            knotOnly = true;
        } else if (std::string(argv[i])=="-simplify") {
            simplify = true;
        }
    }

    std::vector<double> volumes_hyp;
    std::vector<gecko::Grid> grids_hyp;
    std::vector<double> volumes_nonhyp;
    std::vector<gecko::Grid> grids_nonhyp;

    std::vector<double> volumes;
    int nonhyp_count = 0;
    for (int i=0; i<repeats; ++i) {
        std::pair<regina::Triangulation<3>,gecko::Grid> triWithGrid = gecko::sample3Manifolds(n,1,knotOnly,simplify)[0];
        regina::Triangulation<3> tri = triWithGrid.first;
        gecko::Grid grid = triWithGrid.second;
        regina::SnapPeaTriangulation sptrig(tri);
        double sample_volume = 0;
        bool hyp = true;
        for (int i=0; i<100; ++i) {
            sptrig.randomise();
            double vol = sptrig.volume();
            if (vol < 0.5) {
                hyp = false;
                volumes_nonhyp.push_back(vol);
                grids_nonhyp.push_back(grid);
                break;
            } else {
                sample_volume += vol;
            }
        }
        if (hyp) {
            sample_volume = sample_volume / 100.0;
            if (sample_volume < 0.9) {
                hyp = false;
                volumes_nonhyp.push_back(sample_volume);
                grids_nonhyp.push_back(grid);
                break;
            }
            if (hyp) {
                volumes_hyp.push_back(sample_volume);
                grids_hyp.push_back(grid);
            }
        }
        std::cout << "|" << std::flush;
    }
    std::cout << "\n\n";

    double mean;
    double var;
    double median;

    if (volumes_hyp.empty()) {
        mean = 0;
        var = 0;
        median = 0;
    } else {

        mean = 0;
        for (double vol : volumes_hyp) {
            mean += vol;
        }
        mean = mean / volumes_hyp.size();

        var = 0;
        for (double vol : volumes_hyp) {
            var += (vol-mean)*(vol-mean);
        }
        var = var / volumes_hyp.size();

        std::vector<double> volumes_copy = volumes_hyp;
        const auto it_right = volumes_copy.begin() + volumes_copy.size() / 2;
        std::nth_element(volumes_copy.begin(), it_right, volumes_copy.end());
        if (volumes_copy.size() % 2 == 0) {
            const auto it_left = std::max_element(volumes_copy.begin(), it_right);
            median = (*it_left + *it_right) / 2;
        } else {
            median = *it_right;
        }

    }

    std::cout << n << " grid size, " << repeats << " samples\n";
    std::cout << volumes_hyp.size() << " likely hyperbolic ("<< 100.0*volumes_hyp.size()/repeats << "%)\n";
    std::cout << mean << " mean volume of likely hyperbolic samples\n";
    std::cout << median << " median volume of likely hyperbolic samples\n";
    std::cout << var << " variance among likely hyperbolic samples\n";

    std::cout << "  [took " << secondsSince(start) << "s]\n";

    if (saveResults) {
        std::filesystem::create_directories("../sample3_data"); // create storage/output directory (if it doesn't exist)
        std::string baseFilename = "";
        if (knotOnly) {
            baseFilename = "volumes_knot_n"+std::to_string(n)+"_"+std::to_string(repeats);
        } else {
            baseFilename = "volumes_n"+std::to_string(n)+"_"+std::to_string(repeats);
        }
        std::string filename = baseFilename;
        int filenum = 0;
        while (true) {
            filename = baseFilename + "_" + std::to_string(filenum);
            if (!std::filesystem::exists("../sample3_data/"+filename+".txt")) {
                break;
            } else {
                ++filenum;
            }
        }
        std::ofstream writer;
        writer.open("../sample3_data/"+filename+".txt");
        writer << n << " grid size, " << repeats << " samples\n";
        writer << volumes_hyp.size() << " likely hyperbolic ("<< 100.0*volumes_hyp.size()/repeats << "%)\n";
        writer << mean << " mean volume of likely hyperbolic samples\n";
        writer << median << " median volume of likely hyperbolic samples\n";
        writer << var << " variance among likely hyperbolic samples\n";
        writer << "\nlikely hyperbolic: \n\n";
        for (int i=0; i<volumes_hyp.size(); ++i) {
            writer << std::hexfloat << volumes_hyp[i] << " ";
            gecko::printVector(grids_hyp[i].x,writer,false);
            gecko::printVector(grids_hyp[i].o,writer,false);
            writer << "\n";
        }
        writer << "\nlikely non-hyperbolic: \n\n";
        for (int i=0; i<volumes_nonhyp.size(); ++i) {
            writer << std::hexfloat << volumes_nonhyp[i] << " ";
            gecko::printVector(grids_nonhyp[i].x,writer,false);
            gecko::printVector(grids_nonhyp[i].o,writer,false);
            writer << "\n";
        }
        writer.close();
        std::cout << "Saved to " + filename + "\n";
    }

    return 0;
}