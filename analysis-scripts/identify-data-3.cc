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
 * Command line utility which reads the content of a data file as produced by sample3 or exhaustive3
 * and (attempts to) identify the 3-manifolds present.
 * 
 * Will print the number of grids corresponding to each 3-manifold, and one example grid for each.
 * Any manifolds that were unable to be identified are given as isomorphism signatures.
 **/

#include <string>
#include <vector>
#include <chrono>

#include <triangulation/dim3.h>
#include <snappea/snappeatriangulation.h>

#include "../gecko/geckocore.h"

using namespace gecko;

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage:\n ./identify-data-3 filepath [-knotOnly]\n";
        return 1;
    }

    // if true, filters for only knots in the sample
    // this is NOT necessary for a sample that already consists of only knots
    // (and will slightly slow things down unnecessarily in that case)
    bool knotOnly = false;
    for (int i=3; i<argc; ++i) {
        if (std::string(argv[i])=="-knotOnly") {
            knotOnly = true;
        }
    }

    std::vector<std::pair<std::vector<std::string>,std::vector<Grid>>> data = readSampleData(argv[1]);
    std::map<std::string,int> identified;
    std::map<std::string,Grid> gridMap;
    for (auto [sigs,grids] : data) {
        int count = 0;
        Grid first;
        if (knotOnly) {
            for (Grid grid : grids) {
                auto info = linkInfo(grid);
                if (info.empty() || info.size() == 1 && std::get<1>(info[0]).size() == 1) { // <-- knot only
                    // this is a knot
                    if (count == 0) {
                        first = grid;
                    }
                    ++count;
                }
            }
            if (count == 0) {
                continue;
            }
        } else {
            count = grids.size();
            first = grids[0];
        }
        std::string name = "";
        for (int i=0; i<sigs.size()-1; ++i) {
            std::string part = identify3Manifold(sigs[i]);
            if (part == "?") {
                part = sigs[i];
            }
            name = name + part + " # ";
        }
        std::string part = identify3Manifold(sigs[sigs.size()-1]);
        if (part == "?") {
            part = sigs[sigs.size()-1];
        }
        name = name + part;
        identified[name] += count;
        if (!gridMap.contains(name)) {
            gridMap[name] = first;
        }
    }

    std::cout << "Examples:\n\n";

    for (auto [name,count] : identified) {
        std::cout << name << "\n";
        printVector(gridMap[name].x,std::cout);
        printVector(gridMap[name].o,std::cout);
        drawGrid(gridMap[name]);
    }

    std::cout << "\nCounts:\n\n";

    int total = 0;
    for (auto [name,count] : identified) {
        std::cout << name << " : " << count << "\n";
        total  = total+count;
    }
    std::cout << "\nTotal : " << total << "\n";

    return 0;
}
