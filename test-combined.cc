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
 * Runs test cases for the overall grid -> 3-manifold and grid -> 4-manifold procedures
 **/

#include <string>
#include <vector>
#include <chrono>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <snappea/snappeatriangulation.h>
#include <algebra/intersectionform.h>

#include "gecko/geckocore.h"
#include "katie/katie.h"

using namespace gecko;

typedef std::tuple<bool,size_t,long> formTup; // form.odd(), form.rank(), form.signature()

int main() {

    bool verbose = false;
    
    std::vector<std::pair<Grid,std::string>> testsDim3 = 
    {
        {Grid(),"S3"},
        {Grid({0},{0}), "S3"},
        {Grid(inverse({0,1}),inverse({1,0})),"S2 x S1"},
        {Grid(inverse({0,2,1}),inverse({1,0,2})),"S3"},
        {Grid(inverse({2,0,1}),inverse({1,2,0})),"S3"},
        {Grid(inverse({0,1,2,3}),inverse({1,0,3,2})),"S2 x S1 # S2 x S1"},
        {Grid(inverse({0,2,1,3}),inverse({1,0,3,2})),"RP3"},
        {Grid(inverse({2,3,0,1}),inverse({3,1,2,0})),"RP3"},
        {Grid(inverse({0,1,3,2}),inverse({3,2,0,1})),"S2 x S1 # S2 x S1"},
        {Grid(inverse({1,0,2,3}),inverse({2,3,1,0})),"S2 x S1 # S2 x S1"},
        {Grid(inverse({1,0,2,3}),inverse({3,2,0,1})),"S2 x S1 # S2 x S1"},
        {Grid(inverse({0,1,3,2}),inverse({2,3,1,0})),"S2 x S1 # S2 x S1"},
        {Grid(inverse({0,1,2,3}),inverse({2,3,0,1})),"S3"},
        {Grid(inverse({1,0,3,2}),inverse({3,2,1,0})),"S3"},
        {Grid(inverse({0,4,2,1,3}),inverse({1,0,3,4,2})),"S2 x S1"},
        {Grid(inverse({2,4,3,0,1}),inverse({3,1,2,4,0})),"S2 x S1"},
        {Grid(inverse({1,4,0,3,2}),inverse({2,1,3,0,4})),"S2 x S1"},
        {Grid(inverse({3,2,4,0,1}),inverse({1,4,2,3,0})),"S2 x S1"},
        {Grid(inverse({0,4,1,2,3}),inverse({2,0,3,4,1})),"S3"},
        {Grid(inverse({0,1,2,4,3}),inverse({2,3,0,1,4})),"S3"},
        {Grid(inverse({0,1,2,3,4,5}),inverse({1,0,3,2,5,4})),"S2 x S1 # S2 x S1 # S2 x S1"},
        {Grid(inverse({0,1,2,3,4,5,6,7,8,9,10,11}),inverse({1,0,3,2,5,4,7,6,9,8,11,10})),"S2 x S1 # S2 x S1 # S2 x S1 # S2 x S1 # S2 x S1 # S2 x S1"},
        {Grid(inverse({0,1,2,5,4,3}),inverse({5,4,3,0,1,2})),"S2 x S1 # S2 x S1 # S2 x S1"},
        {Grid(inverse({1,0,3,2,5,4}),inverse({3,5,1,4,0,2})),"S2 x S1 # S2 x S1 # S2 x S1"},
        {Grid(inverse({0,2,1,3,5,4}),inverse({1,0,2,4,3,5})),"S3"},
        {Grid(inverse({0,2,1,5,3,4}),inverse({1,0,2,4,5,3})),"S3"},
        {Grid(inverse({2,0,1,3,5,4}),inverse({1,2,0,4,3,5})),"S3"},
        {Grid(inverse({2,0,1,5,3,4}),inverse({1,2,0,4,5,3})),"S3"},
        {Grid(inverse({0,5,2,1,4,3}),inverse({1,0,3,5,2,4})),"S3"},
        {Grid(inverse({0,5,4,1,2,3}),inverse({1,0,3,5,4,2})),"S3"},
        {Grid(inverse({0,3,5,2,4,1}),inverse({1,0,4,5,2,3})),"S3"},
        {Grid(inverse({0,2,1,3,5,4,6,8,7}),inverse({1,0,2,4,3,5,7,6,8})),"S3"},
        {Grid(inverse({2,0,1,3,5,4,6,8,7}),inverse({1,2,0,4,3,5,7,6,8})),"S3"},
        {Grid(inverse({0,2,1,5,3,4,6,8,7}),inverse({1,0,2,4,5,3,7,6,8})),"S3"},
        {Grid(inverse({0,2,1,3,5,4,8,6,7}),inverse({1,0,2,4,3,5,7,8,6})),"S3"},
        {Grid(inverse({0,2,1,3,5,4,8,6,7,9,11,10}),inverse({1,0,2,4,3,5,7,8,6,10,9,11})),"S3"},
        {Grid(inverse({0,4,2,1,8,5,3,7,6}),inverse({1,0,3,4,2,6,8,5,7})),"S3"},
        {Grid(inverse({0,4,2,1,8,7,3,5,6}),inverse({1,0,3,4,2,6,8,7,5})),"S3"},
        {Grid(inverse({0,4,2,1,5,3,6,8,7}),inverse({1,0,3,4,2,5,7,6,8})),"S3"},
        {Grid(inverse({0,4,2,1,5,3,8,6,7}),inverse({1,0,3,4,2,5,7,8,6})),"S3"},
        {Grid(inverse({0,5,4,1,2,3,6,8,7}),inverse({1,0,3,5,4,2,7,6,8})),"S3"},
        {Grid(inverse({6,7,8,0,1,2,3,4,5}),inverse({8,1,6,2,3,4,5,7,0})),"S3"},
        {Grid({1,2,3,5,0,4,6},{4,0,1,2,3,6,5}),"S3/P120"},
        {Grid({0,1,2,3,4,6,5},{2,4,0,5,1,3,6}), "S3"}, // link where highlighting runs into another quadricolour
        {Grid({0,1,3,2,4,5},{4,2,1,5,0,3}), "S3"}, // link which results in an extra "wrong" graph quadricolour
        {Grid({0,1,5,3,2,4,6},{4,2,1,6,5,0,3}),"S3"}, // link which broke things when highlighted undercrossings were added before over
        {Grid({0,1,5,6,2,4,3},{4,2,1,3,5,0,6}),"S3"}, // link same as above with some x/o's swapped, which also broke things with the side of the region for 1-handles
        {Grid({0,13,6,7,9,4,3,2,5,11,8,14,10,12,1},{1,0,12,13,6,7,5,4,3,2,10,9,8,14,11}),"S3"}, // link where we shouldn't really be able to always add highlighted overcrossings first
        {Grid({0,2,4,1,3},{1,0,3,4,2}),"S3"}, // next 8: contain "double twists" which result in an extra wrong graph quadricolour
        {Grid({0,3,1,4,2},{1,0,4,2,3}),"S3"},
        {Grid({0,3,2,4,1},{2,0,4,1,3}),"S3"},
        {Grid({0,4,2,1,3},{1,3,0,4,2}),"S3"},
        {Grid({0,1,3,5,2,4},{1,2,0,4,5,3}), "S3"},
        {Grid({0,4,2,5,1,3},{4,1,0,3,5,2}), "S3"},
        {Grid({3,0,2,1,5,4},{5,1,0,4,3,2}), "S3"},
        {Grid({0,2,1,5,4,3},{1,0,4,3,2,5}), "S3"}
    };

    std::vector<std::tuple<Grid,std::string,std::string>> failsDim3;

    std::vector<std::pair<Grid,std::vector<formTup>>> testsDim4SC = 
    {
        {Grid(),{{false,0,0}}},
        {Grid({0},{0}), {{false,0,0}}},
        {Grid(inverse({0,2,1}),inverse({1,0,2})),{{true,1,1}}},
        {Grid(inverse({2,0,1}),inverse({1,2,0})),{{true,1,-1}}},
        {Grid(inverse({0,1,2,3}),inverse({2,3,0,1})),{{false,2,0},{false,0,0},{false,0,0}}},
        {Grid(inverse({1,0,3,2}),inverse({3,2,1,0})),{{false,2,0},{false,0,0},{false,0,0}}},
        {Grid(inverse({0,4,1,2,3}),inverse({2,0,3,4,1})),{{true,2,0},{false,0,0}}},
        {Grid(inverse({0,1,2,4,3}),inverse({2,3,0,1,4})),{{true,2,0},{false,0,0}}},
        {Grid(inverse({0,2,1,3,5,4}),inverse({1,0,2,4,3,5})),{{true,2,2}}},
        {Grid(inverse({0,2,1,5,3,4}),inverse({1,0,2,4,5,3})),{{true,2,0}}},
        {Grid(inverse({2,0,1,3,5,4}),inverse({1,2,0,4,3,5})),{{true,2,0}}},
        {Grid(inverse({2,0,1,5,3,4}),inverse({1,2,0,4,5,3})),{{true,2,-2}}},
        {Grid(inverse({0,5,2,1,4,3}),inverse({1,0,3,5,2,4})),{{true,2,2}}},
        {Grid(inverse({0,5,4,1,2,3}),inverse({1,0,3,5,4,2})),{{true,2,0}}},
        {Grid(inverse({0,3,5,2,4,1}),inverse({1,0,4,5,2,3})),{{true,2,0}}},
        {Grid(inverse({0,2,1,3,5,4,6,8,7}),inverse({1,0,2,4,3,5,7,6,8})),{{true,3,3}}},
        {Grid(inverse({2,0,1,3,5,4,6,8,7}),inverse({1,2,0,4,3,5,7,6,8})),{{true,3,1}}},
        {Grid(inverse({0,2,1,5,3,4,6,8,7}),inverse({1,0,2,4,5,3,7,6,8})),{{true,3,1}}},
        {Grid(inverse({0,2,1,3,5,4,8,6,7}),inverse({1,0,2,4,3,5,7,8,6})),{{true,3,1}}},
        {Grid(inverse({0,2,1,3,5,4,8,6,7,9,11,10}),inverse({1,0,2,4,3,5,7,8,6,10,9,11})),{{true,4,2}}},
        {Grid(inverse({0,4,2,1,8,5,3,7,6}),inverse({1,0,3,4,2,6,8,5,7})),{{true,3,3}}},
        {Grid(inverse({0,4,2,1,8,7,3,5,6}),inverse({1,0,3,4,2,6,8,7,5})),{{true,3,1}}},
        {Grid(inverse({0,4,2,1,5,3,6,8,7}),inverse({1,0,3,4,2,5,7,6,8})),{{true,3,3}}},
        {Grid(inverse({0,4,2,1,5,3,8,6,7}),inverse({1,0,3,4,2,5,7,8,6})),{{true,3,1}}},
        {Grid(inverse({0,5,4,1,2,3,6,8,7}),inverse({1,0,3,5,4,2,7,6,8})),{{true,3,1}}},
        {Grid(inverse({6,7,8,0,1,2,3,4,5}),inverse({8,1,6,2,3,4,5,7,0})),{{true,2,0},{false,0,0}}},
        {Grid({0,1,2,3,4,6,5},{2,4,0,5,1,3,6}), {{true,3,1},{true,1,1}}}, // link where highlighting runs into another quadricolour
        {Grid({0,1,3,2,4,5},{4,2,1,5,0,3}), {{false,2,0},{false,0,0}}}, // link which results in an extra "wrong" graph quadricolour
        {Grid({0,1,5,3,2,4,6},{4,2,1,6,5,0,3}),{{true,3,1},{true,1,1},{true,1,1}}}, // link which broke things when highlighted undercrossings were added before over
        {Grid({0,1,5,6,2,4,3},{4,2,1,3,5,0,6}),{{true,3,1},{true,1,1},{true,1,1}}}, // link same as above with some x/o's swapped, which also broke things with the side of the region for 1-handles
        {Grid({0,13,6,7,9,4,3,2,5,11,8,14,10,12,1},{1,0,12,13,6,7,5,4,3,2,10,9,8,14,11}),{{true,4,0},{true,2,0},{true,2,0},{false,0,0}}}, // link where we shouldn't really be able to always add highlighted overcrossings first
        {Grid({0,2,4,1,3},{1,0,3,4,2}), {{true,1,1}}}, // next 8: contain "double twists" which result in an extra wrong graph quadricolour
        {Grid({0,3,1,4,2},{1,0,4,2,3}), {{true,1,1}}},
        {Grid({0,3,2,4,1},{2,0,4,1,3}), {{true,1,1}}},
        {Grid({0,4,2,1,3},{1,3,0,4,2}), {{true,1,1}}},
        {Grid({0,1,3,5,2,4},{1,2,0,4,5,3}), {{true,1,1}}},
        {Grid({0,4,2,5,1,3},{4,1,0,3,5,2}), {{true,1,1}}},
        {Grid({3,0,2,1,5,4},{5,1,0,4,3,2}), {{false,2,0},{false,0,0}}},
        {Grid({0,2,1,5,4,3},{1,0,4,3,2,5}), {{false,2,0},{false,0,0}}}
    };

    std::vector<std::tuple<Grid,formTup,int,std::string>> failsDim4SC;

    std::cout << "------------\nDimension 3:\n------------\n";

    int count = 0;
    for (auto [grid,mfld] : testsDim3) {
        ++count;
        if (verbose) {
            std::cout << "\n";
            drawGrid(grid);
            printVector(grid.x);
            printVector(grid.o);
            std::cout << "\n";
        }
        regina::Triangulation<3> trig = build3Manifold(grid,false,true);
        if (verbose) {std::cout << trig.sig() << "\n\n";}
        if (!trig.isValid()) {
            failsDim3.push_back({grid,mfld,"Invalid"});
        } else if (!trig.isClosed()) {
            failsDim3.push_back({grid,mfld,"Not Closed"});
        } else if (!trig.isConnected()) {
            failsDim3.push_back({grid,mfld,"Disconnected"});
        } else if (!trig.isOriented()) {
            failsDim3.push_back({grid,mfld,"Not Oriented"});
        } else {
            trig.simplify();
            std::string id = identify3Manifold(trig.sig());
            if (verbose) {std::cout << "should be: " << mfld << "\n";}
            if (verbose) {std::cout << "found:     " << id << "\n";}
            if (id != mfld) {
                failsDim3.push_back({grid,mfld,id});
            } else {
                std::cout << "(" << count << "/" << testsDim3.size() << ") success\n";
                continue;
            }
        }
        std::cout << "(" << count << "/" << testsDim3.size() << ") failure <---\n";
    }

    std::cout << "\n\n------------\nDimension 4:\n------------\n";

    count = 0;
    for (auto [grid,tups] : testsDim4SC) {
        ++count;
        if (verbose) {
            std::cout << "\n";
            drawGrid(grid);
            printVector(grid.x);
            printVector(grid.o);
            std::cout << "\n";
        }
        std::vector<regina::Triangulation<4>> trigs = buildAll4Manifolds(grid,false,true);
        std::cout << "(" << count << "/" << testsDim4SC.size() << ")\n";
        for (int i=0; i<trigs.size(); ++i) {
            regina::Triangulation<4> trig = trigs[i];
            formTup tup = tups[i];
            if (verbose) {std::cout << "  " << trig.sig() << "\n\n";}
            if (!trig.isValid()) {
                failsDim4SC.push_back({grid,tup,i,"Invalid"});
            } else if (!trig.isClosed()) {
                failsDim4SC.push_back({grid,tup,i,"Not Closed"});
            } else if (!trig.isConnected()) {
                failsDim4SC.push_back({grid,tup,i,"Disconnected"});
            } else if (!trig.isOriented()) {
                failsDim4SC.push_back({grid,tup,i,"Not Oriented"});
            } else {
                regina::IntersectionForm form = trig.intersectionForm();
                if (verbose) {std::cout << "  should be: " << (std::get<0>(tup) ? "odd" : "even") << ", rank " << std::get<1>(tup) << ", signature " << std::get<2>(tup) << "\n";}
                if (verbose) {std::cout << "  found:     " << (form.odd() ? "odd" : "even") << ", rank " << form.rank() << ", signature " << form.signature() << "\n";}
                if (form.odd() != std::get<0>(tup) || form.rank() != std::get<1>(tup) || form.signature() != std::get<2>(tup)) {
                    failsDim4SC.push_back({grid,tup,i,form.str()});
                } else {
                    std::cout << "  [" << i+1 << "/" << trigs.size() << "] success\n";
                    continue;
                }
            }
            std::cout << "  [" << i+1 << "/" << trigs.size() << "] failure <---\n";
        }
    }

    if (failsDim3.empty()) {
        std::cout << "\nDimension 3: all tests successful!\n";
    } else {
        std::cout << "\nDimension 3: " << failsDim3.size() << " tests failed\n";
    }

    if (failsDim4SC.empty()) {
        std::cout << "\nDimension 4: all tests successful!\n";
    } else {
        std::cout << "\nDimension 4: " << failsDim4SC.size() << " tests failed\n";
    }

    for (auto [grid,mfld,fail] : failsDim3) {
        std::cout << "\n----- Failure -----\n\n";
        drawGrid(grid);
        printVector(grid.x);
        printVector(grid.o);
        std::cout << "\n" << fail << "\n";
        std::cout << "Should be: " << mfld << "\n";
    }

    for (auto [grid,tup,i,fail] : failsDim4SC) {
        std::cout << "\n----- Failure -----\n\n";
        drawGrid(grid);
        printVector(grid.x);
        printVector(grid.o);
        std::cout << "\n" << fail << "\n";
        std::cout << "Should be: " << (std::get<0>(tup) ? "odd" : "even") << ", rank " << std::get<1>(tup) << ", signature " << std::get<2>(tup) << "\n";
    }

    return 0;
}
