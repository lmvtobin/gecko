//
// Katie
// Kirby Diagrams to Graphs and Triangulations
//
// Created by Rhuaidi Antonio Burke on 17/05/24.
// Copyright © 2024 Regina Development Team. All rights reserved.
//

/**
 * Edited by Lucy Tobin as part of GecKo.
 * 
 * Command line utility which recovers the original usage of Katie,
 * with the additional functionality added in Gecko.
 * 
 * In particular, katie-cmd now works for PD codes representing
 * disconnected diagrams, and should generally produce slightly
 * smaller triangulations, especially for dimension 3 and for
 * blackboard framed diagrams.
 **/

#include <iostream>
#include <vector>
#include <tuple>
#include <cctype>
#include <cstring>
#include <string>
#include <sstream>
#include <numeric>
#include <unistd.h>
#include <string>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <link/link.h>

#include "katie.h"

typedef std::vector<std::array<int, 4>> pdcode;

void usage(const char* progName, const std::string& error = std::string()) {
    if (!error.empty()) {
        std::cerr << error << "\n\n";
    }
    
    std::cerr << "Usage:" << std::endl;
    std::cerr << "    " << progName << " \"PD Code\" \"Framing Vector\", "
        " { -3, --dim3 | -4, --dim4 } "
        "[ -g, --graph ]\n" 
        // "[ -d, --debug ]\n"
        "    " << progName << " [ -v, --version | -?, --help ]\n\n";
    std::cerr << "    -3, --dim3    : Build a 3-manifold via integer "
        "Dehn surgery.\n";
    std::cerr << "    -4, --dim4    : Build a 4-manifold by attaching "
        "1- and 2-handles along a decorated link.\n";
    std::cerr << "                    The PD code must be the first argument and wrapped with quotation marks.\n";
    std::cerr << "                    The framing sequence must be the second argument and wrapped with quotation marks.\n";
    std::cerr << "                    Use 'x' or '.' to denote 1-handles within the framing sequence.\n\n";
    std::cerr << "    -g, --graph   : Output an edge-coloured graph, "
        "not an isomorphism signature.\n";
    //std::cerr << "    -r, --real    : Builds the 4-manifold triangulation with real boundary "
    //        "(not ideal or closed).\n";
    //std::cerr << "                    This option is incompatible with the --dim3 flag.\n\n";
    // std::cerr << "    -d, --debug   : Display debug information.\n";
    std::cerr << "    -v, --version : Show which version of Regina "
        "is being used\n";
    std::cerr << "    -?, --help    : Display this help\n\n";

    std::cerr << "Example usage:\n";
    std::cerr << progName << " \"PD: [(4,8,1,9),(9,3,10,4),(1,5,2,6),(6,2,7,3),(7,5,8,10)]\" \"x 0\"\n";
    
    exit(1);
}

int main(int argc, char* argv[]) {
    
    int dimFlag = 4; // Default to build a 4-manifold.
    bool outputGraph = false; // Default to ouptut an isomorphism signature.
    bool realBdry = false; // Default to build closed/ideal triangulation.
    
    // Check for standard arguments:
    for (int i=1; i<argc; ++i) {
        if (strcmp(argv[i], "-?") == 0 || strcmp(argv[i], "--help") == 0)
            usage(argv[0]);
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) {
            if (argc != 2)
                usage(argv[0],
                    "Option --version cannot be used with "
                        "any other arguments.");
            std::cout << PACKAGE_BUILD_STRING << std::endl;
            exit(0);
        }
    }
            
    /*
     START Process PD Code
     */
    std::string rawPDinput;
    
    if (argc < 3) {
        usage(argv[0], "Please provide a PD code and framing sequence.");
    }
    else {
        rawPDinput = argv[1];
    }
    
    /*
     "Sanitise" the raw input string:
     Blank everything that isn't a digit.
     Use sstream to handle 'multidigt' numbers.
     */
    for (char &c : rawPDinput) {
        if (!isdigit(c)) {
            c = ' ';
        }
    }
    
    std::stringstream ssPDC(rawPDinput);
    std::vector<int> rawPDVect;
    int currPDVal;
    while (ssPDC >> currPDVal) {
        rawPDVect.push_back(currPDVal);
    }
    /*
     END Process PD Code
     */
    
    /*
     START Process Framings
     */
    std::string rawFramingInput = argv[2];
    size_t rawFramingSize = rawFramingInput.size();
    
    std::vector<int> framingVector, twoHandleFramings;
    std::vector<bool> isOneHandleVector;
    
    std::istringstream iss(rawFramingInput);
    std::string framingToken;
    while (std::getline(iss, framingToken, ' ')) {
        if (framingToken == "x" || framingToken == ".") {
            framingVector.push_back(0);
            isOneHandleVector.push_back(true);
        } else {
            int framingInt = std::stoi(framingToken);
            framingVector.push_back(framingInt);
            twoHandleFramings.push_back(framingInt);
            isOneHandleVector.push_back(false);
        }
    }
    /*
     END Process Framings
     */

    if (3 <= argc && argc < 8) {
        for (int i=3; i<argc; ++i) {
            if (!strcmp(argv[i], "-3") || !strcmp(argv[i], "--dim3")) {
                dimFlag = 3;
            }
            else if (!strcmp(argv[i], "-4") || !strcmp(argv[i], "--dim4")) {
                dimFlag = 4;
            }
            else if (!strcmp(argv[i], "-g") || !strcmp(argv[i], "--graph")) {
                outputGraph = true;
            }
            else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--real")) {
                realBdry = true;
            }
            // else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--debug")) {
                // printDebugInfo = true;
            // }
            else {
                usage(argv[0], std::string("Invalid option: ") + argv[i]);
            }
        }
    }


    pdcode PDCworking;
        
    /*
     Check if the input PD code has come from the Snappy console.
     The Snappy console indexes the strands from 0, contrary to every other place.
     If the code has come from the Snappy console, bump everything up by 1.
     */
    bool codeFromSnappy = false;
    if (std::find(rawPDVect.begin(),rawPDVect.end(),0) != rawPDVect.end()) {
        codeFromSnappy = true;
    }

    if (codeFromSnappy) {
        for (int i=0; i<rawPDVect.size(); i++) {
            rawPDVect[i]++;
        }
    }
    
    for (int i=0; i<rawPDVect.size(); i+=4) {
        std::array<int, 4> currPDTup;
        for (int j=0; j<4; j++) {
            currPDTup[j] = rawPDVect[i+j];
        }
        PDCworking.push_back(currPDTup);
    }

    std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> info;
    regina::Link link = regina::Link::fromPD(PDCworking.begin(),PDCworking.end());
    std::vector<regina::Link> comps = link.diagramComponents();
    auto indicesPair = link.diagramComponentIndices();
    auto compIndices = indicesPair.first;
    size_t numNonTrivial = indicesPair.second;
    // std::set<int> remainingStrands;
    // for (int i=0; i<link.countComponents(); ++i) {
        // remainingStrands.insert(i);
    // }
    for (int i=0; i<numNonTrivial; ++i) {
        pdcode pdGrouped = comps[i].pdData();
        std::vector<int> pd;
        for (auto cross : pdGrouped) {
            for (int e : cross) {
                pd.push_back(e);
            }
        }
        std::vector<int> framing;
        std::vector<bool> oneHandle;
        for (int j=0; j<link.countComponents(); ++j) {
            auto strand = link.component(j);
            if (compIndices[strand.crossing()->index()] == i) {
                framing.push_back(framingVector[j]);
                oneHandle.push_back(isOneHandleVector[j]);
                // remainingStrands.erase(j);
            }
        }
        info.push_back({pd,framing,oneHandle});
    }
    // for (int i=numNonTrivial; i<comps.size(); ++i) {
    //     std::vector<int> pd = {};
    //     std::vector<int> framing;
    //     std::vector<bool> oneHandle;
    //     int j = *remainingStrands.begin();
    //     remainingStrands.erase(j);
    //     framing.push_back(framingVector[j]);
    //     oneHandle.push_back(isOneHandleVector[j]);
    //     info.push_back({pd,framing,oneHandle});
    // }

    if (dimFlag == 3) {
        // std::tuple<std::vector<int>,std::vector<int>> infoEl = {rawPDVect,framingVector};
        // std::vector<std::tuple<std::vector<int>,std::vector<int>>> info = {infoEl};
        std::vector<std::tuple<std::vector<int>,std::vector<int>>> infoTrunc;
        for (auto tup : info) {
            std::tuple<std::vector<int>,std::vector<int>> firstTwo = {std::get<0>(tup),std::get<1>(tup)};
            infoTrunc.push_back(firstTwo);
        }
        if (outputGraph) {
            std::clog << "\rHere is the edge list of the coloured graph:\n" << std::flush;
            katie::katie3PrintGraph(infoTrunc);    
        } else {
            regina::Triangulation<3> tri = katie::katie3(infoTrunc);
            std::clog << "\rHere is the isomorphism signature:\n" << std::flush;
            std::cout << tri.sig() << "\n";
        }
    } else if (dimFlag == 4) {
        // std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>> infoEl = {rawPDVect,framingVector,isOneHandleVector};
        // std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> info = {infoEl};
        if (outputGraph) {
            std::clog << "\rHere is the edge list of the coloured graph:\n" << std::flush;
            katie::katie4PrintGraph(info);    
        } else {
            regina::Triangulation<4> tri = katie::katie4(info);
            if (realBdry) {
                tri.truncateIdeal();
            }
            std::clog << "\rHere is the isomorphism signature:\n" << std::flush;
            std::cout << tri.sig() << "\n";
        }
    }

}