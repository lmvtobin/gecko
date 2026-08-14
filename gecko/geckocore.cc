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
 * Implementation for geckocore.h
 **/

#include <string>
#include <array>
#include <vector>
#include <set>
#include <map>
#include <stack>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <random>
#include <cmath>
#include <ranges>

#include <stdexcept>

#include <mutex>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example3.h>
#include <triangulation/example4.h>
#include <maths/spec/perm5.h>
#include <census/census.h>
#include <snappea/snappeatriangulation.h>
#include <surface/normalsurfaces.h>
#include <subcomplex/standardtri.h>

#include "geckocore.h"
#include "../katie/katie.h"

namespace gecko {

void printVector(std::vector<int> vec, std::ostream& stream, bool newline) {
    if (vec.empty()) {
        stream << "[]";
    } else if (vec.size()==1) {
        stream << "[" << *vec.begin() << "]";
    } else {
        stream << "[" << *vec.begin();
        for (auto iter = std::next(vec.begin()); iter!=vec.end(); ++iter) {
            stream << "," << *iter;
        }
        stream << "]";
    }
    if (newline) {stream << "\n";}
}

void printVector(std::vector<bool> vec, std::ostream& stream, bool newline) {
    if (vec.empty()) {
        stream << "[]";
    } else if (vec.size()==1) {
        stream << "[" << *vec.begin() << "]";
    } else {
        stream << "[" << *vec.begin();
        for (auto iter = std::next(vec.begin()); iter!=vec.end(); ++iter) {
            stream << "," << *iter;
        }
        stream << "]";
    }
    if (newline) {stream << "\n";}
}

std::vector<int> inverse(std::vector<int> perm) {
    std::vector<int> inv(perm.size(),0);
    for (int i=0; i<perm.size(); ++i) {
        inv[perm[i]]=i;
    }
    return inv;
}

std::vector<int> multiplyPerms(std::vector<int> p1, std::vector<int> p2) {
    std::vector<int> p(p1.size());
    for (int i=0; i<p1.size(); ++i) {
        p[i]=p1[p2[i]];
    }
    return p;
}

bool noCollisions(Grid grid) {
    for (int i=0; i<grid.x.size(); ++i) {
        if (grid.x[i]==grid.o[i]) {
            return false;
        }
    }
    return true;
}

std::vector<int> firstDerangement(int n) {
    std::vector<int> p(n);
    if (n%2 == 0) {
        for (int i=0; i<n/2; ++i) {
            p[2*i]=2*i+1;
            p[2*i+1]=2*i;
        }
    } else {
        for (int i=0; i<(n-3)/2; ++i) {
            p[2*i]=2*i+1;
            p[2*i+1]=2*i;
        }
        p[n-3]=n-2;
        p[n-2]=n-1;
        p[n-1]=n-3;
    }
    return p;
}

Grid firstGrid(int n) {
    std::vector<int> x(n);
    for (int i=0; i<n; ++i) {
        x[i]=i;
    }
    Grid grid = {x,firstDerangement(n)};
    return grid;
}

bool nextGrid(Grid* grid) {
    do {
        bool oDone = !std::next_permutation(grid->o.begin(),grid->o.end());
        if (oDone) {
            bool xDone = !std::next_permutation(grid->x.begin(),grid->x.end());
            if (xDone) {
                grid->o = firstDerangement(grid->x.size());
                return false;
            } else {
                grid->o = grid->x;
            }
        }
    } while (!noCollisions(*grid));
    return true;
}

std::tuple<std::vector<int>,std::vector<int>,std::vector<int>> coordToFollow(Grid grid) {
    std::vector<int> xrow = inverse(grid.x); // xrow[i] = row containing x in column i
    std::vector<int> rfollow;
    std::vector<int> cfollow;
    std::vector<int> partition;
    std::set<int> seen;
    int lowestRemaining = 1;
    int cycleStart = 0;
    int currentRow = 0;
    for (int i=0; i<xrow.size(); ++i) {
        seen.insert(currentRow);
        rfollow.push_back(currentRow);
        cfollow.push_back(grid.x[currentRow]);
        currentRow = xrow[grid.o[currentRow]];
        if (lowestRemaining == currentRow) {
            do {
                ++lowestRemaining;
            }
            while (seen.contains(lowestRemaining));
        }
        if (seen.contains(currentRow)) {
            currentRow = lowestRemaining;
            partition.push_back(i+1-cycleStart);
            cycleStart = i+1;
        } else {
            seen.insert(currentRow);
        }
    }
    return {rfollow,cfollow,partition};
}

Grid followToCoord(std::vector<int> rfollow, std::vector<int> cfollow, std::vector<int> partition) {
    std::vector<int> x(rfollow.size());
    std::vector<int> o(rfollow.size());
    int i=0;
    for (int length : partition) {
        for (int j=0; j<length-1; ++j) {
            x[rfollow[i]] = cfollow[i];
            o[rfollow[i]] = cfollow[i+1];
            ++i;
        }
        x[rfollow[i]] = cfollow[i];
        o[rfollow[i]] = cfollow[i+1-length];
        ++i;
    }
    return Grid(x,o);
}

Grid followToCoord(std::vector<int> rfollow, std::vector<int> cfollow) {
    int n = rfollow.size();
    std::vector<int> x(n);
    std::vector<int> o(n);
    for (int i=0; i<n-1; ++i) {
        x[rfollow[i]] = cfollow[i];
        o[rfollow[i]] = cfollow[i+1];
    }
    x[rfollow[n-1]] = cfollow[n-1];
    o[rfollow[n-1]] = cfollow[0];
    return Grid(x,o);
}

std::vector<std::set<int>> connectedComponents(std::map<int,std::set<int>> graph) {

    std::vector<std::set<int>> components;
    std::set<int> remainingVertices;
    for (auto pair : graph) {
        remainingVertices.insert(pair.first);
    }

    while (!remainingVertices.empty()) {

        std::stack<int> toCheck;
        toCheck.push(*remainingVertices.begin());
        std::set<int> component;
        component.insert(*remainingVertices.begin());
        remainingVertices.erase(remainingVertices.begin());

        while (!toCheck.empty()) {
            int vertex = toCheck.top();
            toCheck.pop();
            for (int neigh : graph[vertex]) {
                if (!component.contains(neigh)) {
                    component.insert(neigh);
                    toCheck.push(neigh);
                    remainingVertices.erase(neigh);
                }
            }
        }

        components.push_back(component);
    
    }

    return components;
}

std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<std::vector<bool>>>> linkInfo(Grid grid, bool noOnehandles, bool canonicalOneHandles, bool debug) {

    int n = grid.x.size();

    if (n < 2) {
        return std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<std::vector<bool>>>>();
    }

    std::vector<int> xrow = inverse(grid.x); // xrow[i] = row containing x in column i
    std::vector<int> ocol = grid.o; // ocol[i] = column containing o in row i

    /*
    First pass traversal along the link
    - find crossings
    - record direction of the second time encountering each crossing
      (so pd code entries can be ordered correctly on second pass)
    - count number of components
    - check for trivially disconnected graphs
    */

    if (debug) {std::clog << "First Pass";}

    std::set<int> xremaining; // note: sets are ordered
    for (int i=0; i<n; ++i) {
        xremaining.emplace(i);
    }
    std::set<std::pair<int,int>> seen;
    std::set<std::pair<int,int>> crossings;
    std::map<std::pair<int,int>,std::array<int,5>> crossingInfo; 
        // first 4 entries are PD code tuple
        // 5th is: 1 if the last time we saw this crossing it was left-to-right / bottom-to-top,
        //         0 if the last time was right-to-left / top-to-bottom
    std::map<int,int> strandToComponent;
    int component = 0;

    while (!xremaining.empty()) {

        if (debug) {std::clog << "\n\nComponent " << component;}

        auto first = xremaining.begin();
        int col = *first;
        xremaining.erase(first);
        std::pair<int,int> start {col,xrow[col]};
        std::pair<int,int> corner {col,xrow[col]};

        // remember first strand and keep track of most recent crossing,
        // to fix when we cycle back to the beginning
        std::pair<int,int> lastCrossing{-1,-1};

        while (true) {

            // go from x -> o horizontally

            std::pair<int,int> nextCorner {ocol[corner.second],corner.second};
            if (nextCorner.first > corner.first) {
                // going right
                if (debug) {std::clog << "\n(" << corner.first << "," << corner.second << ")X > -- ";}
                for (int i=corner.first+1; i<nextCorner.first; ++i) {
                    std::pair<int,int> step{i,corner.second};
                    if (debug) {std::clog << "(" << step.first << "," << step.second << ")";}
                    if (seen.contains(step)) {
                        if (debug) {std::clog << "C";}
                        crossings.insert(step);
                        crossingInfo[step] = std::array<int,5>{-1,-1,-1,-1,1};
                        lastCrossing = step;
                    } else{
                        seen.insert(step);
                    }
                    if (debug) {std::clog << " -- ";}
                }
            } else {
                // going left
                if (debug) {std::clog << "\n(" << corner.first << "," << corner.second << ")X < -- ";}
                for (int i=corner.first-1; i>nextCorner.first; --i) {
                    std::pair<int,int> step{i,corner.second};
                    if (debug) {std::clog << "(" << step.first << "," << step.second << ")";}
                    if (seen.contains(step)) {
                        if (debug) {std::clog << "C";}
                        crossings.insert(step);
                        crossingInfo[step] = std::array<int,5>{-1,-1,-1,-1,0};
                        lastCrossing = step;
                    } else{
                        seen.insert(step);
                    }
                    if (debug) {std::clog << " -- ";}
                } 
            }
            corner = nextCorner;

            // go from o -> x vertically

            nextCorner = {corner.first,xrow[corner.first]};
            if (nextCorner.second > corner.second) {
                // going up
                if (debug) {std::clog << "\n(" << corner.first << "," << corner.second << ")O ^ -- ";}
                for (int i=corner.second+1; i<nextCorner.second; ++i) {
                    std::pair<int,int> step{corner.first,i};
                    if (debug) {std::clog << "(" << step.first << "," << step.second << ")";}
                    if (seen.contains(step)) {
                        if (debug) {std::clog << "C";}
                        crossings.insert(step);
                        crossingInfo[step] = std::array<int,5>{-1,-1,-1,-1,1};
                        lastCrossing = step;
                    } else{
                        seen.insert(step);
                    }
                    if (debug) {std::clog << " -- ";}
                }
            } else {
                // going down
                if (debug) {std::clog << "\n(" << corner.first << "," << corner.second << ")O V -- ";}
                for (int i=corner.second-1; i>nextCorner.second; --i) {
                    std::pair<int,int> step{corner.first,i};
                    if (debug) {std::clog << "(" << step.first << "," << step.second << ")";}
                    if (seen.contains(step)) {
                        if (debug) {std::clog << "C";}
                        crossings.insert(step);
                        crossingInfo[step] = std::array<int,5>{-1,-1,-1,-1,0};
                        lastCrossing = step;
                    } else{
                        seen.insert(step);
                    }
                    if (debug) {std::clog << " -- ";}
                } 
            }
            corner = nextCorner;

            if (corner == start) {
                // we are back to the beginning of the component
                if (debug) {std::clog << "\n(" << corner.first << "," << corner.second << ")X";}
                ++component;
                break;
            } else {
                xremaining.erase(corner.first);
            }
        }
    }

    if (debug) {
        std::clog << "\n\nCrossings + direction of second occurence (1=right/up)\n";
        for (auto c : crossings) {
            std::clog << "(" << c.first << "," << c.second << "):" << crossingInfo[c][4];
            std::clog << " ";
        }
    }

    /*
    Second pass through
    - add pd code entries for first time crossing is seen, and put into correct order
    - record writhe of each component
    - check one-handle validity of each component, exluding linking between one-handles
      (zero writhe + split under- and over-crossings)
    */

    if (debug) {std::clog << "\n\nSecond Pass";}

    std::vector<int> writhes;
    std::vector<bool> validOnehandle; 

    xremaining.clear();
    for (int i=0; i<n; ++i) {
        xremaining.emplace(i);
    }
    int strand = 1;
    component = 0;

    while (!xremaining.empty()) {

        if (debug) {std::clog << "\n\nComponent " << component;}

        writhes.emplace_back(0);
        validOnehandle.emplace_back(true);

        auto first = xremaining.begin();
        int col = *first;
        xremaining.erase(first);
        std::pair<int,int> start{col,xrow[col]};
        std::pair<int,int> corner{col,xrow[col]};

        // remember first strand and keep track of most recent crossing,
        // to fix when we cycle back to the beginning
        int firstStrand = strand;
        std::pair<int,int> lastCrossing{-1,-1};

        int alternateCounter = 0; // number of times over/under alternates around the component
        int underOverTracker = 0; // 1 if last crossing was over, -1 if under, 0 if we are at the start

        std::set<std::pair<int,int>> thisStrandCrossings;

        while (true) {

            // go from x -> o horizontally

            std::pair<int,int> nextCorner {ocol[corner.second],corner.second};
            if (nextCorner.first > corner.first) {
                // going right
                if (debug) {std::clog << "\n(" << corner.first << "," << corner.second << ")X > -- ";}
                for (int i=corner.first+1; i<nextCorner.first; ++i) {
                    std::pair<int,int> step{i,corner.second};
                    if (debug) {std::clog << "(" << step.first << "," << step.second << ")";};
                    if (crossings.contains(step)) {
                        if (debug) {std::clog << "C";}
                        crossingInfo[step][0]=strand;
                        crossingInfo[step][2]=strand+1;
                        if (debug) {
                            std::clog << "{";
                            for (int j=0; j<=2; ++j) {
                                std::clog << crossingInfo[step][j] << ",";
                            }
                            std::clog << crossingInfo[step][3] << "}";
                        }
                        if (crossingInfo[step][1] == -1) {
                            // this is the first time seeing this crossing
                            thisStrandCrossings.insert(step);
                            crossingInfo[step][4] = 1;
                        } else if (thisStrandCrossings.contains(step)) {
                            // second time seeing, and it is a self-crossing
                            validOnehandle[component] = false;
                            if (crossingInfo[step][4] == 1) {
                                // over-strand goes bottom-to-top
                                if (debug) {std::clog << "[-w]";}
                                --writhes[component];
                            } else {
                                // over-strand goes top-to-bottom
                                if (debug) {std::clog << "[+w]";}
                                ++writhes[component];
                            }
                        }
                        if (underOverTracker == 1) {
                            if (debug) {std::clog << "[alt]";}
                            ++alternateCounter;
                        }
                        underOverTracker = -1;
                        ++strand;
                        lastCrossing = step;
                        seen.insert(step);
                    }
                    if (debug) {std::clog << " -- ";}
                }
            } else {
                // going left
                if (debug) {std::clog << "\n(" << corner.first << "," << corner.second << ")X < -- ";}
                for (int i=corner.first-1; i>nextCorner.first; --i) {
                    std::pair<int,int> step{i,corner.second};
                    if (debug) {std::clog << "(" << step.first << "," << step.second << ")";};
                    if (crossings.contains(step)) {
                        if (debug) {std::clog << "C";}
                        crossingInfo[step][0]=strand;
                        crossingInfo[step][2]=strand+1;
                        if (debug) {
                            std::clog << "{";
                            for (int j=0; j<=2; ++j) {
                                std::clog << crossingInfo[step][j] << ",";
                            }
                            std::clog << crossingInfo[step][3] << "}";
                        }
                        if (crossingInfo[step][1] == -1) {
                            // this is the first time seeing this crossing
                            thisStrandCrossings.insert(step);
                            crossingInfo[step][4] = 0;
                        } else if (thisStrandCrossings.contains(step)) {
                            // second time seeing, and it is a self-crossing
                            validOnehandle[component] = false;
                            if (crossingInfo[step][4] == 1) {
                                // over-strand goes bottom-to-top
                                if (debug) {std::clog << "[+w]";}
                                ++writhes[component];
                            } else {
                                // over-strand goes top-to-bottom
                                if (debug) {std::clog << "[-w]";}
                                --writhes[component];
                            }
                        }
                        if (underOverTracker == 1) {
                            if (debug) {std::clog << "[alt]";}
                            ++alternateCounter;
                        }
                        underOverTracker = -1;
                        ++strand;
                        lastCrossing = step;
                        seen.insert(step);
                    }
                    if (debug) {std::clog << " -- ";}
                } 
            }
            corner = nextCorner;

            // go from o -> x vertically

            nextCorner = {corner.first,xrow[corner.first]};
            if (nextCorner.second > corner.second) {
                // going up
                if (debug) {std::clog << "\n(" << corner.first << "," << corner.second << ")X ^ -- ";}
                for (int i=corner.second+1; i<nextCorner.second; ++i) {
                    std::pair<int,int> step{corner.first,i};
                    if (debug) {std::clog << "(" << step.first << "," << step.second << ")";};
                    if (crossings.contains(step)) {
                        if (debug) {std::clog << "C";}
                        if (crossingInfo[step][4] == 1) {
                            // under-strand goes left-to-right
                            crossingInfo[step][1]=strand;
                            crossingInfo[step][3]=strand+1;
                        } else {
                            // under-strand goes right-to-left
                            crossingInfo[step][1]=strand+1;
                            crossingInfo[step][3]=strand;
                        }
                        if (debug) {
                            std::clog << "{";
                            for (int j=0; j<=2; ++j) {
                                std::clog << crossingInfo[step][j] << ",";
                            }
                            std::clog << crossingInfo[step][3] << "}";
                        }
                        if (crossingInfo[step][0] == -1) {
                            // this is the first time seeing this crossing
                            thisStrandCrossings.insert(step);
                            crossingInfo[step][4] = 1;
                        } else if (thisStrandCrossings.contains(step)) {
                            // second time seeing, and it is a self-crossing
                            validOnehandle[component] = false;
                            if (crossingInfo[step][4] == 1) {
                                // under-strand goes left-to-right
                                if (debug) {std::clog << "[-w]";}
                                --writhes[component];
                            } else {
                                // under-strand goes right-to-left
                                if (debug) {std::clog << "[+w]";}
                                ++writhes[component];
                            }
                        }
                        if (underOverTracker == -1) {
                            if (debug) {std::clog << "[alt]";}
                            ++alternateCounter;
                        }
                        underOverTracker = 1;
                        ++strand;
                        lastCrossing = step;
                        seen.insert(step);
                    }
                    if (debug) {std::clog << " -- ";}
                }
            } else {
                // going down
                if (debug) {std::clog << "\n(" << corner.first << "," << corner.second << ")X ^ -- ";}
                for (int i=corner.second-1; i>nextCorner.second; --i) {
                    std::pair<int,int> step{corner.first,i};
                    if (debug) {std::clog << "(" << step.first << "," << step.second << ")";};
                    if (crossings.contains(step)) {
                        if (debug) {std::clog << "C";}
                        if (crossingInfo[step][4] == 1) {
                            // under-strand goes left-to-right
                            crossingInfo[step][1]=strand+1;
                            crossingInfo[step][3]=strand;
                        } else {
                            // under-strand goes right-to-left
                            crossingInfo[step][1]=strand;
                            crossingInfo[step][3]=strand+1;
                        }
                        if (debug) {
                            std::clog << "{";
                            for (int j=0; j<=2; ++j) {
                                std::clog << crossingInfo[step][j] << ",";
                            }
                            std::clog << crossingInfo[step][3] << "}";
                        }
                        if (crossingInfo[step][0] == -1) {
                            // this is the first time seeing this crossing
                            thisStrandCrossings.insert(step);
                            crossingInfo[step][4] = 0;
                        } else if (thisStrandCrossings.contains(step)) {
                            // second time seeing, and it is a self-crossing
                            validOnehandle[component] = false;
                            if (crossingInfo[step][4] == 1) {
                                // under-strand goes left-to-right
                                if (debug) {std::clog << "[+w]";}
                                ++writhes[component];
                            } else {
                                // under-strand goes right-to-left
                                if (debug) {std::clog << "[-w]";}
                                --writhes[component];
                            }
                        }
                        if (underOverTracker == -1) {
                            if (debug) {std::clog << "[alt]";}
                            ++alternateCounter;
                        }
                        underOverTracker = 1;
                        ++strand;
                        lastCrossing = step;
                        seen.insert(step);
                    }
                    if (debug) {std::clog << " -- ";}
                }
            }
            corner = nextCorner;

            if (corner == start) {
                // we are back to the beginning of the component
                if (debug) {std::clog << "\n(" << corner.first << "," << corner.second << ")X";}
                if (lastCrossing.first == -1) {
                    // this component has no crossings (an isolated circle)
                    if (debug) {std::clog << "\ncomponent is an isolated circle\n";}
                    if (debug) {std::clog << "component consists of strand " << firstStrand << "\n";}
                    strandToComponent[firstStrand]=component;
                    strand++;
                } else {
                    // the current value of strand is firstStrand + number of strands
                    // the last crossing we saw used this for the pd code when it
                    // should have used firstStrand since we've returned to the beginning
                    // fix it
                    for (int i=0; i<=5; ++i) {
                        if (crossingInfo[lastCrossing][i] == strand) {
                            crossingInfo[lastCrossing][i] = firstStrand;
                        }
                    }
                    if (debug) {
                        std::clog << "\nfix last crossing: ";
                        std::clog << "(" << lastCrossing.first << "," << lastCrossing.second << ") {";
                        for (int j=0; j<=2; ++j) {
                            std::clog << crossingInfo[lastCrossing][j] << ",";
                        }
                        std::clog << crossingInfo[lastCrossing][3] << "}\n";
                    }
                    if (debug) {std::clog << "component consists of strands ";}
                    for (int i=firstStrand; i<strand; ++i) {
                        if (debug) {std::clog << i << " ";}
                        strandToComponent[i] = component;
                    }
                    // valid one-handle must have 0 writhe, at least one each under and over-crossing
                    // and switch under/over exactly twice cyclically
                    // (alternateCounter=1 implies the second alternation is at the end of the cycle)
                    if (debug) {std::clog << "\nwrithe = " << writhes[component] << "\nalternated over/under " << alternateCounter << " times\n";}
                    if (writhes[component]!=0 || alternateCounter<1 || alternateCounter>2) {
                        if (debug) {std::clog << "\nnot a valid one-handle";}
                        validOnehandle[component] = false;
                    }
                }
                ++component;
                break;
            } else {
                xremaining.erase(corner.first);
            }
        }
    }

    int numComponents = component;

    if (debug) {
        std::clog << "\n\nWrithes: ";
        printVector(writhes,std::clog);
        std::clog << "Valid 1-handles: ";
        printVector(validOnehandle,std::clog);
    }

    /*
    Separate into connected components
    */

    if (debug) {std::clog << "\nBuilding connectivity graph on " << numComponents << " link components\n";}

    std::map<int,std::set<int>> graph; // vertices are link components, edge if they overlap
    std::vector<std::set<int>> graphComponents;

    if (numComponents > 1) {

        for (int i=0; i<numComponents; ++i) {
            graph[i] = {};
        }
        for (auto c : crossings) {
            int comp1 = strandToComponent[crossingInfo[c][0]];
            int comp2 = strandToComponent[crossingInfo[c][1]];
            if (comp1 != comp2) {
                graph[comp1].insert(comp2);
                graph[comp2].insert(comp1);
            }
        }

        if (debug) {
            for (int i=0; i<numComponents; ++i) {
                if (graph[i].empty()) {
                    std::clog << i << ":() ";
                } else if (graph[i].size()==1) {
                    std::clog << i << ":(" << *graph[i].begin() << ") ";
                } else {
                    std::clog << i << ":(" << *graph[i].begin();
                    for (auto iter = std::next(graph[i].begin()); iter!=graph[i].end(); ++iter) {
                        std::clog << "," << *iter;
                    }
                    std::clog << ") ";
                }
            }
            std::clog << "\n";
        }

        graphComponents = connectedComponents(graph);

    } else {

        graph[0] = std::set<int>();
        graphComponents.push_back({0});

    }

    if (debug) {std::clog << "\nGraph has " << graphComponents.size() << " components\n";}

    /*
    Check 1-handles do not cross each other
    - run through each valid 1-handle component in traversal order as above,
      and check all other 1-handles it crosses
    - if two one-handles cross, the one earlier in the order stays a 1-handle
      and the later one becomes a 0-framed 2-handle
    */

    std::vector<std::vector<bool>> oneHandleVecs;

    if (noOnehandles) {

        oneHandleVecs.emplace_back(std::vector<bool>(numComponents,false));

    } else {

        if (canonicalOneHandles) {

            if (debug) {std::clog << "\nChecking for linked 1-handles\n\n";}

            for (int i=0; i<numComponents; ++i) {
                if (validOnehandle[i]) {
                    if (debug) {std::clog << "component " << i << " is a 1-handle\n";}
                    for (int neigh : graph[i]) {
                        if (validOnehandle[neigh] && neigh>i) {
                            validOnehandle[neigh] = false;
                            if (debug) {std::clog << "component " << neigh << " is linked with 1-handle component " << i << "\n";}
                        }
                    }
                }
            }

            if (debug) {std::clog << "1-handles: "; printVector(validOnehandle,std::clog);}

            oneHandleVecs.push_back(validOnehandle);

        } else {

            if (debug) {std::clog << "\nFinding all combinations of unlinked 1-handles\n\n";}

            int numPossibleOneHandles = 0;
            for (bool valid : validOnehandle) {
                if (valid) {
                    ++numPossibleOneHandles;
                }
            }

            // build all possible combinations of potential 1-handles
            // (all possible boolean strings of length = numPossibleOneHandles, embedded
            //  at the appropriate positions in a vector of length numComponents)

            std::vector<std::vector<bool>> combos;

            for (int i=0; i<std::pow(2,numPossibleOneHandles); ++i) {
                combos.emplace_back(std::vector<bool>(numComponents,false));
            }

            int blocksize = 1;
            for (int i=numComponents-1; i>=0; --i) {
                if (validOnehandle[i]) {
                    bool b = false;
                    int counter = 0;
                    for (int j=0; j<combos.size(); ++j) {
                        combos[j][i]=b;
                        counter++;
                        if (counter==blocksize) {
                            b = !b;
                            counter = 0;
                        }
                    }
                    blocksize *= 2;
                }
            }
            
            // check which are valid

            for (std::vector<bool> combo : combos) {
                if (debug) {printVector(combo,std::clog);}
                bool valid = true;
                for (int i=0; i<combo.size(); ++i) {
                    if (combo[i]) {
                        for (int neigh : graph[i]) {
                            if (combo[neigh] && neigh>i) {
                                valid = false;
                                if (debug) {std::clog << "invalid, components " << i << " and " << neigh << " cross\n\n";}
                                break;
                            }
                        }
                        if (!valid) {
                            break;
                        }
                    }
                }
                if (valid) {
                    oneHandleVecs.push_back(combo);
                    if (debug) {std::clog << "valid\n\n";}
                }
            }


        }

    }

    /*
    Compile info by component
    */

    if (debug) {std::clog << "\nComponents:\n";}

    std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<std::vector<bool>>>> fullInfo;

    for (std::set<int> graphComp : graphComponents) {

        std::vector<int> compPD;
        std::vector<int> compWrithes;
        std::vector<std::vector<bool>> compOneHandleCombos(oneHandleVecs.size());

        for (auto c : crossings) {
            if (graphComp.contains(strandToComponent[crossingInfo[c][0]])) {
                for (int i=0; i<=3; ++i) {
                    compPD.push_back(crossingInfo[c][i]);
                }
            }
        }

        // pd code should have entries reduced so they are 0,1,2,...
        std::vector<int> sorted(compPD.begin(),compPD.end());
        std::sort(sorted.begin(),sorted.end());
        std::map<int,int> rescale;
        int i=0;
        int lastSeen=-1;
        for (int j : sorted) {
            if (j > lastSeen) {
                lastSeen = j;
                rescale[j]=i;
                ++i;
            }
        }

        for (int i=0; i<compPD.size(); ++i) {
            compPD[i]=rescale[compPD[i]];
        }

        for (int i : graphComp) {
            compWrithes.push_back(writhes[i]);
            for (int j=0; j<oneHandleVecs.size(); ++j) {
                compOneHandleCombos[j].push_back(oneHandleVecs[j][i]);
            }
        }

        fullInfo.push_back({compPD,compWrithes,compOneHandleCombos});

        if (debug) {
            std::clog << "\nPD: ";
            printVector(compPD,std::clog);
            std::clog << "Writhes: ";
            printVector(compWrithes,std::clog);
            std::clog << "Valid 1-handle combinations:\n";
            for (auto compOneHandle : compOneHandleCombos) {
                std::clog << "  ";
                printVector(compOneHandle,std::clog);
            }
            std::clog << "\n";
        }

    }

    return fullInfo;

}

std::vector<std::tuple<std::vector<int>,std::vector<int>>> linkInfoNoOneHandles(Grid grid, bool debug) {
    std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<std::vector<bool>>>> fullInfo = linkInfo(grid,true,false,debug);
    std::vector<std::tuple<std::vector<int>,std::vector<int>>> output;
    for (auto comp : fullInfo) {
        output.push_back({std::get<0>(comp),std::get<1>(comp)});
    }
    return output;
}

std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> linkInfoCanonical(Grid grid, bool debug) {
    std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<std::vector<bool>>>> fullInfo = linkInfo(grid,false,true,debug);
    std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> output;
    for (auto comp : fullInfo) {
        output.push_back({std::get<0>(comp),std::get<1>(comp),std::get<2>(comp)[0]});
    }
    return output;
}

regina::Triangulation<3> build3Manifold(Grid grid, bool simplify, bool orient) {

    auto output = linkInfoNoOneHandles(grid);

    if (output.empty()) {
        if (simplify) {
            // return the minimal 3-sphere triangulation
            return regina::Example<3>::threeSphere();
        } else {
            // for consistency return the minimal 3-sphere ~gem~ with the same standard choice
            // for gluing permutations depending on whether we are orienting
            if (orient) {
                return regina::Triangulation<3>::fromGluings(2, {
                    { 0, 0, 1, {1,0,2,3} }, { 0, 1, 1, {1,0,2,3} },
                    { 0, 2, 1, {1,0,2,3} }, { 0, 3, 1, {1,0,2,3} }});
            } else {
                return regina::Triangulation<3>::fromGluings(2, {
                    { 0, 0, 1, {0,1,2,3} }, { 0, 1, 1, {0,1,2,3} },
                    { 0, 2, 1, {0,1,2,3} }, { 0, 3, 1, {0,1,2,3} }});
            }
        }
    }

    regina::Triangulation<3> trig;
    try {
        trig = katie::katie3(output,orient);
    } catch (std::logic_error const& le) {
        std::cerr << "An error occured in Katie3\n";
        std::cerr << "grid       "; printVector(grid.x,std::cerr);
        std::cerr << "           "; printVector(grid.o,std::cerr);
        for (auto [pdcode,writhes] : output) {
            std::cerr << "pd code    "; printVector(pdcode,std::cerr);
            std::cerr << "writhes    "; printVector(writhes,std::cerr);
        }
        throw le;
    }

    if (simplify) {trig.simplify();}
    return trig;

}

regina::Triangulation<4> build4ManifoldCanonical(Grid grid, bool noOneHandles, bool simplify, bool orient) {

    auto output = linkInfoCanonical(grid);

    if (output.empty()) {
        // return the minimal 4-sphere gem, which is also a minimal triangulation
        if (orient) {
            return regina::Triangulation<4>::fromGluings(2, {
                { 0, 0, 1, {1,0,2,3,4} }, { 0, 1, 1, {1,0,2,3,4} },
                { 0, 2, 1, {1,0,2,3,4} }, { 0, 3, 1, {1,0,2,3,4} },
                { 0, 4, 1, {1,0,2,3,4} }});
        } else {
            return regina::Triangulation<4>::fromGluings(2, {
                { 0, 0, 1, {0,1,2,3,4} }, { 0, 1, 1, {0,1,2,3,4} },
                { 0, 2, 1, {0,1,2,3,4} }, { 0, 3, 1, {0,1,2,3,4} },
                { 0, 4, 1, {0,1,2,3,4} }});
        }
    }

    if (noOneHandles) {
        for (int i=0; i<output.size(); ++i) {
            size_t s = std::get<0>(output[i]).size();
            if (s==0) {
                output[i] = {std::get<0>(output[i]),std::get<1>(output[i]),std::vector<bool>(1,false)};
            } else {
                output[i] = {std::get<0>(output[i]),std::get<1>(output[i]),std::vector<bool>(s,false)};
            }
        }
    }

    regina::Triangulation<4> trig;
    try {
        trig = katie::katie4(output,orient);
    } catch (std::logic_error const& le) {
        std::cerr << "An error occured in Katie4\n";
        std::cerr << "grid       "; printVector(grid.x,std::cerr);
        std::cerr << "           "; printVector(grid.o,std::cerr);
        for (auto [pdcode,writhes,isOneHandle] : output) {
            std::cerr << "pd code    "; printVector(pdcode,std::cerr);
            std::cerr << "writhes    "; printVector(writhes,std::cerr);
            std::cerr << "onehandles "; printVector(isOneHandle,std::cerr);
        }
        throw le;
    }

    if (simplify) {trig.simplify();}
    return trig;

}

std::vector<regina::Triangulation<4>> buildAll4Manifolds(Grid grid, bool simplify, bool orient) {

    std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<std::vector<bool>>>> info = linkInfo(grid);
    std::vector<regina::Triangulation<4>> trigs;

    if (info.empty()) {
        // return the minimal 4-sphere gem, which is also a minimal triangulation
        if (orient) {
            trigs.push_back(regina::Triangulation<4>::fromGluings(2, {
                { 0, 0, 1, {1,0,2,3,4} }, { 0, 1, 1, {1,0,2,3,4} },
                { 0, 2, 1, {1,0,2,3,4} }, { 0, 3, 1, {1,0,2,3,4} },
                { 0, 4, 1, {1,0,2,3,4} }}));
            return trigs;
        } else {
            trigs.push_back(regina::Triangulation<4>::fromGluings(2, {
                { 0, 0, 1, {0,1,2,3,4} }, { 0, 1, 1, {0,1,2,3,4} },
                { 0, 2, 1, {0,1,2,3,4} }, { 0, 3, 1, {0,1,2,3,4} },
                { 0, 4, 1, {0,1,2,3,4} }}));
            return trigs;
        }
    }

    for (int oneHandleSlice=0; oneHandleSlice < std::get<2>(info[0]).size(); ++oneHandleSlice) {

        // take a slice of the possible 1-handle vectors, with the same slice index for each
        // component of the connected sum

        std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> sliceInfo;
        for (auto [pdcode,writhes,oneHandleAssignments] : info) {
            sliceInfo.push_back({pdcode,writhes,oneHandleAssignments[oneHandleSlice]});
        }

        regina::Triangulation<4> trig;
        try {
            trig = katie::katie4(sliceInfo,orient);
        } catch (std::logic_error const& le) {
            std::cerr << "An error occured in Katie4\n";
            std::cerr << "grid       "; printVector(grid.x,std::cerr);
            std::cerr << "           "; printVector(grid.o,std::cerr);
            for (auto [pdcode,writhes,isOneHandle] : sliceInfo) {
                std::cerr << "pd code    "; printVector(pdcode,std::cerr);
                std::cerr << "writhes    "; printVector(writhes,std::cerr);
                std::cerr << "onehandles "; printVector(isOneHandle,std::cerr);
            }
            throw le;
        }

        if (simplify) {trig.simplify();}
        trigs.push_back(trig);
    }

    return trigs;

}

/**
 * Fisher-Yates Shuffle
 **/
std::vector<int> randomPermutation(int n) {
    static std::random_device rd;
    static std::mt19937 mt(rd());
    static std::map<int,std::uniform_int_distribution<>> randomToK;

    for (int i=1; i<n; ++i) {
        std::uniform_int_distribution<> dist(0,i);
        randomToK[i] = dist;
    }

    std::vector<int> perm(n);

    for (int i=0; i<n; ++i) {
        perm[i] = i;
    }

    for (int i=n-1; i>0; --i) {
        int j = randomToK[i](mt);
        std::swap(perm[i],perm[j]);
    }

    return perm;
}

std::vector<int> randomDerangement(int n) {
    while (true) {
        std::vector<int> der = randomPermutation(n);
        bool valid = true;
        for (int i=0; i<n; ++i) {
            if (der[i] == i) {
                valid = false;
                break;
            }
        }
        if (valid) {
            return der;
        }
    }
}

Grid randomGrid(int n) {
    std::vector<int> x = randomPermutation(n);
    std::vector<int> xinvo = randomDerangement(n);
    return Grid(x,multiplyPerms(x,xinvo));
}

Grid randomKnot(int n) {
    std::vector<int> rfollow(n,0);
    std::vector<int> tail = randomPermutation(n-1);
    for (int i=1; i<n; ++i) {
        rfollow[i] = tail[i-1]+1;
    }
    std::vector<int> cfollow = randomPermutation(n);
    return followToCoord(rfollow,cfollow);
}

std::vector<std::pair<regina::Triangulation<3>,Grid>> sample3Manifolds(int n, int repeats, bool knotOnly, bool simplify, bool verbose) {
    std::vector<std::pair<regina::Triangulation<3>,Grid>> samples(repeats);
    for (int i=0; i<repeats; ++i) {
        Grid grid;
        if (knotOnly) {
            grid = randomKnot(n);
        } else {
            grid = randomGrid(n);
        }
        if (verbose) {std::cout << "grid "; printVector(grid.x,std::cout);
                      std::cout << "     "; printVector(grid.o,std::cout);}
        regina::Triangulation<3> trig = build3Manifold(grid,simplify);
        if (verbose) {std::cout << " -> " << trig.isoSig() << "\n";}
        samples[i]={trig,grid};
    }
    return samples;
}

std::vector<std::pair<regina::Triangulation<3>,Grid>> sampleNearVolume(double targetvol, int repeats, bool knotOnly, bool simplify, bool verbose) {
    // estimated median volume for grid sizes 10-70
    // note these are only rough estimates: taken from a single sample of 10000 at sizes 10-30,35,40,50,60,70
    // plus piecewise linear interpolation, and rounded to 4 significant figures
    // nonetheless they should be more than close enough for a proof of concept
    // TODO: more accurate estimates of these values
    std::vector<double> medians;
    if (knotOnly) {
        medians = {1.610,1.941,2.521,2.712,3.123,3.394,3.809,4.086,4.809,5.538,
                   5.888,6.426,7.361,8.458,9.510,10.94,12.22,13.65,15.73,17.65,
                   19.85,22.57,25.30,28.02,30.75,33.47,36.79,40.11,43.43,46.75,
                   50.07,54.12,58.17,62.23,66.28,70.34,74.39,78.44,82.49,86.55,
                   90.60,95.34,100.1,104.8,109.6,114.3,119.0,123.8,128.5,133.3,
                   138.0,143.3,148.6,153.9,159.2,164.5,169.7,175.0,180.3,185.6,
                   190.9};
    } else {
        medians = {2.030,2.022,2.115,2.793,3.066,3.123,3.495,4.060,4.078,5.025,
                   5.057,5.506,6.511,7.202,7.920,9.117,9.765,11.21,12.13,13.93,
                   15.45,16.52,19.18,21.93,24.76,27.67,30.67,33.75,36.91,40.16,
                   43.49,46.90,50.40,53.98,57.64,61.39,65.21,69.13,73.12,77.20,
                   81.36,85.61,89.93,94.35,98.84,103.4,108.1,112.8,117.6,122.6,
                   127.6,132.6,137.8,143.0,148.4,153.8,159.3,164.8,170.5,176.2,
                   182.0};
    }

    int n = 70;
    if (targetvol < medians[0]) {
        n = 10;
    } else {
        for (int i=1;i<medians.size();++i) {
            if (targetvol < medians[i]) {
                if (medians[i]-targetvol > targetvol-medians[i-1]) {
                    n = i + 9;
                } else {
                    n = i + 10;
                }
                break;
            }
        }
    }

    if (verbose) {std::cout << "sampling from " << (knotOnly ? "grid knots" : "grids") << " of size " << n << ", with median volume ~" << medians[n-10] << "\n";}

    std::vector<std::pair<regina::Triangulation<3>,Grid>> samples;
    while (samples.size() < repeats) {
        if (verbose) {std::cout << "|" << std::flush;}
        Grid grid;
        if (knotOnly) {
            grid = randomKnot(n);
        } else {
            grid = randomGrid(n);
        }
        regina::Triangulation<3> trig = build3Manifold(grid);
        regina::SnapPeaTriangulation strig(trig);
        double vol = strig.volume();
        if (vol > 0.9*targetvol && vol < 1.1*targetvol) {
            samples.push_back({trig,grid});
            if (verbose) {std::cout << " found! volume = " << vol << "\n";}
        }
    }
    return samples;
}

std::vector<std::pair<regina::Triangulation<4>,Grid>> sample4Manifolds(int n, int repeats, bool noOnehandles, bool closed, bool verbose) {
    std::vector<std::pair<regina::Triangulation<4>,Grid>> samples(repeats);
    for (int i=0; i<repeats; ++i) {
        while (true) {
            Grid grid = randomGrid(n);
            regina::Triangulation<4> trig = build4ManifoldCanonical(grid,noOnehandles);
            if (closed) {
                if (!trig.isClosed()) {
                    if (verbose) {std::cout << "found non-closed triangulation, retrying...\n";}
                    continue;
                }
            }
            if (verbose) {std::cout << "grid "; printVector(grid.x,std::cout);
                          std::cout << "     "; printVector(grid.o,std::cout);}
            if (verbose) {std::cout << " -> " << trig.isoSig() << "\n";}
            samples[i]={trig,grid};
            break;
        }
    }
    return samples;
}

std::string groupToManifold(std::string groupName) {
    if (groupName=="Z") {
        return "S2 x S1";
    } else if (groupName == "Z_2") {
        return "RP3";
    } else if (groupName.starts_with("Z_")) {
        std::string pString = groupName.substr(2);
        int p = std::stoi(pString);
        if (p<5) {
            return "L("+pString+",1)";
        } else {
            return "L("+pString+",?)";
        }
    } else if (groupName.starts_with("Free(")) {
        int k = std::stoi(groupName.substr(5,groupName.size()-6));
        std::string name = "S2 x S1";
        for (int i=0; i<k-1; ++i) {
            name = name + " # S2 x S1";
        }
        return name;
    } else if (groupName.starts_with("FreeProduct(")) {
        std::string component;
        std::string name;
        for (int i=12; i<groupName.size()-1; ++i) {
            char c = groupName[i];
            if (c==',') {
                name = name + groupToManifold(component) + " # ";
                component.clear();
            } else if (c!=' ') {
                component.push_back(c);
            }
        }
        name.erase(name.size()-3,3); // remove last " # "
        return name;
    } else {
        return groupName;
    }
}

std::string identify3Manifold(std::string sig) {
    
    std::string name;
    regina::Triangulation<3> trig(sig);

    trig.simplify();

    // sphere?

    if (trig.isSphere()) {
        return "S3";
    }

    // standard triangulation?

    auto recog = regina::StandardTriangulation::recognise(trig);
    if (recog) {
        name = recog->manifold()->name();
        return name;
    }

    // in census?

    std::list<regina::CensusHit> lookup = regina::Census::lookup(trig);
    if (!lookup.empty()) {
        name = lookup.front().name();
        int index = name.find(" : #");
        name.erase(index,name.size()); // trim numbering from census entry name
        if (!name.empty()) {
            return name;
        }
    }

    // connected sum?

    std::vector<regina::Triangulation<3>> summands = trig.summands();
    if (summands.size()>1) {
        for (regina::Triangulation<3> summand : summands) {
            summand.simplify();
            name = name + identify3Manifold(summand.isoSig()) + " # ";
        }
        name.erase(name.size()-3,3); // remove last " # "
        return name;
    }

    // recognise fundamental group

    std::string groupName = trig.group().recogniseGroup();

    if (!groupName.empty()) {
        return groupToManifold(groupName);
    }

    return "?";

}

bool isSplitDiagram(Grid grid) {
    // TODO: there are clearly more efficient ways to do this
    return (linkInfo(grid,true).size() > 1);
}

bool isSplitLink(Grid grid) {

    // note: PD codes can't detect zero-crossing unknot components, but we'll
    // catch them with linkInfo anyway

    auto info = linkInfo(grid); 

    if (info.empty()) {
        // this is a trivial unknot diagram
        return false;
    }

    if (info.size() != 1) {
        // this is a split diagram
        return true;
    }

    if (std::get<1>(info[0]).size() == 1) {
        // this is a knot
        return false;
    }
    
    std::vector<int> pd = std::get<0>(info[0]);
    std::vector<std::array<int,4>> pdFormat(pd.size()/4);
    for (int i=0; i<pd.size()/4; ++i) {
        for (int j=0; j<4; ++j) {
            pdFormat[i][j] = pd[4*i+j]+1;
        }
    }

    regina::Link link = regina::Link::fromPD(pdFormat.begin(),pdFormat.end());
    regina::Triangulation<3> complement = link.complement();

    // TODO: easy preconditions which confirm not split (e.g. linking number)

    // link is split iff complement has a vertex normal surface which is a sphere
    // and separates two boundary components
    // I believe cutAlong() should always produce real boundaries for newly created
    // boundary comps, so components should contain one of the original boundaries
    // iff they have an ideal vertex
    // TODO: verify this is actually true
    regina::NormalSurfaces surfaces(complement,regina::NormalCoords::Standard);
    for (const regina::NormalSurface& surf : surfaces) {
        if (surf.eulerChar()==2) {
            auto comps = surf.cutAlong().components();
            bool splits = true;
            int count = 0;
            for (auto comp : comps) {
                ++count;
                if (count > 2) {
                    splits = false;
                    break;
                }
                if (!comp->isIdeal()) {
                    splits = false;
                    break;
                }
            }
            if (count == 1) {
                splits = false;
            }
            if (splits) {
                return true;
            }
        }
    }

    return false;
}

std::vector<std::pair<std::vector<std::string>,std::vector<Grid>>> readSampleData(std::string file) {
    std::vector<std::pair<std::vector<std::string>,std::vector<Grid>>> data;
    std::ifstream storageReader;
    storageReader.open(file);
    for (std::string line; std::getline(storageReader,line); ) { // for each line in file
        std::vector<std::string> sigs;
        auto pt = line.begin();

        while (*pt != ':') {
            std::string sig;
            while (*pt != ' ') {
                sig = sig + *pt;
                pt++;
            }
            // finished reading a sig
            sigs.push_back(sig);
            pt++;
        }

        // finished reading all sigs, grids begin after next occurence of ":"
        pt++;
        while (*pt != ':') {
            pt++;
        }
        pt+=2;

        // at first character '[' of first grid
        std::vector<Grid> grids;
        while (pt != line.end()) {
            pt++;
            std::vector<int> x;
            while (*pt != '[') {
                std::string intStr;
                while (*pt != ',' && *pt != ']') {
                    intStr = intStr + *pt;
                    pt++;
                }
                x.emplace_back(std::stoi(intStr));
                pt++;
            }
            pt++; 
            std::vector<int> o;
            while (*pt != ' ') {
                std::string intStr;
                while (*pt != ',' && *pt != ']') {
                    intStr = intStr + *pt;
                    pt++;
                }
                o.emplace_back(std::stoi(intStr));
                pt++;
            }
            Grid grid(x,o);
            grids.push_back(grid);
            pt++;
        }

        data.push_back({sigs,grids});
    }

    return data;
}

// TODO: the std::array<int,3> should probably be a formTup, i.e. match the format the entries were written from
std::vector<std::pair<formTup,std::vector<Grid>>> readFormData(std::string file) {
    std::vector<std::pair<formTup,std::vector<Grid>>> data;
    std::ifstream storageReader;
    storageReader.open(file);
    for (std::string line; std::getline(storageReader,line); ) { // for each line in file
        std::array<std::string,3> entries;
        int i = 0;
        auto pt = line.begin();
        while (*pt != ':') {
            std::string entry;
            while (*pt != ' ') {
                entry = entry + *pt;
                ++pt;
            }
            // finished reading a form entry
            entries[i]=entry;
            ++pt;
            ++i;
        }
        formTup form = {std::stoi(entries[0]),std::stoul(entries[1]),std::stol(entries[2])};

        // finished reading form, grids begin after next occurence of ":"
        ++pt;
        while (*pt != ':') {
            ++pt;
        }
        pt+=2;

        // at first character '[' of first grid
        std::vector<Grid> grids;
        while (pt != line.end()) {
            ++pt;
            std::vector<int> x;
            while (*pt != '[') {
                std::string intStr;
                while (*pt != ',' && *pt != ']') {
                    intStr = intStr + *pt;
                    ++pt;
                }
                x.emplace_back(std::stoi(intStr));
                ++pt;
            }
            ++pt; 
            std::vector<int> o;
            while (*pt != ' ') {
                std::string intStr;
                while (*pt != ',' && *pt != ']') {
                    intStr = intStr + *pt;
                    ++pt;
                }
                o.emplace_back(std::stoi(intStr));
                ++pt;
            }
            Grid grid(x,o);
            grids.push_back(grid);
            ++pt;
        }

        data.push_back({form,grids});
    }

    return data;
}

std::string gridPicture(Grid grid) {
    std::vector<std::string> matrix;
    for (int r=0; r<grid.x.size(); ++r) {
        std::string row(grid.x.size()*2,' ');
        int x = grid.x[r];
        int o = grid.o[r];
        row[2*x] = 'x';
        row[2*o] = 'o';
        for (int c=2*std::min(x,o)+1; c<2*std::max(x,o); ++c) {
            row[c] = '-';
        }
        matrix.push_back(row);
    }
    std::vector<int> xcol = inverse(grid.x);
    std::vector<int> ocol = inverse(grid.o);
    for (int c=0; c<grid.x.size(); ++c) {
        int x = xcol[c];
        int o = ocol[c];
        for (int r=std::min(x,o)+1; r<std::max(x,o); ++r) {
            matrix[r][2*c] = '|';
        }
    }
    std::string picture = "* ";
    for (int c=0; c<grid.x.size(); ++c) {
        picture += "* ";
    }
    picture += "*\n";
    for (std::string row : matrix | std::views::reverse) {
        picture = picture + "* " + row + "*\n";
    }
    picture += "* ";
    for (int c=0; c<grid.x.size(); ++c) {
        picture += "* ";
    }
    picture += "*\n";
    return picture;
}

void drawGrid(Grid grid, std::ostream& stream) {
    stream << gridPicture(grid) << "\n";
}

Grid logo() {
    return Grid({{1,0,3,2,10,9,11,18,12,5,4,6,8,17,19,21,20,13,7,14,16,15,22},{8,2,1,6,9,11,3,10,0,4,7,5,13,12,20,18,21,22,15,19,14,16,17}});
}

std::string logoPictureCompact() {
    return 
        "* * * * * * * * * * * * * * * * * * * * *\n"
        "*                           o---------x *\n"
        "*     G e c K o         x-o |  #  w w | *\n"
        "*                     o-|-x |       w | *\n"
        "*       Grid          x-|---|---o     | *\n"
        "*       Kirby     x-----o   |   |   # | *\n"
        "*                 |     x---|---|-----o *\n"
        "*             x---o     |   |   | x-o   *\n"
        "*             |     o---|---x o-|-|-x   *\n"
        "*         o-x |     |   |     | x-o     *\n"
        "*       x-|-|-o x---|---o     |         *\n"
        "*       o-x |   |   |     o---x         *\n"
        "* o---------|---|---x     |             *\n"
        "* |         |   |   o-----x             *\n"
        "* |     o---|---|---|-x                 *\n"
        "* |     |   |   | x-|-o        By       *\n"
        "* |   x-|---o   | o-x         Lucy      *\n"
        "* | o-|-x       |             Tobin     *\n"
        "* x-|-o         |                       *\n"
        "*   x-----------o                       *\n"
        "* * * * * * * * * * * * * * * * * * * * *\n";
}

void drawLogoCompact(std::ostream& stream) {
    stream << logoPictureCompact() << "\n";
}

} // namespace