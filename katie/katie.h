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
 * This is an update to / expansion on Rhuaidi Burke's Katie software, which
 * implements the algorithm of Casali and Cristofori for constructing gems
 * from Kirby diagrams.
 * 
 * The original is available at 
 * github.com/raburke/Comp4Top
 * 
 * The key updates are:
 * - now a set of callable functions rather than a command line utility,
 *   importable from this header file // #include <katie.h> //
 * - can perform gem connected sums to deal with disconnected diagrams,
 *   the input format is a vector containing the data for each connected
 *   component
 * - expanded pre-processing which now comes with theoretical guarantees
 *   of correctness (as described in my PhD thesis), suitable for running
 *   on randomly generated input
 * - the 3- and 4-dimensional constructions are now separate functions
 *   katie3() and katie4(), to increase efficiency for dimension 3
 * - fixed some bugs where incorrect graph quadricolours were used,
 *   or the highlighting procedure stopped short slightly too early
 **/

#ifndef __KATIE_H
#define __KATIE_H

#include <iostream>
#include <map>
#include <vector>
#include <iterator>
#include <array>
#include <tuple>
#include <set>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <stack>
#include <string>
#include <sstream>
#include <numeric>
#include <unistd.h>
#include <string>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <link/link.h>

namespace katie {

/**
 * Build a triangulation of a 3-manifold from a Kirby diagram (integer surgery presentation).
 * 
 * The result is a gem. If orient is false, all gluings of the triangulation will use identity permutations. 
 * If orient is true, the result will be an oriented triangulation consistent with the orientation specified
 *   by the Kirby diagram, and all gluings will use the permutation 1023
 * 
 * param info: PD code and framing vector for each connected component of the diagram, in format as per gecko::linkInfoNoOneHandles()
 * param orient: return an oriented triangulation, consistent with the orientation specified by the Kirby diagram
 * 
 * throws: std::logic_error, including std::invalid_argument
 **/
regina::Triangulation<3> katie3(std::vector<std::tuple<std::vector<int>,std::vector<int>>> info, bool orient=true);

/**
 * Build a triangulation of a 4-manifold from a Kirby diagram.
 * 
 * The result is a gem. If orient is false, all gluings of the triangulation will use identity permutations. 
 * If orient is true, the result will be an oriented triangulation consistent with the orientation specified
 *   by the Kirby diagram, and all gluings will use the permutation 10234
 * 
 * param info: PD code, framing vector, and 1-handle assignments for each connected component of the diagram, in format as per gecko::linkInfo()
 * param orient: return an oriented triangulation, consistent with the orientation specified by the Kirby diagram
 * 
 * throws: std::logic_error, including std::invalid_argument
 **/
regina::Triangulation<4> katie4(std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> info, bool orient=true);

/**
 * Alternate form of katie4() which instead treats its input as a vector of Kirby diagrams for orientable (not oriented) manifolds,
 * and performs all 2^(#diagrams-1) potentially-distinct choices of how to do the connected sums
 * 
 * For connected components of a Kirby diagram, this does not make much topological sense, but it is preserved here in case it
 * becomes useful at some point
 * 
 * param info: PD code, framing vector, and 1-handle assignments for each diagram (connected component), in format as per gecko::linkInfo()
 * param orient: return an oriented triangulation, consistent with the orientation specified by the Kirby diagram
 * 
 * throws: std::logic_error, including std::invalid_argument
 **/
std::vector<regina::Triangulation<4>> katie4AllSums(std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> info, bool orient=true);

/**
 * Print to std::cout a description of the arc-coloured graph for a 3-manifold built from a Kirby diagram.
 * 
 * This is included to allow katie-cmd to print the graph without exposing katie's internal graph<d> class.
 * 
 * param info: PD code and framing vector for each connected component of the diagram, in format as per gecko::linkInfoNoOneHandles()
 * param orient: return an oriented triangulation, consistent with the orientation specified by the Kirby diagram
 * 
 * throws: std::logic_error, including std::invalid_argument
 **/
void katie3PrintGraph(std::vector<std::tuple<std::vector<int>,std::vector<int>>> info);

/**
 * Print to std::cout a description of the arc-coloured graph for a 4-manifold built from a Kirby diagram.
 * 
 * This is included to allow katie-cmd to print the graph without exposing katie's internal graph<d> class.
 * 
 * param info: PD code, framing vector, and 1-handle assignments for each connected component of the diagram, in format as per gecko::linkInfo()
 * param orient: return an oriented triangulation, consistent with the orientation specified by the Kirby diagram
 * 
 * throws: std::logic_error, including std::invalid_argument
 **/
void katie4PrintGraph(std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> info);

} // namespace

#endif