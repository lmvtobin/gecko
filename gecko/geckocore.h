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
 * Core GecKo functionality, for working with grids, building the corresponding
 * 3- and 4-manifolds, and making samples and exhaustive enumerations 
 **/

#ifndef __GECKO_CORE_H
#define __GECKO_CORE_H

#include <string>
#include <array>
#include <vector>
#include <set>
#include <map>
#include <stack>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>

#include <triangulation/dim4.h>
#include <census/census.h>

typedef std::tuple<bool,size_t,long> formTup;

namespace gecko {

struct Grid {
    std::vector<int> x;
    std::vector<int> o;
};

/**
 * Print an integer vector to console or a given output stream, as a comma-separated list in square brackets
 * 
 * param vec: input vector
 * param stream: output stream
 * param newline: print a newline character at end of output
 **/
void printVector(std::vector<int> vec, std::ostream& stream=std::cout, bool newline=true);

/**
 * Print a boolean vector to console or a given output stream, as a comma-separated list in square brackets
 * 
 * param vec: input vector
 * param stream: output stream
 * param newline: print a newline character at end of output
 **/
void printVector(std::vector<bool> vec, std::ostream& stream=std::cout, bool newline=true);

/**
 * Compute the inverse of a permutation
 * 
 * param perm: permutation to be inverted
**/
std::vector<int> inverse(std::vector<int> perm);

/**
 * Multiply two permutations
 * Uses function composition order: p[x]=p1[p2[x]]
 * 
 * params p1,p2: permutations to be multiplied
**/
std::vector<int> multiplyPerms(std::vector<int> p1, std::vector<int> p2);

/**
 * Returns whether the x and o permutations of the grid have no common entries
 * 
 * param grid: input grid
**/
bool noCollisions(Grid grid);

/**
 * Return the lexicographically smallest derangement on {0,1,...,n-1}
 * 
 * param n: grid size
**/
std::vector<int> firstDerangement(int n);

/**
 * Return the lexicographically smallest grid of size n with x < o
 * 
 * param n: grid size
**/
Grid firstGrid(int n);

/** 
 * Changes a grid into the next grid by lexicographic order of permutations, with x < o
 * Returns true if such a grid exists
 * If no such grid exists, returns false and changes the grid to the first in this order
 * 
 * note: return behaviour is modeled on the behaviour of std::next_permutation
 * 
 * param grid: input grid
**/
bool nextGrid(Grid* grid);

/**
 * Convert a grid in standard representation to "follow-around" representation
 * follow-around is as follows:
 * - order X's as they are encountered traversing around the link, starting with
 *   the X in the lowest row not yet encountered and moving horizontally first,
 *   until the component is fully traversed, and repeat
 * - first two outputs rfollow,cfollow are permutations, such that the ith X
 *   in this order has coordinates (rfollow[i],cfollow[i])
 * - third output records the size of each component (number of X's) in order
 *   encountered, forming a partition of n (grid size)
 * 
 * param grid: input grid
 **/
std::tuple<std::vector<int>,std::vector<int>,std::vector<int>> coordToFollow(Grid grid);

/**
 * Convert a grid in follow-around representation to standard representation
 * 
 * params: as in output of coordToFollow()
 **/
Grid followToCoord(std::vector<int> rfollow, std::vector<int> cfollow, std::vector<int> partition);

/**
 * Specialisation of followToCoord() when the grid is a knot, with only one component
 * (i.e. when partition={n} is a partition into one element of size n)
 * 
 * params: as in first two outputs of coordToFollow()
 **/
Grid followToCoord(std::vector<int> rfollow, std::vector<int> cfollow);

/**
 * Find the connected components of a graph
 * Returns a vector, each entry is the subset of vertices corresponding to a component
 * 
 * param graph: input graph with vertices labelled by integers, as a map from vertex to set of neighbours
 **/
std::vector<std::set<int>> connectedComponents(std::map<int,std::set<int>> graph);

/**
 * For each connected component of the grid (as a graph), compute the three pieces of information needed to feed into katie
 * 1. pd code (concatenated as a single vector)
 * 2. writhe of each component (framing vector for blackboard framing)
 * 3. a vector containing all possible valid assignments of 1-handles
 *    (each combination is a boolean vector of true for 1-handles, false for 2-handles)
 * 
 * Components are ordered left-to-right by the leftmost column they occupy
 * 
 * If canonicalOneHandles or noOneHandles is true, the third vector will instead contain only a single entry, which is:
 *  noOneHandles: all components are 2-handles
 *  canonicalOneHandles && !noOneHandles: the lexicographically largest vector
 *      (i.e. iterate in left-to-right order and always choose to include as a 1-handle if valid)
 * 
 * note: the whole output is an empty vector if and only if the grid size is strictly less than 2 (to be interpreted as an empty link)
 * note: the pd code returned is an empty vector if and only if that component is a single 0-crossing unknot
 * 
 * note: if there are multiple connected components, 1-handle vectors should be taken in slices to get all possibilities
 *       e.g. take connected sum of output[0][0],output[0][1],output[0][2][i] and output[1][0],output[1][1],output[1][2][i] for each i
 *       (taking all pairs of third indices will cause unnecessary duplicates)
 * 
 * param grid: imput grid
 * param noOneHandles: skip 1-handle computations and return only the assignment consisting of only 2-handles
 * param canonicalOneHandles: return only one canonical assignment of 1-handles
 * param debug: print debugging info to cerr
**/
std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<std::vector<bool>>>> linkInfo(Grid grid, bool noOneHandles=false, bool canonicalOneHandles=false, bool debug=false);

/**
 * Condensed output of linkInfo() for ease of use in the case when 1-handles are ignored
 * Returns only pd code and writhe for each component
 * 
 * param grid: input grid
 * param debug: print debugging info to cerr
**/
std::vector<std::tuple<std::vector<int>,std::vector<int>>> linkInfoNoOneHandles(Grid grid, bool debug=false);

/**
 * Condensed output of linkInfo() for ease of use in the case of canonical 1-handle assignment
 * Returns only pd code, writhe and 1-handle vector for each component
 * (instead of wrapping 1-handles inside a vector of size 1 as in the full linkInfo())
 * 
 * param grid: input grid
 * param debug: print debugging info to cerr
**/
std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> linkInfoCanonical(Grid grid, bool debug=false);

/**
 * Build a triangulation of the 3-manifold with this grid as its Kirby / integer surgery diagram
 * If the grid is a split diagram, perform connected sums in the canonical way specified in Katie
 * 
 * If simplify is false, the result is a gem, and orient determines the gluing permutations used as specified in katie::katie3(), i.e.
 * - if orient is false, all gluings will use identity permutations
 * - if orient is true, the result will be the correct oriented triangulation as specified by katie::katie3() (with all permutations 1023)
 * If simplify is true, orientation will be preserved, but not the nice gem structure of consistent permutations
 * 
 * param grid: input grid
 * param simpify: attempt to simplify the triangulation with Regina after constructing
 * param orient: return an oriented triangulation, consistent with the orientation specified by the Kirby diagram
 **/
regina::Triangulation<3> build3Manifold(Grid grid, bool simplify=true, bool orient=true);

/**
 * Build a triangulation of the 4-manifold with this grid as its Kirby diagram, with a single canonical choice of 1-handle assignment
 * If the grid is a split diagram, perform connected sums in the canonical way specified in Katie
 * 
 * By default uses the canonical 1-handle assignment from linkInfo(), i.e. the lexicographically largest 1-handle assignment vector
 * If noOneHandles is true, interprets all link components as 2-handles
 * 
 * If simplify is false, the result is a gem, and orient determines the gluing permutations used as specified in katie::katie3(), i.e.
 * - if orient is false, all gluings will use identity permutations
 * - if orient is true, the result will be the correct oriented triangulation as specified by katie::katie3() (with all permutations 1023)
 * If simplify is true, orientation will be preserved, but not the nice gem structure of consistent permutations
 * 
 * param grid: input grid
 * param noOneHandles: interpret all components as 2-handles
 * param simpify: attempt to simplify the triangulation with Regina after constructing
 * param orient: return an oriented triangulation, consistent with the orientation specified by the Kirby diagram
 **/
regina::Triangulation<4> build4ManifoldCanonical(Grid grid, bool noOneHandles=false, bool simplify=true, bool orient=true);

/**
 * Build all triangulations of 4-manifolds with this grid as its Kirby diagram, over all valid assignments of 1-handles
 * 
 * If simplify is false, the results are gems, and orient determines the gluing permutations used as specified in katie::katie3(), i.e.
 * - if orient is false, all gluings will use identity permutations
 * - if orient is true, the result will be the correct oriented triangulation as specified by katie::katie3() (with all permutations 1023)
 * If simplify is true, orientation will be preserved, but not the nice gem structure of consistent permutations
 * 
 * param grid: input grid
 * param simpify: attempt to simplify the triangulations with Regina after constructing
 * param orient: return oriented triangulations, consistent with the orientation specified by the Kirby diagram
 **/
std::vector<regina::Triangulation<4>> buildAll4Manifolds(Grid grid, bool simplify=true, bool orient=true);

/**
* Generate a uniformly random permutation of 0,1,...,n-1
**/
std::vector<int> randomPermutation(int n);

/**
* Generate a uniformly random derangement of 0,1,...,n-1, a permutation with no fixed points
**/
std::vector<int> randomDerangement(int n);

/**
* Generate a uniformly random grid of size n by n
**/
Grid randomGrid(int n);

/**
* Generate a uniformly random grid of size n by n with only one link component
**/
Grid randomKnot(int n);

/**
* Generate a specified number of random grids and construct the corresponding 3-manifold triangulations
* Returns a vector of pairs (triangulation,grid)
* 
* param n: size of the grids
* param repeats: number of grids to sample
* param knotOnly: only sample from grids which have one link component
* param verbose: print samples to cout
**/
std::vector<std::pair<regina::Triangulation<3>,Grid>> sample3Manifolds(int n, int repeats, bool knotOnly=false, bool simplify=false, bool verbose=false);

/**
* Generate random 3-manifold triangulations with hyperbolic volume between 0.9 and 1.1 times the given target volume,
* from grids of an appropriate grid size
* Returns a vector of pairs (triangulation,grid)
* 
* note: target volume should be between 2 and 200
* 
* param targetvol: target volume
* param repeats: number of triangulations to sample
* param knotOnly: only sample from grids which have one link component
* param verbose: print samples to cout
**/
std::vector<std::pair<regina::Triangulation<3>,gecko::Grid>> sampleNearVolume(double targetvol, int repeats, bool knotOnly=false, bool simplify=false, bool verbose=false);

/**
* Generate a specified number of random grids and construct the corresponding 4-manifold triangulations
* Returns a vector of pairs (triangulation,grid)
* 
* param n: size of the grids
* param repeats: number of grids to sample
* param noOneHandles: interpret all components as 2-handles
* param closed: only sample closed manifold (that is, Kirby diagrams which give closed manifolds without needing any 3-handles)
* param verbose: print samples to cout
**/
std::vector<std::pair<regina::Triangulation<4>,Grid>> sample4Manifolds(int n, int repeats, bool noOneHandles=false, bool closed=false, bool verbose=false);

/**
 * Attempt to identify a 3-manifold from a fundamental group which has already been identified and named by Regina
 * If unsuccessful, return back the given group name
 * 
 * param groupName: name of the group, in format of the output of regina::Triangulation<3>.group().recogniseGroup()
 **/
std::string groupToManifold(std::string groupName);

/**
 * Attempt to identify a 3-manifold from a triangulation, given as an isomorphism signature
 * Uses Regina's census lookup, connected sum decomposition, and fundamental group recognition
 * 
 * If not prime, return summands separated by " # ", e.g. "RP3 # S2 x S1"
 * 
 * If for a summand the group could be identified but the manifold could not (e.g. lens spaces where census/simplification failed),
 * use the group name instead for that component (e.g. "RP3 # Z_7")
 * 
 * If neither manifold nor group could be identified, use a "?" for that component (e.g. "RP3 # ?")
 * 
 * warning: not thread safe due to Regina's census lookup
 * 
 * param sig: isomorphism signature of a 3-manifold triangulation
 **/
std::string identify3Manifold(std::string sig);

/**
 * Identify if a grid diagram is split, i.e. the underlying 4-regular graph has more than one connected component
 * 
 * param grid: input grid
 **/
bool isSplitDiagram(Grid grid);

/**
 * Identify if a grid diagram represents a split link
 * 
 * param grid: input grid
 **/
bool isSplitLink(Grid grid);

/**
 * Read in data saved in a text file, in the format the programs sample3 and exhaustive3 save in
 * 
 * Format is one line for each triangulation, of the form
 * "sig1 sig2 ... sigk : number_of_grids : grid1 grid2 grid3 ... gridm"
 * each sig (representing one connected summand) is a Regina isomorphism signature 
 * each grid is two comma-separated lists (x and o) with square brackets and no spaces between them
 * e.g. "[0,2,1][1,0,2]"
 * 
 * Output is a vector with one entry for each sig, containing a pair (list_of_sigs,list_of_grids)
 * 
 * param file: path to file
 **/
std::vector<std::pair<std::vector<std::string>,std::vector<Grid>>> readSampleData(std::string file);

/**
 * Read in data saved in a text file, in the format the program exhaustive4 saves in for closed 
 * simply connected 4-manifolds, using intersection forms
 * 
 * Format is one line for each triangulation, of the form
 * "form.odd() form.rank() form.signature() : number_of_grids : grid1 grid2 grid3 ... gridm"
 * the three values specifying the intesection form are all integers, to be interpreted as a
 * bool, size_t, and long respectively (i.e. a formTup)
 * each grid is two comma-separated lists (x and o) with square brackets and no spaces between them
 * e.g. "[0,2,1][1,0,2]"
 * 
 * Output is a vector with one entry for each intersection form, containing a pair (form,list_of_grids)
 * 
 * param file: path to file
 **/
std::vector<std::pair<formTup,std::vector<Grid>>> readFormData(std::string file);

/**
 * Make an ASCII art representation of a grid, as a string
 * 
 * param grid: input grid
 **/
std::string gridPicture(Grid grid);

/**
 * Print an ASCII art representation of a grid to console, or a given output stream
 * 
 * param grid: input grid
 * param stream: output stream
 **/
void drawGrid(Grid grid, std::ostream& stream=std::cout);

/**
 * The GecKo logo as a grid
 **/
Grid logo();

/**
 * The GecKo logo drawn in ASCII art more compactly with decorations, not quite as a grid
 **/
std::string logoPictureCompact();

/**
 * Print the ASCII GecKo logo from logoPictureCompact() to console, or a given output stream
 **/
void drawLogoCompact(std::ostream& stream=std::cout);

} // namespace

#endif