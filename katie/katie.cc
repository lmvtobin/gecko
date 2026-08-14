//
// Katie
// Kirby Diagrams to Graphs and Triangulations
//
// Created by Rhuaidi Antonio Burke on 17/05/24.
// Copyright © 2024 Regina Development Team. All rights reserved.
//

/**
 * Edited by Lucy Tobin as part of GecKo.
 * Implementation for katie.h
 **/

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
#include <bitset>
#include <ranges>

#include <chrono>
#include <stdexcept>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <link/link.h>

#include "katie.h"

bool printDebugInfo = false;

typedef std::vector<std::array<int, 4>> pdcode;

std::ostream& operator<<(std::ostream& os, const std::array<int, 4>& arr) {
    os << "(" << arr[0] << ", " << arr[1] << ", " << arr[2] << ", " << arr[3] << ")";
    return os;
}

struct node {
    int nodeID = -1;
    int strand = -1; // This comes from the PD code
    int subgraphComponent = -1; // This is a component ID w.r.t to the subgraphs.
};

const node emptyNode = node{-1,-1,-1};

struct edge {
    node n1 {};
    node n2 {};
    int colour {};
};

std::ostream& operator<<(std::ostream& os, const node& n) {
    os << "(" << n.nodeID << ", " << n.subgraphComponent << ")";
    return os;
}

bool operator ==(const node& x, const node& y) {
    return ((x.nodeID == y.nodeID) && (x.strand == y.strand) && (x.subgraphComponent == y.subgraphComponent));
}

bool operator <(const node& x, const node& y) {
    return std::tie(x.nodeID, x.strand, x.subgraphComponent) < std::tie(y.nodeID, y.strand, y.subgraphComponent);
}

long getIndex(std::vector<node> n, node K) {
    // 2022 Rewrite: (Doc) Given a list of nodes n, and a node K, return the index of K in n.
    long ans = -1;
    
    auto it = std::find_if(n.begin(),n.end(), [&tNode = K](const node& cNode) -> bool {return tNode == cNode;});
    
    if (it != n.end()) {
        long index = it - n.begin();
        ans = index;
    }
    else {
        ans = -1;
    }
    
    return ans;
}

template <typename T>
bool contains(std::vector<T> a, std::vector<T> b) {
    for (const auto& a_element : a) {
        if (std::find(b.begin(),b.end(), a_element) == b.end()) {
            return false;
        }
    }
    return true;
}

bool isCurl(const regina::StrandRef &ref) {
    
    bool ans = false;
    
    long refCrossingIndex = ref.crossing()->index();
    long nextRefCrossingIndex = ref.next().crossing()->index();
    long prevRefCrossingIndex = ref.prev().crossing()->index();

    // Doing this very "safely"...
    if (refCrossingIndex == nextRefCrossingIndex) {
        ans = true;
    }
    if (refCrossingIndex == prevRefCrossingIndex) {
        ans = true;
    }
    
    return ans;
}

template <int dim>
class graph {
    std::map<node, std::array<node,dim+1>> adjList;
    node orientation; // LT: distinguish a node which defines an orientation / bipartition class
    // TODO (LT): clean up code to prevent possibility of this being unassigned

public:

    node getOrientation() const {
        return orientation;
    }

    void setOrientation(node n) {
        orientation = n;
    }
    
    std::map<node, std::array<node,dim+1>> adjacencyList() const {
        return adjList;
    }
    
    void fromAdjacencyList(const std::map<node, std::array<node,dim+1>>& graphData) {
        adjList = graphData;
    }
        
    void addEdge(edge e, bool force=false) {
        // LT: now checks we don't already have an edge of this colour at either endpoint
        // unless we specify we should ignore this and force addition (useful e.g. in connected sum)
        if (force || adjList[e.n1][e.colour] == emptyNode && adjList[e.n2][e.colour] == emptyNode) {
            adjList[e.n1][e.colour] = e.n2;
            adjList[e.n2][e.colour] = e.n1;
        }
    }
    
    void addEdges(const std::vector<edge>& el) {
        for (const auto& e : el) {
            addEdge(e);
        }
    }
    
    std::vector<node> nodes() const {
        std::vector<node> nodeList;
        
        for (const auto& [key,val] : adjList) {
            nodeList.emplace_back(key);
        }
        
        return nodeList;
    }
    
    std::vector<edge> edges() const {
        std::vector<edge> edgeList;
        
        for (const auto& [node, nbrs] : adjList) {
            for (int i=0; i<dim+1; i++) {
                if ((node < nbrs[i]) && (node.nodeID != 0)) {
                    edgeList.push_back({node,nbrs[i],i});
                }
            }
        }
        
        return edgeList;
    }

     std::vector<edge> edgesDirected() const {
        std::vector<edge> edgeList;
        
        for (const auto& [node, nbrs] : adjList) {
            for (int i=0; i<dim+1; i++) {
                if (node.nodeID != 0) {
                    edgeList.push_back({node,nbrs[i],i});
                }
            }
        }
        
        return edgeList;
    }
    
    void disjoint_union(graph<dim> h, int &nextID) {
        std::vector<edge> hEdges = h.edges();
        
        int currentID = nextID;
        nextID++;
        
        for (auto e : hEdges) {
            node n1, n2;
            n1 = e.n1;
            n2 = e.n2;
            n1.subgraphComponent = currentID;
            n2.subgraphComponent = currentID;

            int col = e.colour;
            
            adjList[n1][col] = n2;
            adjList[n2][col] = n1;
        }
    }

    void disjoint_union_external(graph<dim> h, int lastID) {
        std::vector<edge> hEdges = h.edges();
        
        for (auto e : hEdges) {
            node n1, n2;
            n1 = e.n1;
            n2 = e.n2;
            n1.subgraphComponent = n1.subgraphComponent + lastID;
            n2.subgraphComponent = n2.subgraphComponent + lastID;

            int col = e.colour;
            
            adjList[n1][col] = n2;
            adjList[n2][col] = n1;
        }
    }
    
    size_t size() {
        return adjList.size();
    }
    
    void print() {
        for (const auto& [n,nbrs] : adjList) {
            for (int i=0; i<dim+1; i++) {
                if ((n < nbrs[i]) && (n.nodeID != 0)) {
                    std::cout << "[" << n << "," << nbrs[i] << "," << i << "],\n";
                }
            }
        }
    }
    
    void printNodes() {
        for (const auto& [key,val] : adjList) {
            std::cout << key << "\n";
        }
        std::cout << std::endl;
    }
    
    void pdSub(const pdcode& code) {
        /*
         Debating whether to make this a
         standalone function which operates
         on a graph, rather than a method
         of the graph class...
         */
        
        node newNbr;
        
        for (auto [currNode, nbrs] : adjList) {
            for (int i=0; i<4; i++) {
                if (nbrs[i].nodeID != 0) {
                    edge subbedEdge;
                    switch (nbrs[i].strand) {
                        case 0:
                            break;
                        case 1:
                            newNbr = {nbrs[i].nodeID, code[currNode.subgraphComponent][0],nbrs[i].subgraphComponent};
                            adjList.erase(nbrs[i]);
                            adjList[currNode][i] = newNbr;
                            adjList[newNbr][i] = currNode;
                            break;
                        case 2:
                            newNbr = {nbrs[i].nodeID, code[currNode.subgraphComponent][1],nbrs[i].subgraphComponent};
                            adjList.erase(nbrs[i]);
                            adjList[currNode][i] = newNbr;
                            adjList[newNbr][i] = currNode;
                            break;
                        case 3:
                            newNbr = {nbrs[i].nodeID, code[currNode.subgraphComponent][2],nbrs[i].subgraphComponent};
                            adjList.erase(nbrs[i]);
                            adjList[currNode][i] = newNbr;
                            adjList[newNbr][i] = currNode;
                            break;
                        case 4:
                            newNbr = {nbrs[i].nodeID, code[currNode.subgraphComponent][3],nbrs[i].subgraphComponent};
                            adjList.erase(nbrs[i]);
                            adjList[currNode][i] = newNbr;
                            adjList[newNbr][i] = currNode;
                            break;
                    }
                }
            }
        }
    }
    
    std::vector<std::pair<node,node>> fuseList() {
        /*
         Debating whether to make this a
         standalone function which operates
         on a graph, rather than a method
         of the graph class...
         */

        /*
         Let N_i = (c_i, n_i, s_i), V_j = (c_j, n_j, s_j)
         Criteria in the if statement below are as follows:
         1. Avoids duplicate pairs (works because elements are ordered).
         2. Only operate on "outer" nodes ("internal" nodes denoted via s = 0).
         3. c_i ≠ c_j (different "components")
         4. s_i = s_j (same strand/PD element)
         5. n_j mod 4 = (5 - (n_i mod 4)) mod 4.
         */
        
        std::vector<std::pair<node,node>> result;
        
        for (auto const& [n1,nbrs1] : adjList) {
            for (auto const& [n2,nbrs2] : adjList) {
                if (
                    (n1.subgraphComponent < n2.subgraphComponent) &&
                    (n1.strand !=0) && (n2.strand != 0) &&
                    (n1.strand == n2.strand) &&
                    ((n1.nodeID)%4 == (5-((n2.nodeID)%4))%4)
                    ) {
                        result.emplace_back(n1,n2);
                }
            }
        }

        // LT: this is specifically to deal with the case when the link has only one crossing
        // in this case, we do need to allow fusing within the same "component"

        if (result.empty()) {
            for (auto const& [n1,nbrs1] : adjList) {
                for (auto const& [n2,nbrs2] : adjList) {
                    if (
                        (n1 != n2) &&
                        (n1.strand !=0) && (n2.strand != 0) &&
                        (n1.strand == n2.strand) &&
                        ((n1.nodeID)%4 == (5-((n2.nodeID)%4))%4)
                        ) {
                            result.emplace_back(n1,n2);
                    }
                }
            }
        }
        
        return result;
    }
    
    void fuse(node n1, node n2) {
        std::array<node,dim+1> n1nbrs, n2nbrs;
        n1nbrs = adjList[n1];
        n2nbrs = adjList[n2];
        
        adjList.erase(n1);
        adjList.erase(n2);
        
        for (int i=0; i<dim+1; i++) {
            adjList[n1nbrs[i]][i] = n2nbrs[i];
            adjList[n2nbrs[i]][i] = n1nbrs[i];
        }
        
        adjList.erase({0,0,0});
    }
    
    void addQuadriEdges(const std::vector<std::array<node,4>>& quadriVect) {
        if (printDebugInfo) {
            std::clog << "Adding quadricolour edges..." << std::endl;
        }
        for (const auto& quadri : quadriVect) {
            if (printDebugInfo) {std::clog << "  adding (" << getIndex(nodes(),quadri[0]) << "," << getIndex(nodes(),quadri[1]) << ") [" << quadri[0] << "," << quadri[1] << "]\n";}
            adjList[quadri[0]][4] = quadri[1];
            adjList[quadri[1]][4] = quadri[0];
            
            if (printDebugInfo) {std::clog << "  adding (" << getIndex(nodes(),quadri[2]) << "," << getIndex(nodes(),quadri[3]) << ") [" << quadri[2] << "," << quadri[3] << "]\n";}
            adjList[quadri[2]][4] = quadri[3];
            adjList[quadri[3]][4] = quadri[2];
            
            node P4 = adjList[quadri[3]][1];
            node P5 = adjList[quadri[0]][1];

            if (printDebugInfo) {std::clog << "  adding (" << getIndex(nodes(),P4) << "," << getIndex(nodes(),P5) << ") [" << P4 << "," << P5 << "]\n";}            
            adjList[P4][4] = P5;
            adjList[P5][4] = P4;
        }
    }
    
    void addDoubleOneEdges() {
        if (printDebugInfo) {
            std::clog << "Adding doubled 1-coloured edges..." << std::endl;
        }
        for (const auto& [key,nbrs] : adjList) {
            if (
                (key < nbrs[1]) &&
                (adjList[key][4] == emptyNode) &&
                (adjList[nbrs[1]][4] == emptyNode)
                ) {
                    edge e = {key,nbrs[1],4};
                    if (printDebugInfo) {std::clog << "  adding (" << getIndex(nodes(),key) << "," << getIndex(nodes(),nbrs[1]) << ") [" << key << "," << nbrs[1] << "]\n";}
                    addEdge(e);
            }
        }
    }
    
    void addOneHandleMarkerEdges(std::vector<std::pair<node,node>> markerNodePairs) {
        if (printDebugInfo) {
            std::clog << "Adding 1-handle marked edges..." << std::endl;
        }
        for (const auto& pair : markerNodePairs) {
            edge e = {pair.first,pair.second,4};
            if (printDebugInfo) {std::clog << "  adding (" << getIndex(nodes(),pair.first) << "," << getIndex(nodes(),pair.second) << ") [" << pair.first << "," << pair.second << "]\n";}
            addEdge(e);
        }
        if (printDebugInfo) {
            std::clog << "Successfully added 1-handle marked edges!" << std::endl;
        }
    }

    /* LT: new function, this is a sanity check to verify the marked edges correspond to valid rho_3 switches
     the point of the highlighting procedure is to (implicitly) do rho_2 switches which make these rho_3
     switches valid, so after highlighting we double check this worked
    */
    bool checkOneHandleRhoThrees(std::vector<std::pair<node,node>> markerNodePairs) {
        if (printDebugInfo) {
            std::clog << "Confirming 1-handle marked edges are valid rho_3 switches..." << std::endl;
        }
        bool allValid = true;
        int i = 0;
        for (const auto& pair : markerNodePairs) {
            for (int c : std::vector({0,2,3})) {
                bool valid = true;
                node current = pair.first;
                while (true) {
                    current = adjList[current][c];
                    if (current == pair.second) {
                        break;
                    }
                    if (adjList[current][4] == emptyNode) {
                        current = adjList[current][1];
                    } else {
                        current = adjList[current][4];
                    }
                    if (current == pair.first) {
                        valid = false;
                        break;
                    }
                }
                if (!valid) {
                    allValid = false;
                    if (printDebugInfo) {
                        std::clog << "  1-handle " << i << ", colours (" << c << ",4): error, marked vertices in different residues!" << std::endl;
                    }
                }
                else if (printDebugInfo) {
                    std::clog << "  1-handle " << i << ", colours (" << c << ",4): same residue" << std::endl;
                }
            }
            ++i;
        }
        return allValid;
    }
    
    void addHighlightEdges(std::vector<std::vector<regina::StrandRef>> highlightCrossings) {
        if (printDebugInfo) {
            std::clog << "Adding highlight edges..." << std::endl;
        }
        std::vector<node> allNodes = nodes();
        
        std::vector<std::vector<node>> highlightOverNodes, highlightUnderNodes, highlightCurlNodes;
        
        for (const auto& highlightVect : highlightCrossings) {
            for (const auto& ref : highlightVect) {
                std::vector<node> currCrossingNodes;
                for (const auto& n : allNodes) {
                    if (n.subgraphComponent == ref.crossing()->index()) {
                        currCrossingNodes.emplace_back(n);
                    }
                }
                if ((ref.strand() == 0) && !(isCurl(ref))) {
                    highlightUnderNodes.emplace_back(currCrossingNodes);
                }
                else if ((ref.strand() == 1) && !(isCurl(ref))) {
                    highlightOverNodes.emplace_back(currCrossingNodes);
                }
                else if (isCurl(ref)) {
                    highlightCurlNodes.emplace_back(currCrossingNodes);
                }

            }
        }

        // Over
        if (printDebugInfo) {std::clog << "  Overcrossings...\n";}
        for (const auto& vect : highlightOverNodes) {
            for (const auto& x : vect) {
                for (const auto& y : vect) {
                    if (x < y) {
                        if (
                            ((x.nodeID == 1) && (y.nodeID == 2)) ||
                            ((x.nodeID == 5) && (y.nodeID == 6))
                            ) {
                            edge e = {x,y,4};
                            if (adjList[x][4] != emptyNode || adjList[y][4] != emptyNode) {
                                if (printDebugInfo) {std::clog << "  not adding (" << getIndex(allNodes,x) << "," << getIndex(allNodes,y) << ") [" << x << "," << y << "] already added an edge to an endpoint\n";}
                            } else {
                                if (printDebugInfo) {std::clog << "  adding (" << getIndex(allNodes,x) << "," << getIndex(allNodes,y) << ") [" << x << "," << y << "]\n";}
                                addEdge(e);
                            }
                        }
                    }
                }
            }
        }
        
        // Under
        if (printDebugInfo) {std::clog << "  Undercrossings...\n";}
        for (const auto& vect : highlightUnderNodes) {
            for (const auto& x : vect) {
                for (const auto& y : vect) {
                    if (x < y) {
                        if (
                            ((x.nodeID == 1) && (y.nodeID == 6)) ||
                            ((x.nodeID == 2) && (y.nodeID == 5)) ||
                            ((x.nodeID == 3) && (y.nodeID == 4)) ||
                            ((x.nodeID == 7) && (y.nodeID == 8))
                            ) {
                            edge e = {x,y,4};
                            if (adjList[x][4] != emptyNode || adjList[y][4] != emptyNode) {
                                if (printDebugInfo) {std::clog << "  not adding (" << getIndex(allNodes,x) << "," << getIndex(allNodes,y) << ") [" << x << "," << y << "] already added an edge to an endpoint\n";}
                            } else {
                                if (printDebugInfo) {std::clog << "  adding (" << getIndex(allNodes,x) << "," << getIndex(allNodes,y) << ") [" << x << "," << y << "]\n";}
                                addEdge(e);
                            }
                        }
                    }
                }
            }
        }

        // Curl
        if (printDebugInfo) {std::clog << "  Curls...\n";}
        for (const auto& vect : highlightCurlNodes) {
            for (const auto& x : vect) {
                for (const auto& y : vect) {
                    if (x < y) {
                        if ((adjList[x][4] == emptyNode) && (adjList[y][4] == emptyNode)) {
                            if (
                                ((x.nodeID == 1) && (y.nodeID == 4)) ||
                                ((x.nodeID == 2) && (y.nodeID == 3))
                                ) {
                                edge e = {x,y,4};
                                if (adjList[x][4] != emptyNode || adjList[y][4] != emptyNode) {
                                    if (printDebugInfo) {std::clog << "  not adding (" << getIndex(allNodes,x) << "," << getIndex(allNodes,y) << ") [" << x << "," << y << "] already added an edge to an endpoint\n";}
                                } else {
                                    if (printDebugInfo) {std::clog << "  adding (" << getIndex(allNodes,x) << "," << getIndex(allNodes,y) << ") [" << x << "," << y << "]\n";}
                                    addEdge(e);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // if (printDebugInfo) {
        //     std::clog << "Successfully added highlight edges!" << std::endl;
        // }

    }
    
    void addRemainderEdges() {
        if (printDebugInfo) {
            std::clog << "Adding remainder edges..." << std::endl;
        }
        std::vector<node> allNodes = nodes();
        for (const auto& x : allNodes) {
            if (!(x == emptyNode)) {
                if (adjList[x][4] == emptyNode) {
                    node y = x;
                    int i;
                    int j = 0;
                    do {
                        if (j > nodes().size()) {
                            throw std::logic_error("katie appears to have hit an infinite loop in addRemainderEdges()");
                        }
                        i = 4*(j%2)+(j+1)%2;
                        bool isEdge = false;
                        for (edge e : edges()) {
                            if (e.colour == i && (e.n1 == y && e.n2 == adjList[y][i] || e.n2 == y && e.n1 == adjList[y][i])) {
                                isEdge = true;
                            }
                        }
                        y = adjList[y][i];
                        j+=1;
                    } while (!(adjList[y][4] == emptyNode));
                    edge e = {x,y,4};
                    if (printDebugInfo) {std::clog << "  adding (" << getIndex(allNodes,x) << "," << getIndex(allNodes,y) << ") [" << x << "," << y << "]\n";}
                    addEdge(e);
                }
            }
        }
        // if (printDebugInfo) {
        //     std::clog << "Successfully added remainder edges!" << std::endl;
        // }
    }
    
    void DEBUGremainingNoCol4Nodes() {
        std::vector<node> allNodes = nodes();
        int counter = 0;
        for (const auto& n : allNodes) {
            if (adjList[n][4] == emptyNode) {
                counter++;
//                std::clog << n << " has no colour 4 neighbour.\n";
            }
        }
        std::clog << "Remaining nodes without a colour 4 edge: " << counter << std::endl;
    }
    
    void cleanup() {
        adjList.erase(emptyNode);
    }

    /* LT: new function
     we do the connected sum with:
     - the node 0-adjacent to the orientation node in this graph
     - either the orientation node or its 0-neighbour in the other graph, depending on reversing / non reversing
     this ensures the orientation node is never removed from the graph: it still defines the orientation
    */
    void connected_sum(graph<dim> h, bool reverse, int &nextID) {
        node n1 = adjList[orientation][0];
        node n2 = h.getOrientation();
        if (reverse) {
            n2 = h.adjacencyList()[n2][0];
        }
        int oldID = nextID;
        disjoint_union_external(h,nextID);
        for (node n : nodes()) {
            if (n.nodeID == n2.nodeID && n.strand == n2.strand && n.subgraphComponent == n2.subgraphComponent+oldID) {
                n2 = n;
                break;
            }
        }
        for (int i=0; i<dim+1; ++i) {
            addEdge({adjList[n1][i],adjList[n2][i],i},true);
        }
        adjList.erase(n1);
        adjList.erase(n2);
        cleanup();
    }
};

std::vector<std::array<node, 4>> findGraphQuadricolours(const graph<4>& G) {
    std::vector<std::array<node, 4>> ans;
    
    std::map<node, std::array<node, 5>> adjList = G.adjacencyList();
    
    std::array<node, 2> tmp;
    
    for (const auto& [n, nbrs] : adjList) {
        if (!(n == emptyNode)) {
            tmp[0] = adjList[n][0];
            tmp[1] = adjList[n][3];
            if (adjList[tmp[0]][1] == adjList[tmp[1]][2]) {
                std::array<node, 4> currentQuadri;
                currentQuadri[0] = n;
                currentQuadri[1] = tmp[0];
                currentQuadri[2] = adjList[tmp[0]][1];
                currentQuadri[3] = tmp[1];
                ans.emplace_back(currentQuadri);
            }
        }
    }
    
    return ans;
}

void walkAroundLink(regina::Link lnk) {
    std::cout << "Debug link walkaround:" << std::endl;
    for (const auto& comp : lnk.components()) {
        auto ref = comp;
        do {
            std::cout << ref << ", ";
            ref = ref.next();
        } while (ref != comp);
        std::cout << std::endl;
    }
}

std::vector<std::pair<regina::StrandRef,regina::StrandRef>> findLinkQuadriPairs(const regina::StrandRef &twoHandle) {
    std::vector<std::pair<regina::StrandRef,regina::StrandRef>> result;
    auto currentRef = twoHandle;

    if (isCurl(currentRef) && currentRef.next().next().crossing()->index() == currentRef.crossing()->index()) {
        // LT: here the component is just a single curl
        // if not dealt with separately, this process would find it has a quadricolour with itself
        // instead return empty vector
        return result;
    }

    do {
        auto next = currentRef.next();
        if (isCurl(currentRef)) {
            // The current crossing is a curl, and the next one is a curl of the same sign.
            if (isCurl(next)) {
                if (next.crossing()->index() == currentRef.crossing()->index()) {
                    auto next2 = next.next();
                    // LT: fixed this to compare crossing signs directly
                    // was previously comparing under/over strands, which fails when curls are same sign but different sides
                    if ((isCurl(next2)) && (next2.crossing()->sign() == currentRef.crossing()->sign())) {
                        result.emplace_back(currentRef,next2);
                    }
                }
            // The current crossing is a curl, and the next one is an undercrossing.
            } else if (next.strand() == 0) {
                result.emplace_back(currentRef,next);
            } 
        }
        else {
            // The current crossing is an undercrossing and the next one is a curl.
            if ((currentRef.strand() == 0) && isCurl(next)) {
                result.emplace_back(next,currentRef);
            } 
        }
        currentRef = currentRef.next();
    } while (currentRef != twoHandle);
    
    return result;
}

std::vector<int> pdCodeXTypes(const pdcode& code) {
    /*
     Assigns an identifier to each crossing of the link
     based on the PD code tuple, distinguishing between
     a "true" crossing and the four different PD code
     tuples that can arise from a curl.
     */
    std::vector<int> result;

    for (const auto& x : code) {
        if (x[2]==x[3]) {
            // (a,b,x,x) Positive
            result.emplace_back(1);
        }
        else if (x[0]==x[1]) {
            // (x,x,c,d) Positive
            result.emplace_back(2);
        }
        else if (x[1]==x[2]) {
            // (a,x,x,d) Negative
            result.emplace_back(3);
        }
        else if (x[0]==x[3]) {
            // (x,b,c,x) Negative
            result.emplace_back(4);
        }
        else {
            // regular crossing
            result.emplace_back(0);
        }
    }

    return result;
}

std::vector<int> pdCodeOrientations(const pdcode& code) {
    std::array<int, 4> eovInit = {0,0,0,0};
    
    std::array<int, 4> negative = {1,1,-1,-1};
    std::array<int, 4> positive = {1,-1,-1,1};

    
    long pdLength = code.size();
    int numberOfStrands = 2*pdLength;

    std::vector<std::array<int, 4>> extendedOrientationVector(pdLength,eovInit);
        
    std::vector<std::vector<bool>> visited(pdLength,std::vector<bool>(4,false));
    std::vector<int> seenStrands;
    
    int i = 0, j = 0;
    int currentStrand = code[i][j];
    int count = 1;
    
    while (!visited[i][j]) {
        
        bool carry = false;
        int carryRow = 0; // LT: weird behaviour when not initialised here, leading to segfault

        visited[i][j] = true;
        seenStrands.emplace_back(currentStrand);

        if (count%2 == 1) {
            extendedOrientationVector[i][j] = 1;
        }
        else {
            extendedOrientationVector[i][j] = -1;
        }
        count++;
        
        j = (j+2)%4;
        
        currentStrand = code[i][j];
        visited[i][j] = true;
        seenStrands.emplace_back(currentStrand);
        
        if (count%2 == 1) {
            extendedOrientationVector[i][j] = 1;
        }
        else {
            extendedOrientationVector[i][j] = -1;
        }
        count++;
        
        if (std::count(seenStrands.begin(),seenStrands.end(),currentStrand) == 2) {
            carry = true;
            for (int row=0; row<pdLength; row++) {
                if (!visited[row][0]) {
                    currentStrand = code[row][0];
                    carryRow = row;
                    break;
                }
            }
        }
        
        int nextI = -1, nextJ = -1;
        if (!carry) {
            for (int row=0; row<pdLength; row++) {
                for (int col=0; col<4; col++) {
                    if (!visited[row][col] && code[row][col] == currentStrand) {
                        nextI = row;
                        nextJ = col;
                        break;
                    }
                }
                if (nextI != -1) {
                    break;
                }
            }
        }
        else {
            nextI = carryRow;
            nextJ = 0;
        }
        
        i = nextI;
        j = nextJ;

        if (nextI == -1 && nextJ == -1) {
            break;
        }
    }

    std::vector<int> orientations;

    for (const auto& x : extendedOrientationVector) {
//      DEBUG: Print the current EOV tuple.
       // std::cout << x << std::endl;
        if (x == positive) {
            orientations.push_back(1);
        }
        else if (x == negative) {
            orientations.push_back(-1);
        }
    }

    return orientations;
}

std::vector<std::pair<int,int>> pdCodeXTypeOrientations(const pdcode& code) {
    /*
     Each element in this list is a pair consisting of:
        1.  The crossing type of the current crossing --
            "true" crossing or curl.
        2.  The orientation of the current crossing.
     */
    std::vector<std::pair<int,int>> result;

    std::vector<int> pdcXTypes = pdCodeXTypes(code);
    std::vector<int> pdcOrientations = pdCodeOrientations(code);
    
    for (size_t i=0; i<code.size(); i++) {
        result.emplace_back(pdcXTypes[i],pdcOrientations[i]);
    }
    
    return result;
}

template <int dimP1>
std::vector<std::tuple<int,int,int>> gluingList(graph<dimP1> G) {
    std::vector<std::tuple<int,int,int>> result;
    std::vector<node> nodes = G.nodes();
    std::vector<edge> edges = G.edges();
    result.reserve(edges.size());
    for (const auto& e : edges) {
        result.emplace_back(getIndex(nodes,e.n1),getIndex(nodes,e.n2),e.colour);
    }
    return result;
}

template <int dimP1>
void printGluingList(graph<dimP1> G) {
    std::vector<node> nodes = G.nodes();
    std::vector<edge> edges = G.edgesDirected();
    for (const auto& e : edges) {
        // Print commas after everything except the last line :D
        if (&e != &edges.back()) {
            std::cout << "[" << getIndex(nodes,e.n1) << ", " << getIndex(nodes,e.n2) << ", " << e.colour << "],\n";
        }
        else {
            std::cout << "[" << getIndex(nodes,e.n1) << ", " << getIndex(nodes,e.n2) << ", " << e.colour << "]\n";
        }
    }
}

graph<3> katie3Graph(std::vector<int> rawPDVect, std::vector<int> framingVector, int &nextIDGlobal) {

    /* LT: preprocessing for empty PD codes / circles 
     */
    if (rawPDVect.empty()) {
        // empty PD code means a circle, replace it with something we can work with
        // this a circle with four cancelling curls (two +ve in a row, two -ve in a row)
        // which is the smallest possible that won't need any extra pre-processing
        if (printDebugInfo) {std::clog << "Component is a circle, replacing with a circle with four cancelling curls\n\n";}
        rawPDVect = {1,8,2,1,3,2,4,3,4,6,5,5,6,8,7,7};
        framingVector = {0};
    }

    int nextID = 0;
    
    pdcode PDCworking;
        
    /*
     START Process PD Code
     */
    
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
    /*
     END Process PD Code
     */
    
    /*
     START Process Framings
     */

    std::vector<int> twoHandleFramings;
    for (int i=0; i<framingVector.size(); ++i) {
        twoHandleFramings.push_back(framingVector[i]);
    }
    /*
     END Process Framings
     */

    regina::Link linkObjWorking = regina::Link::fromPD(PDCworking.begin(),PDCworking.end());
    
    size_t numberOfComponents = linkObjWorking.countComponents();
        
    /*
     Dedicated vectors containing references to 1- and 2-handles.
     WARNING/NOTE:  These vectors will be of size "numberOfOneHandles"
                    and "numberOfTwoHandles" respectively.
                    So indexing these vectors is done w.r.t these sizes as well.
                    This could be a potential "vector" for "mismatched index"
                    errors later on down the track, so keep these ones in mind.
     */
    std::vector<regina::StrandRef> twoHandleComponentRefs;
    for (int i=0; i<numberOfComponents; i++) {
        twoHandleComponentRefs.emplace_back(linkObjWorking.component(i));
    }
    
    size_t numberOfTwoHandles;
    numberOfTwoHandles = twoHandleComponentRefs.size();
    
    /*
     Dedicated vectors containing the crossing indices of the 1- and 2-handles.
     */
    std::vector<std::set<int>> twoHandleCrossingIndices;
    for (const auto& twoHandleRef : twoHandleComponentRefs) {
        std::set<int> currentTwoHandleCrossingIndices;
        auto currentTwoHandleRef = twoHandleRef;
        do {
            currentTwoHandleCrossingIndices.emplace(currentTwoHandleRef.crossing()->index());
            currentTwoHandleRef = currentTwoHandleRef.next();
        } while (currentTwoHandleRef != twoHandleRef);
        twoHandleCrossingIndices.emplace_back(currentTwoHandleCrossingIndices);
    }
    
    // init debugging
    if (printDebugInfo) {
        std::clog << "PD code:  " << linkObjWorking.pd() << std::endl;
        std::clog << "Framings: ";
        for (const auto& x : framingVector) {
            std::clog << x << ", ";
        }
        std::clog << std::endl << std::endl;
        
        std::clog << "There are " << numberOfTwoHandles << " 2-handles." << std::endl;
 
        std::clog << "2-handle crossing indices:\n";
        for (const auto& x : twoHandleCrossingIndices) {
            for (const auto& y : x) {
                std::clog << y << ", ";
            }
            std::clog << std::endl;
        }

        std::clog << "\n";
    }
    // end init debugging
    
    /*
     START Framing Procedure
     */
    std::vector<long> twoHandleWrithes;
    for (const auto& twoHandle : twoHandleComponentRefs) {
        twoHandleWrithes.emplace_back(linkObjWorking.writheOfComponent(twoHandle));
    }
    
    for (int i=0; i<numberOfTwoHandles; i++) {
        long currentWrithe = twoHandleWrithes[i];
        int currentFraming = twoHandleFramings[i];
        regina::StrandRef twoHandle = twoHandleComponentRefs[i];
        
        if (currentWrithe > currentFraming) {
            if (printDebugInfo) {
                std::clog << "Self-framing component " << i << " (--)\n";
            }
            do {
                linkObjWorking.r1(twoHandle, 0 /* left */, -1);
                --currentWrithe;
            } while (currentWrithe != currentFraming);
        }
        else if (currentWrithe < currentFraming) {
            if (printDebugInfo) {
                std::clog << "Self-framing component " << i << " (++)\n";
            }
            do {
                linkObjWorking.r1(twoHandle, 0 /* left */, 1);
                ++currentWrithe;
            } while (currentWrithe != currentFraming);
        } else if (printDebugInfo) {std::clog << "Component " << i << " is already self-framed\n";}
    }

    if (printDebugInfo) {std::clog << "\n";}

    /*
     LT: ensure all components have:
     (1) at least one undercrossing
         (no disconnected 0-coloured arcs)
     (2) at least one overcrossing
         (single connected 0-residue)
     (3) if it has >1 curl, at least one non-curl undercrossing or two curls of the same sign in a row
         (dipole moves to produce the curl subgraphs are valid)
     If not, add two cancelling curls or four cancelling curls (with two of the same sign in a row)
     as necessary.
     */
    for (int i=0; i<numberOfTwoHandles; ++i) {

        const regina::StrandRef& twoHandle = twoHandleComponentRefs[i];
        regina::StrandRef currentRef = twoHandle;
        bool hasOvercrossing = false;
        bool curlSimplificationValid = false;
        int countCurls = 0;
        regina::StrandRef lastCurl = twoHandle; // arbitrary starting assigment to avoid non-assignment issues

        do {
            if (currentRef.strand()==0) {
                // undercrossing
                if (!isCurl(currentRef)) {
                    // non-curl undercrossing
                    curlSimplificationValid = true;
                } else {
                    // curl
                    ++countCurls;
                    lastCurl = currentRef;
                    hasOvercrossing = true;
                    regina::StrandRef twoAhead = currentRef.next().next();
                    if (twoAhead.strand()==0 && isCurl(twoAhead)) {
                        // two curls of the same sign in a row
                        curlSimplificationValid = true;
                        break;
                    }
                }
            } else {
                // overcrossing
                hasOvercrossing = true;
            }
            if (hasOvercrossing && curlSimplificationValid) {
                break;
            }
            currentRef=currentRef.next();
        } while (currentRef != twoHandle);

        // note: if we broke from the loop early, countCurls may not be correct
        // but we only break if curl simplification is valid, so how many curls there are won't matter

        if (!curlSimplificationValid) {
            // note: 
            // - if the component has no curls, we still need to do this
            //   in that case, being here would imply it has no undercrossings, so we would have to add two curls
            //   and so would actually still need four to make the simplification valid
            // - if there is exactly one curl, we are in a special case where we don't need to do anything,
            //   since the first curl simplification is valid for free, and the curl is an over- and under-crossing
            if (countCurls == 0) {
                linkObjWorking.r1(twoHandle, 0, 1);
                linkObjWorking.r1(twoHandle, 0, -1);
                linkObjWorking.r1(twoHandle, 0, -1);
                linkObjWorking.r1(twoHandle, 0, 1);
                if (printDebugInfo) {
                    std::clog << "Component " << i << " does not have a non-curl undercrossing or any curls, adding four cancelling curls...\n";
                }
            } else if (countCurls > 1) {
                linkObjWorking.r1(lastCurl, 0, -lastCurl.crossing()->sign());
                linkObjWorking.r1(lastCurl, 0, lastCurl.crossing()->sign());
                if (printDebugInfo) {
                    std::clog << "Component " << i << " does not have a non-curl undercrossing or double curl, adding two cancelling curls...\n";
                }
            } else if (printDebugInfo) {std::clog << "Component " << i << " is good (has exactly one curl)\n";}
        } else if (!hasOvercrossing) {
            linkObjWorking.r1(twoHandle, 0, 1);
            linkObjWorking.r1(twoHandle, 0, -1);
            if (printDebugInfo) {
                std::clog << "Component " << i << " has only undercrossings, adding two cancelling curls...\n";
            }
        } else if (printDebugInfo) {std::clog << "Component " << i << " is good\n";}
    }
    
    /*
     END Framing Procedure
     */

    if (printDebugInfo) {std::clog << "\nPD code post framing and pre-processing:\n" << linkObjWorking.pd() << "\n";}
    
    // Debugging stage 2
    for (int i=0; i<numberOfComponents; i++) {
        long currentWrithe = linkObjWorking.writheOfComponent(i);
        if (currentWrithe != framingVector[i]) {
            throw std::logic_error("Something went wrong during the framing process");
        }
    }
    // end debugging
    
    /*
     Since the framing procedure changes the link,
     we need to recompute relevant link data so that
     indices, etc. match between the link object and
     graph objects generated later on.
     */
    regina::Link linkObjFramed = regina::Link(linkObjWorking.pd());
    linkObjWorking = linkObjFramed;
    
    // Recompute reference vectors
    twoHandleComponentRefs.clear();
    for (int i=0; i<numberOfComponents; i++) {
        twoHandleComponentRefs.emplace_back(linkObjWorking.component(i));
    }
    
    // Recompute crossing indices
    twoHandleCrossingIndices.clear();
    for (const auto& twoHandleRef : twoHandleComponentRefs) {
        std::set<int> currentTwoHandleCrossingIndices;
        auto currentTwoHandleRef = twoHandleRef;
        do {
            currentTwoHandleCrossingIndices.emplace(currentTwoHandleRef.crossing()->index());
            currentTwoHandleRef = currentTwoHandleRef.next();
        } while (currentTwoHandleRef != twoHandleRef);
        twoHandleCrossingIndices.emplace_back(currentTwoHandleCrossingIndices);
    }
    
    if (printDebugInfo) {
        std::clog << std::endl;
        std::clog << "Post-recompute walk-around:\n";
        walkAroundLink(linkObjWorking);
        std::clog << std::endl;
    }
    
    graph<3> posCross, negCross, posCurlA, posCurlB, negCurlA, negCurlB;
    graph<3> boundaryGraph;

    node n1  = { 1,0}, n2  = { 2,0}, n3  = { 3,0}, n4  = { 4,0};
    node n5  = { 5,0}, n6  = { 6,0}, n7  = { 7,0}, n8  = { 8,0};
    node n9  = { 9,1}, n10 = {10,1}, n11 = {11,1}, n12 = {12,1};
    node n13 = {13,2}, n14 = {14,2}, n15 = {15,2}, n16 = {16,2};
    node n17 = {17,3}, n18 = {18,3}, n19 = {19,3}, n20 = {20,3};
    node n21 = {21,4}, n22 = {22,4}, n23 = {23,4}, n24 = {24,4};

    node pca5 = {5,1}, pca6  = { 6,1}, pca7  = { 7,1}, pca8  = { 8,1};
    node pca9 = {9,2}, pca10 = {10,2}, pca11 = {11,2}, pca12 = {12,2};

    node pcb5 = {5,4}, pcb6  = { 6,4}, pcb7  = { 7,4}, pcb8  = { 8,4};
    node pcb9 = {9,3}, pcb10 = {10,3}, pcb11 = {11,3}, pcb12 = {12,3};

    node nca5 = {5,1}, nca6  = { 6,1}, nca7  = { 7,1}, nca8  = { 8,1};
    node nca9 = {9,4}, nca10 = {10,4}, nca11 = {11,4}, nca12 = {12,4};

    node ncb5 = {5,2}, ncb6  = { 6,2}, ncb7  = { 7,2}, ncb8  = { 8,2};
    node ncb9 = {9,3}, ncb10 = {10,3}, ncb11 = {11,3}, ncb12 = {12,3};
    
    std::vector<edge> posCrossEdgeList = {
        {n1, n6, 0},    {n1,n16,1}, {n1,n8,2},  {n1,n2,3},
        {n2, n5, 0},    {n2,n13,1}, {n2,n3,2},
        {n3, n11,0},    {n3,n12,1},             {n3,n8,3},
        {n4, n10,0},    {n4,n9, 1}, {n4,n5,2},  {n4,n7,3},
                        {n5,n24,1},             {n5,n6,3},
                        {n6,n21,1}, {n6,n7,2},
        {n7, n19,0},    {n7,n20,1},
        {n8, n18,0},    {n8,n17,1},
        {n14,n23,0},
        {n15,n22,0}
    };
    
    std::vector<edge> negCrossEdgeList = {
        {n1, n6, 0},    {n1,n24,1}, {n1,n8,2},  {n1,n2,3},
        {n2, n5, 0},    {n2,n21,1}, {n2,n3,2},
        {n3, n19,0},    {n3,n20,1},             {n3,n8,3},
        {n4, n18,0},    {n4,n17,1}, {n4,n5,2},  {n4,n7,3},
                        {n5,n16,1},             {n5,n6,3},
                        {n6,n13,1}, {n6,n7,2},
        {n7, n11,0},    {n7,n12,1},
        {n8, n10,0},    {n8,n9, 1},
        {n14,n23,0},
        {n15,n22,0}
    };
    
    std::vector<edge> posCurlAEdgeList = {
        {n1,pca6, 0},   {n1,pca9, 1},   {n1,n2,2},  {n1,n4,3},
        {n2,pca7, 0},   {n2,pca8, 1},               {n2,n3,3},
        {n3,pca10,0},   {n3,pca5, 1},   {n3,n4,2},
        {n4,pca11,0},   {n4,pca12,1}
    };
    
    std::vector<edge> posCurlBEdgeList = {
        {n1,pcb6, 0},   {n1,pcb9, 1},   {n1,n2,2},  {n1,n4,3},
        {n2,pcb7, 0},   {n2,pcb8, 1},               {n2,n3,3},
        {n3,pcb10,0},   {n3,pcb5, 1},   {n3,n4,2},
        {n4,pcb11,0},   {n4,pcb12,1}
    };
    
    std::vector<edge> negCurlAEdgeList = {
        {n1,nca6, 0},   {n1,nca5, 1},   {n1,n2,2},  {n1,n4,3},
        {n2,nca7, 0},   {n2,nca12,1},               {n2,n3,3},
        {n3,nca10,0},   {n3,nca9, 1},   {n3,n4,2},
        {n4,nca11,0},   {n4,nca8, 1}
    };
    
    std::vector<edge> negCurlBEdgeList = {
        {n1,ncb6, 0},   {n1,ncb5, 1},   {n1,n2,2},  {n1,n4,3},
        {n2,ncb7, 0},   {n2,ncb12,1},               {n2,n3,3},
        {n3,ncb10,0},   {n3,ncb9, 1},   {n3,n4,2},
        {n4,ncb11,0},   {n4,ncb8, 1}
    };
    
    posCross.addEdges(posCrossEdgeList);
    negCross.addEdges(negCrossEdgeList);
    posCurlA.addEdges(posCurlAEdgeList);
    posCurlB.addEdges(posCurlBEdgeList);
    negCurlA.addEdges(negCurlAEdgeList);
    negCurlB.addEdges(negCurlBEdgeList);
    
    pdcode pdCodeMain = linkObjWorking.pdData();
    
    std::vector<std::pair<int,int>> pdcXOTypes = pdCodeXTypeOrientations(pdCodeMain);

    /* 
     LT: the way the nodes have been labelled, n1 will unfortunately be in different
     bipartite classes in crossings and in curls
     set the orientation node based on which of these the first subgraph encountered is
     (there is a binary choice here, this is the right one so that orientation agrees
     with how regina calculates intersection forms)
     */
    if (pdcXOTypes[0].first == 0) {
        boundaryGraph.setOrientation({2,0,nextID});
    } else {
        boundaryGraph.setOrientation({1,0,nextID});
    }
    
    long totalCrossingCounter = 1;
    for (const auto& pair : pdcXOTypes) {
        ++totalCrossingCounter;
        if ((pair.first == 0) && (pair.second == 1)) {
            if (printDebugInfo) {
                std::clog << "Building positive crossing graph...\n";
            }
            boundaryGraph.disjoint_union(posCross,nextID);
        }
        else if ((pair.first == 0) && (pair.second == -1)) {
            if (printDebugInfo) {
                std::clog << "Building positive negative graph...\n";
            }
            boundaryGraph.disjoint_union(negCross,nextID);
        }
        else if (pair.first == 1) {
            if (printDebugInfo) {
                std::clog << "Building positive curl (type A) graph...\n";
            }
            boundaryGraph.disjoint_union(posCurlA,nextID);
        }
        else if (pair.first == 2) {
            if (printDebugInfo) {
                std::clog << "Building positive curl (type B) graph...\n";
            }
            boundaryGraph.disjoint_union(posCurlB,nextID);
        }
        else if (pair.first == 3) {
            if (printDebugInfo) {
                std::clog << "Building negative curl (type A) graph...\n";
            }
            boundaryGraph.disjoint_union(negCurlA,nextID);
        }
        else if (pair.first == 4) {
            if (printDebugInfo) {
                std::clog << "Building negative curl (type B) graph...\n";
            }
            boundaryGraph.disjoint_union(negCurlB,nextID);
        }
    }
    
    /* LT: we need to not do this pdSub in the case where the whole link is just a single 
     circle with one curl.
    */
    if (pdCodeMain.size() > 1) {
        boundaryGraph.pdSub(pdCodeMain);
    }
    
    std::vector<std::pair<node,node>> bdryGfuseList = boundaryGraph.fuseList();

    for (auto pair : bdryGfuseList) {
        boundaryGraph.fuse(pair.first,pair.second);
    }

    boundaryGraph.cleanup();

    nextIDGlobal += nextID;

    return boundaryGraph;
    
}

graph<4> katie4Graph(std::vector<int> rawPDVect, std::vector<int> framingVector, std::vector<bool> isOneHandleVector, int &nextIDGlobal) {

    /* LT: preprocessing for empty PD codes / circles 
     */
    if (rawPDVect.empty()) {
        // empty PD code means a circle, replace with an equivalent PD code we can work with
        if (isOneHandleVector.empty() || !isOneHandleVector[0]) {
            // single 0-framed circle 2-handle
            // replace it with a circle with four cancelling curls (two +ve in a row, two -ve in a row)
            // this is the smallest possible which already has a quadricolour
            if (printDebugInfo) {std::clog << "Component is a circle 2-handle, replacing with a circle with four cancelling curls\n\n";}
            rawPDVect = {1,8,2,1,3,2,4,3,4,6,5,5,6,8,7,7};
            framingVector = {0};
            isOneHandleVector = {0};
        } else {
            // single 1-handle
            // replace it with a 1-handle, linked to a 2-handle, linked to another (cancelling) 1-handle
            // where the 2-handle is a circle with two cancelling curls, to create a quadricolour
            // to my knowledge, this is the smallest possible which already has a quadricolour + is in good position
            rawPDVect = {10,1,3,2,4,3,5,4,5,7,6,6,2,9,1,10,7,11,8,12,11,9,12,8};
            framingVector = {0,0,0};
            isOneHandleVector = {1,0,1};
        }
    }

    int nextID = 0;
    
    pdcode PDCworking;
        
    // /*
    //  START Process PD Code
    //  */
    
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
    /*
     END Process PD Code
     */
    
    /*
     START Process Framings
     */

    std::vector<int> twoHandleFramings;
    for (int i=0; i<framingVector.size(); ++i) {
        if (!isOneHandleVector[i]) {
            twoHandleFramings.push_back(framingVector[i]);
        }
    }
    /*
     END Process Framings
     */

    regina::Link linkObjWorking = regina::Link::fromPD(PDCworking.begin(),PDCworking.end());
    
    size_t numberOfComponents = linkObjWorking.countComponents();
    
    bool existOneHandles = std::any_of(isOneHandleVector.begin(),isOneHandleVector.end(),[](bool isOneHandle){return isOneHandle == true;});
    
    /*
     Dedicated vectors containing references to 1- and 2-handles.
     WARNING/NOTE:  These vectors will be of size "numberOfOneHandles"
                    and "numberOfTwoHandles" respectively.
                    So indexing these vectors is done w.r.t these sizes as well.
                    This could be a potential "vector" for "mismatched index"
                    errors later on down the track, so keep these ones in mind.
     */
    std::vector<regina::StrandRef> oneHandleComponentRefs, twoHandleComponentRefs;
    for (int i=0; i<numberOfComponents; i++) {
        if (isOneHandleVector[i]) {
            oneHandleComponentRefs.emplace_back(linkObjWorking.component(i));
        }
        else {
            twoHandleComponentRefs.emplace_back(linkObjWorking.component(i));
        }
    }
    
    size_t numberOfOneHandles, numberOfTwoHandles;
    numberOfOneHandles = oneHandleComponentRefs.size();
    numberOfTwoHandles = twoHandleComponentRefs.size();
    
    /*
     Dedicated vectors containing the crossing indices of the 1- and 2-handles.
     */
    std::vector<std::set<int>> oneHandleCrossingIndices, twoHandleCrossingIndices;
    for (const auto& oneHandleRef : oneHandleComponentRefs) {
        std::set<int> currentOneHandleCrossingIndices;
        auto currentOneHandleRef = oneHandleRef;
        do {
            currentOneHandleCrossingIndices.emplace(currentOneHandleRef.crossing()->index());
            currentOneHandleRef = currentOneHandleRef.next();
        } while (currentOneHandleRef != oneHandleRef);
        oneHandleCrossingIndices.emplace_back(currentOneHandleCrossingIndices);
    }
    for (const auto& twoHandleRef : twoHandleComponentRefs) {
        std::set<int> currentTwoHandleCrossingIndices;
        auto currentTwoHandleRef = twoHandleRef;
        do {
            currentTwoHandleCrossingIndices.emplace(currentTwoHandleRef.crossing()->index());
            currentTwoHandleRef = currentTwoHandleRef.next();
        } while (currentTwoHandleRef != twoHandleRef);
        twoHandleCrossingIndices.emplace_back(currentTwoHandleCrossingIndices);
    }
    
    /*
     Matrix consisting of StrandRefs for crossings of 2-handles
     which "intersect" 1-handles. That is, if we have a crossing like
     
                   | <-- 2-handle
                ---|---*----
                   |  /|\
                       |____ 1-handle
     
     then this matrix contains an entry for that crossing,
     indexed w.r.t the 2-handle(s).
     */
    std::vector<std::vector<regina::StrandRef>> oneTwoCommons;
    for (const auto& twoHandleRef : twoHandleComponentRefs) {
        std::vector<regina::StrandRef> currentCommons;
        auto currentTwoHandleRef = twoHandleRef;
        do {
            for (const auto& oneHandle : oneHandleCrossingIndices) {
                for (const auto& oneHandleCrossingIndx : oneHandle) {
                    if (currentTwoHandleRef.crossing()->index() == oneHandleCrossingIndx) {
                        currentCommons.emplace_back(currentTwoHandleRef);
                    }
                }
            }
            currentTwoHandleRef = currentTwoHandleRef.next();
        } while (currentTwoHandleRef != twoHandleRef);
        oneTwoCommons.emplace_back(currentCommons);
    }
    
    // init debugging
    if (printDebugInfo) {
        std::clog << "PD code:     " << linkObjWorking.pd() << std::endl;
        std::clog << "Framings:    ";
        for (const auto& x : framingVector) {
            std::clog << x << ", ";
        }
        std::clog << std::endl;
        std::clog << "One-handle?: ";
        for (const auto& x : isOneHandleVector) {
            std::clog << x << ", ";
        }
        std::clog << std::endl << std::endl;
        
        if (existOneHandles) {
            std::clog << "There are " << numberOfOneHandles << " 1-handles, and " << numberOfTwoHandles << " 2-handles." << std::endl;
        }
        else {
            std::clog << "There are no 1-handles, and " << numberOfTwoHandles << " 2-handles." << std::endl;
        }
        if (existOneHandles) {
            std::clog << "1-handle crossing indices:\n";
            for (const auto& x : oneHandleCrossingIndices) {
                for (const auto& y : x) {
                    std::clog << y << ", ";
                }
                std::clog << std::endl;
            }
        }
        std::clog << "2-handle crossing indices:\n";
        for (const auto& x : twoHandleCrossingIndices) {
            for (const auto& y : x) {
                std::clog << y << ", ";
            }
            std::clog << std::endl;
        }
        
        if (existOneHandles) {
            std::clog << "1/2 Commons:\n";
            for (const auto& x : oneTwoCommons) {
                for (const auto& y : x) {
                    std::clog << y << ", ";
                }
                std::clog << std::endl;
            }
        }

        std::clog << "\n";
    }
    // end init debugging
    
    /*
     START Framing Procedure
     */

    std::vector<long> oneHandleWrithes, twoHandleWrithes;
    for (const auto& twoHandle : twoHandleComponentRefs) {
        twoHandleWrithes.emplace_back(linkObjWorking.writheOfComponent(twoHandle));
    }
    for (const auto& oneHandle : oneHandleComponentRefs) {
        oneHandleWrithes.emplace_back(linkObjWorking.writheOfComponent(oneHandle));
    }
    /*
     While we're at it, check the writhes of any 1-handles.
     If they aren't 0, then this could indicate that the user
     has drawn the 1-handle in a "non-standard" way (i.e. as
     not a proper unknot), so we should alert the user and bail.
     */
    std::vector<int> badOneHandleComponentIndices;
    for (int i=0; i<numberOfOneHandles; i++) {
        auto currentOneHandleWrithe = oneHandleWrithes[i];
        if (currentOneHandleWrithe != 0) {
            badOneHandleComponentIndices.emplace_back(i);
        }
    }
    if (!badOneHandleComponentIndices.empty()) {
        throw std::invalid_argument("Input contains 1-handles with nonzero writhes");
    }

    // LT: record 2-handles for which we have already guaranteed existence of a quadricolour
    std::set<regina::StrandRef> guaranteedQuadricolours;
    // LT: record 2-handles for which we have already guaranteed at least two curls or one overcrossing
    std::set<regina::StrandRef> guaranteedDipoleValid;

    std::vector<regina::StrandRef> r1FramingSites;
    for (int i=0; i<numberOfTwoHandles; i++) {
        auto currentTwoHandle = twoHandleComponentRefs[i];
        auto currentCommons = oneTwoCommons[i];
        /*
         ^^^ Shouldn't be any "mismatched index" type errors here,
         since twoHandleComponentRefs is obviously indexed w.r.t 2-handles,
         but so is oneTwoCommons. So indexing them like this should be fine.
         */
        if (twoHandleWrithes[i] == twoHandleFramings[i]) {
            // LT: this 2-handle is already self-framed, site doesn't matter
            r1FramingSites.emplace_back(currentTwoHandle);
            continue;
        }
        bool found = false;
        for (const auto& commonCrossing : currentCommons) {
            /*
             See if this 2-handle intersects any 1-handles.
             If it does, make sure that the next crossing as we
             travel along the 2-handle is the next common 1-handle
             intersection crossing. Stick the R1 curl between these two:
                            
                        | <-- current 2-handle
                 ...----|----*---... <-- 1-handle (under 2-handle strand)
                        |
                        | <-- Stick R1 curls here.
                        |
                 ...---------*---... <-- same 1-handle (over 2-handle strand)
                        |
             
             */
            if (std::find(currentCommons.begin(),currentCommons.end(),commonCrossing.next() ) != currentCommons.end() ) {
                r1FramingSites.emplace_back(commonCrossing);
                found = true;
                guaranteedQuadricolours.insert(currentTwoHandle);
                break;
            }
            /*
             LT: I don't think any 2-handle intersecting a 1-handle necessarily intersects it under/over consecutively like this.
             Nonetheless, this seems to just be an opimisation trick
             TODO: look more closely at this + construction of oneTwoCommons
            */
        }
        if (found) {
            continue;
        }

        /*
         LT: if framing is off by at least two, we're guaranteed two
         curls of the same type in a row. Just put them at the start.
        */
        if (twoHandleWrithes[i] > twoHandleFramings[i]+1 || twoHandleWrithes[i] < twoHandleFramings[i]-1) {
            r1FramingSites.emplace_back(currentTwoHandle);
            found = true;
            guaranteedQuadricolours.insert(currentTwoHandle);
            guaranteedDipoleValid.insert(currentTwoHandle);
            continue;
        }
        /*
         LT:
         Here, we know we will add exactly one curl to self-frame.
         Try to place the r1 site directly after either an
         undercrossing or a curl of the same sign as the one we
         will be adding, to guarantee a quadricolour.
        */
        auto currentRef = currentTwoHandle;
        int signToAdd = twoHandleWrithes[i] > twoHandleFramings[i] ? -1 : 1;
        do {
            if (!isCurl(currentRef)) {
                if (currentRef.strand()==0) {
                    // found a non-curl undercrossing
                    r1FramingSites.emplace_back(currentRef);
                    found = true;
                    guaranteedQuadricolours.insert(currentTwoHandle);
                    break;
                }
            } else {
                if (currentRef.crossing()->sign() == signToAdd) {
                    // found a curl of the required sign
                    r1FramingSites.emplace_back(currentRef);
                    found = true;
                    guaranteedQuadricolours.insert(currentTwoHandle);
                    break;
                }
            }
            currentRef=currentRef.next();
        } while (currentRef != currentTwoHandle);

        if (found) {
            continue;
        }

        /*
         Lucy:
         There are no (non-curl) undercrossings or curls of the 
         correct sign, just put it at the start and check for a 
         quadricolour later. 
        */ 
        r1FramingSites.emplace_back(currentTwoHandle);

    }
    
    for (int i=0; i<numberOfTwoHandles; i++) {
        long currentWrithe = twoHandleWrithes[i];
        int currentFraming = twoHandleFramings[i];
        regina::StrandRef r1FramingSiteRef = r1FramingSites[i];
        
        if (currentWrithe > currentFraming) {
            if (printDebugInfo) {
                std::clog << "Self-framing 2-handle " << i << " (--" << currentWrithe-currentFraming << ")\n";
            }
            do {
                linkObjWorking.r1(r1FramingSiteRef, 0 /* left */, -1);
                --currentWrithe;
            } while (currentWrithe != currentFraming);
        }
        else if (currentWrithe < currentFraming) {
            if (printDebugInfo) {
                std::clog << "Self-framing 2-handle " << i << " (++" << currentFraming-currentWrithe << ")\n";
            }
            do {
                linkObjWorking.r1(r1FramingSiteRef, 0 /* left */, 1);
                ++currentWrithe;
            } while (currentWrithe != currentFraming);
        } else if (printDebugInfo) {std::clog << "2-handle " << i << " is already self-framed\n";}
    }

    if (printDebugInfo) {std::clog << "\n";}
    
     /*
     LT: modified (hopefully better) procedure for ensuring "link quadricolours"
     (a curl adjacent to an undercrossing or curl of the same sign)
      For any 2-handles without link quadricolours we search for either:
     (a) a curl, or
     (b) a non-curl undercrossing
     and place a pair of cancelling curls after it to guarantee a
     quadricolour. In case (a), we must make sure the curl directly
     after matches the sign of the existing curl.
     If neither case occurs (i.e. component has only over-crossings),
     just add four cancelling curls to get two of the same sign in
     a row.
    */
    for (int i=0; i<numberOfTwoHandles; ++i) {
        const auto& twoHandle = twoHandleComponentRefs[i];
        if (!guaranteedQuadricolours.contains(twoHandle)) {
            std::vector<std::pair<regina::StrandRef,regina::StrandRef>> tempQuadriCheck = findLinkQuadriPairs(twoHandle);
            if (tempQuadriCheck.empty()) {
                if (printDebugInfo) {std::clog << "Adding extra curls to 2-handle " << i << " to create a quadricolour...\n";}
                auto currentRef = twoHandle;
                regina::StrandRef quadriSite;
                bool foundQuadriSite = false;
                int curlSign = 0;
                do {
                    if (isCurl(currentRef)) {
                        if (currentRef.crossing()->index() != currentRef.next().crossing()->index()) {
                            // second occurence of a curl
                            quadriSite = currentRef;
                            foundQuadriSite = true;
                            curlSign = currentRef.crossing()->sign();
                            break;
                        }
                    } else {
                        if (currentRef.strand()==0) {
                            // non-curl undercrossing
                            quadriSite = currentRef;
                            foundQuadriSite = true;
                            break;
                        }
                    }
                    currentRef=currentRef.next();
                } while (currentRef != twoHandle);
                if (foundQuadriSite) {
                    if (curlSign == 0) {
                        linkObjWorking.r1(quadriSite, 0, 1);
                        linkObjWorking.r1(quadriSite, 0, -1);
                    } else {
                        linkObjWorking.r1(quadriSite, 0, -curlSign); // note: these always get added at quadriSite, so they appear in reverse order
                        linkObjWorking.r1(quadriSite, 0, curlSign);
                    }
                    if (printDebugInfo) {std::clog << "  added two cancelling curls\n";}
                } else {
                    // component has no curls or non-curl under-strands
                    // just add four extra curls to guarantee two of a type in a row
                    linkObjWorking.r1(twoHandle, 0, 1);
                    linkObjWorking.r1(twoHandle, 0, 1);
                    linkObjWorking.r1(twoHandle, 0, -1);
                    linkObjWorking.r1(twoHandle, 0, -1);
                    if (printDebugInfo) {std::clog << "  added four cancelling curls\n";} 
                }
            } else if (printDebugInfo) {std::clog << "2-handle " << i << " already has a quadricolour\n";}
        } 
    }

    /*
     LT check if there are any components consisting of only 
     undercrossings + one curl.
     (this may invalidate an implicit dipole cancellation)
     If there are, add two cancelling curls at the beginning.
    */
    for (int i=0; i<numberOfTwoHandles; ++i) {
        const auto& twoHandle = twoHandleComponentRefs[i];
        if (!guaranteedDipoleValid.contains(twoHandle)) {
            auto currentRef = twoHandle;
            std::optional<size_t> firstCurlIndex;
            bool valid = false;
            do {
                if (isCurl(currentRef)) {
                    if (firstCurlIndex.has_value()) {
                        if (currentRef.crossing()->index() != firstCurlIndex.value()) {
                            valid = true;
                            break;
                        }
                    } else {
                        firstCurlIndex = currentRef.crossing()->index();
                    }
                } else if (currentRef.strand()==1) {
                    valid = true;
                    break;
                }
                currentRef=currentRef.next();
            } while (currentRef != twoHandle);
            if (!valid) {
                linkObjWorking.r1(twoHandle, 0, 1);
                linkObjWorking.r1(twoHandle, 0, -1);
                if (printDebugInfo) {std::clog << "Component " << i << " has no overcrossings and <2 curls, adding two cancelling curls...\n";}
            }
        }
    }

    if (printDebugInfo) {
        int i=0;
        for (const auto& twoHandle : twoHandleComponentRefs) {
            std::vector<std::pair<regina::StrandRef,regina::StrandRef>> tempQuadriCheck = findLinkQuadriPairs(twoHandle);
            if (tempQuadriCheck.empty()) {
                std::clog << "ERROR: component " << i << " still does not have a quadricolour\n";
            }
            ++i;
        }
        std::clog << "\nPD code post framing and pre-processing:\n" << linkObjWorking.pd() << "\n";
    }

    /*
     END Framing Procedure
     */
    
    // Debugging stage 2
    // std::clog << "Writhes:\n";
    for (int i=0; i<numberOfComponents; i++) {
        // std::clog << "Component " << i << ": ";
        long currentWrithe = linkObjWorking.writheOfComponent(i);
        if (isOneHandleVector[i]) {
            // std::clog << "1-handle (" << currentWrithe << ")\n";
        }
        else {
            // std::clog << "2-handle, writhe " << currentWrithe;
            if (currentWrithe != framingVector[i]) {
                throw std::logic_error("Something went wrong during the framing process");
            }
            else {
                // std::clog << std::endl;
            }
        }
    }
    // end debugging
    
    /*
     Since the framing procedure changes the link,
     we need to recompute relevant link data so that
     indices, etc. match between the link object and
     graph objects generated later on.
     */
    regina::Link linkObjFramed = regina::Link(linkObjWorking.pd());
    linkObjWorking = linkObjFramed;
    
    // Recompute reference vectors
    oneHandleComponentRefs.clear();
    twoHandleComponentRefs.clear();
    for (int i=0; i<numberOfComponents; i++) {
        if (isOneHandleVector[i]) {
            oneHandleComponentRefs.emplace_back(linkObjWorking.component(i));
        }
        else {
            twoHandleComponentRefs.emplace_back(linkObjWorking.component(i));
        }
    }
    
    // Recompute crossing indices
    oneHandleCrossingIndices.clear();
    twoHandleCrossingIndices.clear();
    for (const auto& oneHandleRef : oneHandleComponentRefs) {
        std::set<int> currentOneHandleCrossingIndices;
        auto currentOneHandleRef = oneHandleRef;
        do {
            currentOneHandleCrossingIndices.emplace(currentOneHandleRef.crossing()->index());
            currentOneHandleRef = currentOneHandleRef.next();
        } while (currentOneHandleRef != oneHandleRef);
        oneHandleCrossingIndices.emplace_back(currentOneHandleCrossingIndices);
    }
    for (const auto& twoHandleRef : twoHandleComponentRefs) {
        std::set<int> currentTwoHandleCrossingIndices;
        auto currentTwoHandleRef = twoHandleRef;
        do {
            currentTwoHandleCrossingIndices.emplace(currentTwoHandleRef.crossing()->index());
            currentTwoHandleRef = currentTwoHandleRef.next();
        } while (currentTwoHandleRef != twoHandleRef);
        twoHandleCrossingIndices.emplace_back(currentTwoHandleCrossingIndices);
    }
    
    // Recompute oneTwoCommons
    oneTwoCommons.clear();
    for (const auto& twoHandleRef : twoHandleComponentRefs) {
        std::vector<regina::StrandRef> currentCommons;
        auto currentTwoHandleRef = twoHandleRef;
        do {
            for (const auto& oneHandle : oneHandleCrossingIndices) {
                for (const auto& oneHandleCrossingIndx : oneHandle) {
                    if (currentTwoHandleRef.crossing()->index() == oneHandleCrossingIndx) {
                        currentCommons.emplace_back(currentTwoHandleRef);
                    }
                }
            }
            currentTwoHandleRef = currentTwoHandleRef.next();
        } while (currentTwoHandleRef != twoHandleRef);
        oneTwoCommons.emplace_back(currentCommons);
    }
    
    if (printDebugInfo) {
        std::clog << std::endl;
        std::clog << "Post-recompute walk-around:\n";
        walkAroundLink(linkObjWorking);
    }

    /*
     START 1-Handle Marked Crossings
     */
    /*
     Assuming 1-handle is traversed counter-clockwise:
     ⚪︎.first is the "leftmost" crossing, and
     ⚪︎.second is the "rightmost" crossing.
     
     TODO: Differentiate traversal directions. Until then insist 1-handles are drawn counter-clockwise, else "garbage in = garbage out".
     */
    std::vector<std::pair<regina::StrandRef,regina::StrandRef>> oneHandleMarkedCrossingRefs;
    for (const auto& oneHandle : oneHandleComponentRefs) {
        std::pair<regina::StrandRef,regina::StrandRef> currentPair;
        auto currentRef = oneHandle;
        // std::clog << "oneHandle " << currentRef.strand() << "\n";
        do {
            // std::clog << "strand " << currentRef.strand() << "\n";
            // std::clog << "next " << currentRef.next().strand() << "\n";
            if ((currentRef.strand() == 0) && (currentRef.next().strand() == 1)) {
                currentPair.first = currentRef;
            }
            if ((currentRef.strand() == 1) && (currentRef.next().strand() == 0)) {
                currentPair.second = currentRef.next();
            }
            currentRef = currentRef.next();
        } while (currentRef != oneHandle);
        oneHandleMarkedCrossingRefs.emplace_back(currentPair);
    }

    if (printDebugInfo) {
        if (existOneHandles) {
            std::clog << "1-Handle Marked Crossings:\n";
            for (const auto& pair : oneHandleMarkedCrossingRefs) {
                std::clog << pair.first << ", " << pair.second << std::endl;
            }
            std::clog << std::endl;
        }
    }

    /*
     END 1-Handle Marked Crossings
     */
    
    /*
     START Link Quadricolour Search
     */
    // LT: Make a list containing the quadricolour we want to use for each handle
    // previously this contained a loop which didn't actually do anything, so
    // now we just take the first quadricolour we find on each strand
    std::vector<std::pair<regina::StrandRef,regina::StrandRef>> quadriPairRefs(numberOfTwoHandles);
    for (int i=0; i<numberOfTwoHandles; i++) {
        auto twoHandle = twoHandleComponentRefs[i];
        std::vector<std::pair<regina::StrandRef,regina::StrandRef>> currentQuadricolourList = findLinkQuadriPairs(twoHandle);
        quadriPairRefs[i] = currentQuadricolourList.front();
    }
    if (printDebugInfo) {
        std::clog << "Quadricolour references:\n";
        for (const auto& pair : quadriPairRefs) {
            std::clog << pair.first << ", " << pair.second << "\n";
        }
        std::clog << std::endl;
    }
    /*
     END Link Quadricolour Search
     */
    
    /*
     START Highlighting Procedure
     */
    // LT: this now just highlights as much as it possibly can, i.e. the entire strand except
    // the undercrossing / previous curl used as part of the quadricolour
    // this leaves just one segment not highlighted on each 2-handle, and cuts everything open
    // previously, we sometimes ran into issues where the 1-handle marked points ended up
    // in the same region on only one side, and we picked the wrong side for them
    // this may look like overkill, but it's a trivial difference in time which saves
    // us doing more complicated checks to fix this
    std::vector<std::vector<regina::StrandRef>> highlightCrossings;
    std::vector<bool> walkOppDirVec;
    if (existOneHandles) {
        for (int i=0; i<numberOfTwoHandles; i++) {
            std::vector<regina::StrandRef> currHighlighted;

            bool walkOppositeDirection = false;
            
            auto currentTwoHandle = oneTwoCommons[i];
            if (!(currentTwoHandle.empty())) {
                
                auto currQuadri = quadriPairRefs[i];
                auto initRef = currQuadri.first;
                auto currQuadriX2 = currQuadri.second;
                
                if (isCurl(initRef)) {
                    if (initRef.next().next() == currQuadriX2) {
                        walkOppositeDirection = true;
                    }
                }
                if (initRef.next() == currQuadriX2) {
                    walkOppositeDirection = true;
                }
                
                auto walken = initRef;
                if (walkOppositeDirection) {
                    if (walken.prev().crossing()->index() == walken.crossing()->index()) {
                        walken = walken.prev().prev();
                    }
                    else {
                        walken = walken.prev();
                    }
                }

                if (walkOppositeDirection) {
                    do {
                        if (isCurl(walken)) {
                            currHighlighted.emplace_back(walken.prev());
                            walken = walken.prev().prev();
                        }
                        else {
                            currHighlighted.emplace_back(walken);
                            walken = walken.prev();
                        }
                    } while (!(walken == initRef || isCurl(walken) && walken.prev() == initRef));
                }
                else {
                    do {
                        if (isCurl(walken)) {
                            currHighlighted.emplace_back(walken);
                            walken = walken.next().next();
                        }
                        else {
                            currHighlighted.emplace_back(walken);
                            walken = walken.next();
                        }
                    }
                    while (!(walken == initRef || isCurl(walken) && walken.next() == initRef));
                }
            }
            
            highlightCrossings.emplace_back(currHighlighted);
        }
    }
    /*
     END Highlighting Procedure
     */
    if (printDebugInfo && existOneHandles) {
        std::clog << "Highlighted crossings:\n";
        for (const auto& twoHandle : highlightCrossings) {
            for (const auto& ref : twoHandle) {
                std::clog << ref << ", ";
            }
            std::clog << std::endl;
        }
        std::clog << std::endl;
    }
    
    graph<4> posCross, negCross, posCurlA, posCurlB, negCurlA, negCurlB;
    graph<4> boundaryGraph;

    node n1  = { 1,0}, n2  = { 2,0}, n3  = { 3,0}, n4  = { 4,0};
    node n5  = { 5,0}, n6  = { 6,0}, n7  = { 7,0}, n8  = { 8,0};
    node n9  = { 9,1}, n10 = {10,1}, n11 = {11,1}, n12 = {12,1};
    node n13 = {13,2}, n14 = {14,2}, n15 = {15,2}, n16 = {16,2};
    node n17 = {17,3}, n18 = {18,3}, n19 = {19,3}, n20 = {20,3};
    node n21 = {21,4}, n22 = {22,4}, n23 = {23,4}, n24 = {24,4};

    node pca5 = {5,1}, pca6  = { 6,1}, pca7  = { 7,1}, pca8  = { 8,1};
    node pca9 = {9,2}, pca10 = {10,2}, pca11 = {11,2}, pca12 = {12,2};

    node pcb5 = {5,4}, pcb6  = { 6,4}, pcb7  = { 7,4}, pcb8  = { 8,4};
    node pcb9 = {9,3}, pcb10 = {10,3}, pcb11 = {11,3}, pcb12 = {12,3};

    node nca5 = {5,1}, nca6  = { 6,1}, nca7  = { 7,1}, nca8  = { 8,1};
    node nca9 = {9,4}, nca10 = {10,4}, nca11 = {11,4}, nca12 = {12,4};

    node ncb5 = {5,2}, ncb6  = { 6,2}, ncb7  = { 7,2}, ncb8  = { 8,2};
    node ncb9 = {9,3}, ncb10 = {10,3}, ncb11 = {11,3}, ncb12 = {12,3};
    
    std::vector<edge> posCrossEdgeList = {
        {n1, n6, 0},    {n1,n16,1}, {n1,n8,2},  {n1,n2,3},
        {n2, n5, 0},    {n2,n13,1}, {n2,n3,2},
        {n3, n11,0},    {n3,n12,1},             {n3,n8,3},
        {n4, n10,0},    {n4,n9, 1}, {n4,n5,2},  {n4,n7,3},
                        {n5,n24,1},             {n5,n6,3},
                        {n6,n21,1}, {n6,n7,2},
        {n7, n19,0},    {n7,n20,1},
        {n8, n18,0},    {n8,n17,1},
        {n14,n23,0},
        {n15,n22,0}
    };
    
    std::vector<edge> negCrossEdgeList = {
        {n1, n6, 0},    {n1,n24,1}, {n1,n8,2},  {n1,n2,3},
        {n2, n5, 0},    {n2,n21,1}, {n2,n3,2},
        {n3, n19,0},    {n3,n20,1},             {n3,n8,3},
        {n4, n18,0},    {n4,n17,1}, {n4,n5,2},  {n4,n7,3},
                        {n5,n16,1},             {n5,n6,3},
                        {n6,n13,1}, {n6,n7,2},
        {n7, n11,0},    {n7,n12,1},
        {n8, n10,0},    {n8,n9, 1},
        {n14,n23,0},
        {n15,n22,0}
    };
    
    std::vector<edge> posCurlAEdgeList = {
        {n1,pca6, 0},   {n1,pca9, 1},   {n1,n2,2},  {n1,n4,3},
        {n2,pca7, 0},   {n2,pca8, 1},               {n2,n3,3},
        {n3,pca10,0},   {n3,pca5, 1},   {n3,n4,2},
        {n4,pca11,0},   {n4,pca12,1}
    };
    
    std::vector<edge> posCurlBEdgeList = {
        {n1,pcb6, 0},   {n1,pcb9, 1},   {n1,n2,2},  {n1,n4,3},
        {n2,pcb7, 0},   {n2,pcb8, 1},               {n2,n3,3},
        {n3,pcb10,0},   {n3,pcb5, 1},   {n3,n4,2},
        {n4,pcb11,0},   {n4,pcb12,1}
    };
    
    std::vector<edge> negCurlAEdgeList = {
        {n1,nca6, 0},   {n1,nca5, 1},   {n1,n2,2},  {n1,n4,3},
        {n2,nca7, 0},   {n2,nca12,1},               {n2,n3,3},
        {n3,nca10,0},   {n3,nca9, 1},   {n3,n4,2},
        {n4,nca11,0},   {n4,nca8, 1}
    };
    
    std::vector<edge> negCurlBEdgeList = {
        {n1,ncb6, 0},   {n1,ncb5, 1},   {n1,n2,2},  {n1,n4,3},
        {n2,ncb7, 0},   {n2,ncb12,1},               {n2,n3,3},
        {n3,ncb10,0},   {n3,ncb9, 1},   {n3,n4,2},
        {n4,ncb11,0},   {n4,ncb8, 1}
    };
    
    posCross.addEdges(posCrossEdgeList);
    negCross.addEdges(negCrossEdgeList);
    posCurlA.addEdges(posCurlAEdgeList);
    posCurlB.addEdges(posCurlBEdgeList);
    negCurlA.addEdges(negCurlAEdgeList);
    negCurlB.addEdges(negCurlBEdgeList);
    
    pdcode pdCodeMain = linkObjWorking.pdData();
    
    std::vector<std::pair<int,int>> pdcXOTypes = pdCodeXTypeOrientations(pdCodeMain);
    
    // LT: the way the nodes have been labelled, n1 will unfortunately be in different
    // bipartite classes in crossings and in curls
    // set the orientation node based on which of these the first subgraph encountered is
    // (there is a binary choice here, this is the right one so that orientation agrees
    // with how regina calculates intersection forms)
    if (pdcXOTypes[0].first == 0) {
        boundaryGraph.setOrientation({2,0,nextID});
    } else {
        boundaryGraph.setOrientation({1,0,nextID});
    }

    long totalCrossingCounter = 1;
    for (const auto& pair : pdcXOTypes) {
        ++totalCrossingCounter;
        if ((pair.first == 0) && (pair.second == 1)) {
            if (printDebugInfo) {
                std::clog << "Building positive crossing graph...\n";
            }
            boundaryGraph.disjoint_union(posCross,nextID);
        }
        else if ((pair.first == 0) && (pair.second == -1)) {
            if (printDebugInfo) {
                std::clog << "Building positive negative graph...\n";
            }
            boundaryGraph.disjoint_union(negCross,nextID);
        }
        else if (pair.first == 1) {
            if (printDebugInfo) {
                std::clog << "Building positive curl (type A) graph...\n";
            }
            boundaryGraph.disjoint_union(posCurlA,nextID);
        }
        else if (pair.first == 2) {
            if (printDebugInfo) {
                std::clog << "Building positive curl (type B) graph...\n";
            }
            boundaryGraph.disjoint_union(posCurlB,nextID);
        }
        else if (pair.first == 3) {
            if (printDebugInfo) {
                std::clog << "Building negative curl (type A) graph...\n";
            }
            boundaryGraph.disjoint_union(negCurlA,nextID);
        }
        else if (pair.first == 4) {
            if (printDebugInfo) {
                std::clog << "Building negative curl (type B) graph...\n";
            }
            boundaryGraph.disjoint_union(negCurlB,nextID);
        }
    }
    
    /* LT: note in the 3-manifold case, we need to not do this pdSub in the case where the whole
     link is just a single circle with one curl.
     In the 4-manifold case we can ignore this, since all 2-handles already have a quadricolour
     (curl-curl or curl-undercrossing) and 1-handles can't have curls.
    */
    boundaryGraph.pdSub(pdCodeMain);

    std::vector<std::pair<node,node>> bdryGfuseList = boundaryGraph.fuseList();

    for (auto pair : bdryGfuseList) {
        boundaryGraph.fuse(pair.first,pair.second);
    }

    boundaryGraph.cleanup();

    std::vector<std::array<node,4>> graphQuadriListAll, graphQuadriListFinal(numberOfTwoHandles);
    graphQuadriListAll = findGraphQuadricolours(boundaryGraph);
    if (printDebugInfo) {
        std::clog << std::endl;
        std::clog << "Graph quadricolours:\n";
        for (const auto& quadri : graphQuadriListAll) {
            std::clog << "nodes ";
            std::clog << getIndex(boundaryGraph.nodes(),quadri[0]) << ", " << getIndex(boundaryGraph.nodes(),quadri[1]) << ", " << getIndex(boundaryGraph.nodes(),quadri[2]) << ", " << getIndex(boundaryGraph.nodes(),quadri[3]);
            std::clog << " [nodes ";
            std::clog << quadri[0].nodeID << ", " << quadri[1].nodeID << ", " << quadri[2].nodeID << ", " << quadri[3].nodeID;
            std::clog << " in crossings ";
            std::clog << quadri[0].subgraphComponent << ", " << quadri[1].subgraphComponent << ", " << quadri[2].subgraphComponent << ", " << quadri[3].subgraphComponent << "]\n";
        }
        std::clog << std::endl;
    }
    
    /*
      LT: this now checks that the ORDERED pairs of crossings vs subgraph components match,
      where previously it was checking unordered pairs (as sets)
      this deals with some weird cases where an extra graph quadricolour shows up in addition
      to the one we want, e.g. in a "double twist"
            o---x
      o-----|-x |
      |   x-|-|-o
      | o-|-x |  
      | x-o   |  
      x-------o 
      using the wrong quadricolour in this situation messes up 1-handle highlighting
    */
    /*
      TODO (LT): should find a way to do away with findGraphQuadricolours() entirely, and pick
      out the graph quadricolour we want based on the linkQuadri we already know, since
      everything relies on us using this one specific graph quadricolour
    */
    for (int i=0; i<numberOfTwoHandles; i++) {
        auto linkQuadri = quadriPairRefs[i];
        std::pair<int,int> linkIndices = {linkQuadri.first.crossing()->index(),linkQuadri.second.crossing()->index()};
        for (const auto& graphQuadri : graphQuadriListAll) {
            if (  graphQuadri[0].subgraphComponent == linkIndices.first &&
                  graphQuadri[1].subgraphComponent == linkIndices.second &&
                  graphQuadri[2].subgraphComponent == linkIndices.first &&
                  graphQuadri[3].subgraphComponent == linkIndices.first) {
                graphQuadriListFinal[i] = graphQuadri;
                break;
            }
        }
    }
    
    if (printDebugInfo) {
        std::clog << "Final graph quadricolours:\n";
        for (const auto& quadri : graphQuadriListFinal) {
            std::clog << "nodes ";
            std::clog << getIndex(boundaryGraph.nodes(),quadri[0]) << ", " << getIndex(boundaryGraph.nodes(),quadri[1]) << ", " << getIndex(boundaryGraph.nodes(),quadri[2]) << ", " << getIndex(boundaryGraph.nodes(),quadri[3]);
            std::clog << " (nodes ";
            std::clog << quadri[0].nodeID << ", " << quadri[1].nodeID << ", " << quadri[2].nodeID << ", " << quadri[3].nodeID;
            std::clog << " in crossings ";
            std::clog << quadri[0].subgraphComponent << ", " << quadri[1].subgraphComponent << ", " << quadri[2].subgraphComponent << ", " << quadri[3].subgraphComponent << ")\n";
        }
        std::clog << std::endl;
    }
    
//        boundaryGraph.DEBUGremainingNoCol4Nodes();
    boundaryGraph.addQuadriEdges(graphQuadriListFinal);

    if (existOneHandles) {
        
        std::vector<std::pair<node,node>> oneHandleMarkedNodes;
        for (const auto& pair : oneHandleMarkedCrossingRefs) {
            int leftComp = pair.first.crossing()->index();
            int rightComp = pair.second.crossing()->index();
            std::pair<node,node> currPair;
            
            int leftOrientation = pdcXOTypes[leftComp].second;
            int rightOrientation = pdcXOTypes[rightComp].second;
            
            if (leftOrientation == 1) {
                currPair.first = {7,0,leftComp};
            }
            if (leftOrientation == -1) {
                currPair.first = {3,0,leftComp};
            }
            if (rightOrientation == 1) {
                currPair.second = {4,0,rightComp};
            }
            if (rightOrientation == -1) {
                currPair.second = {8,0,rightComp};
            }
            oneHandleMarkedNodes.emplace_back(currPair);
        }

        boundaryGraph.addHighlightEdges(highlightCrossings);
        if (!boundaryGraph.checkOneHandleRhoThrees(oneHandleMarkedNodes)) {
            throw std::logic_error("1-handle marker edges are not valid rho_3 switches");
        }
        boundaryGraph.addOneHandleMarkerEdges(oneHandleMarkedNodes);
        boundaryGraph.addDoubleOneEdges();
        boundaryGraph.addRemainderEdges();
    }
    else {
        boundaryGraph.addDoubleOneEdges();
    }
    
    boundaryGraph.cleanup();

    nextIDGlobal += nextID;

    // TODO (LT): clean up naming, this is the 4-gem itself not the boundary graph
    return boundaryGraph;
    
}

inline graph<4> katie4Graph(std::vector<int> rawPDVect, std::vector<int> framingVector, std::vector<bool> isOneHandleVector) {
    int nextID = 0;
    return katie4Graph(rawPDVect,framingVector,isOneHandleVector,nextID);
}

template <int dim>
std::vector<int> bipartiteClasses(graph<dim> g) {
    std::vector<int> sign(g.nodes().size(),0);
    for (int i=0; i<g.nodes().size(); ++i) {
        // LT: this ~should~ always just set sign[0]=1 and break on the first loop
        // (i.e. the orientation node should be the 0th node)
        // but be extra careful in case ordering has gotten messed up somewhere 
        if (g.nodes()[i] == g.getOrientation()) {
            // std::clog << "orientation is node " << i << " / " << g.nodes()[i] << "\n";
            sign[i] = 1;
            break;
        }
    }
    int countDone = 1;
    while (countDone < g.nodes().size()) {
        for (const auto& [fromIndx,toIndx,facet] : gluingList(g)) {
            if (sign[fromIndx] != 0 && sign[toIndx] == 0) {
                sign[toIndx] = -sign[fromIndx];
                ++countDone;
            } else if (sign[fromIndx] == 0 && sign[toIndx] != 0) {
                sign[fromIndx] = -sign[toIndx];
                ++countDone;
            }
        }
    }
    return sign;
}

regina::Triangulation<3> katie::katie3(std::vector<std::tuple<std::vector<int>,std::vector<int>>> info, bool orient) {

    int nextID = 0;
    if (printDebugInfo) {std::clog << "This link has " << info.size() << " split components\n\n";}
    if (printDebugInfo) {std::clog << "--------------------\n\nBuilding graph for component 0 ...\n\n";}
    graph<3> gem = katie3Graph(std::get<0>(info[0]),std::get<1>(info[0]),nextID);
    int lastID = nextID;
    if (printDebugInfo) {std::clog << "\n... size " << gem.nodes().size() << "\n\n";}

    for (auto [rawPDVect,framingVector] : info | std::views::drop(1)) {
        if (printDebugInfo) {std::clog << "--------------------\n\nBuilding graph for next component ...\n\n";}
        graph<3> newGem = katie3Graph(rawPDVect,framingVector,nextID);
        if (printDebugInfo) {std::clog << "\n... size " << newGem.nodes().size() << "\n\n";} 
        if (printDebugInfo) {std::clog << "Connected sum with existing graph ...\n\n";}
        gem.connected_sum(newGem,false,lastID);
        lastID = nextID;
    }

    if (printDebugInfo) {std::clog << "--------------------\n\nFinished building graph, final size " << gem.nodes().size() << "\nBuilding triangulation ...\n";}
    std::vector<std::tuple<int,int,int>> gemGluingList = gluingList(gem);
    std::vector<int> nodeClass = bipartiteClasses(gem);
    regina::Triangulation<3> tri;
    regina::Perm<4> swap01(1,0,2,3);
    tri.newTetrahedra(gem.size());
    for (const auto& [fromIndx,toIndx,facet] : gemGluingList) {
        if (!orient) {
            tri.tetrahedron(fromIndx)->join(facet,tri.tetrahedron(toIndx),regina::Perm<4>());
        } else {
            if (nodeClass[fromIndx] == -1) {
                tri.tetrahedron(fromIndx)->join(swap01[facet],tri.tetrahedron(toIndx),swap01);
            } else {
                tri.tetrahedron(fromIndx)->join(facet,tri.tetrahedron(toIndx),swap01);
            }
        }
    }
    bool allIsWell = tri.isValid();
    if (!allIsWell) {
        std::clog << tri.isoSig() << "\n";
        std::clog << tri.source(regina::Language::Python) << "\n";
        throw std::logic_error("Something went unexpectedly wrong during the construction");
    }
    if (printDebugInfo) {std::clog << "... success!\n\n";}

    return tri;

}

regina::Triangulation<4> katie::katie4(std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> info, bool orient) {

    int nextID = 0;
    if (printDebugInfo) {std::clog << "This link has " << info.size() << " split components\n\n";}
    if (printDebugInfo) {std::clog << "--------------------\n\nBuilding graph for component 0 ...\n\n";}
    graph<4> gem = katie4Graph(std::get<0>(info[0]),std::get<1>(info[0]),std::get<2>(info[0]),nextID);
    int lastID = nextID;
    if (printDebugInfo) {std::clog << "\n... size " << gem.nodes().size() << "\n\n";}
    for (auto [rawPDVect,framingVector,isOneHandleVector] : info | std::views::drop(1)) {
        if (printDebugInfo) {std::clog << "--------------------\n\nBuilding graph for next component ...\n\n";}
        graph<4> newGem = katie4Graph(rawPDVect,framingVector,isOneHandleVector,nextID);
        if (printDebugInfo) {std::clog << "\n... size " << newGem.nodes().size() << "\n\n";} 
        if (printDebugInfo) {std::clog << "Connected sum with existing graph ...\n\n";}
        gem.connected_sum(newGem,false,lastID);
        lastID = nextID;
    }

    if (printDebugInfo) {std::clog << "--------------------\n\nFinished building graph, final size " << gem.nodes().size() << "\nBuilding triangulation ...\n";}
    std::vector<std::tuple<int,int,int>> gemGluingList = gluingList(gem);
    std::vector<int> nodeClass = bipartiteClasses(gem);
    regina::Triangulation<4> fourTri;
    regina::Perm<5> swap01(1,0,2,3,4);
    fourTri.newPentachora(gem.size());
    for (const auto& [fromIndx,toIndx,facet] : gemGluingList) {
        if (!orient) {
            fourTri.pentachoron(fromIndx)->join(facet,fourTri.pentachoron(toIndx),regina::Perm<5>());
        } else {
            if (nodeClass[fromIndx] == -1) {
                fourTri.pentachoron(fromIndx)->join(swap01[facet],fourTri.pentachoron(toIndx),swap01);
            } else {
                fourTri.pentachoron(fromIndx)->join(facet,fourTri.pentachoron(toIndx),swap01);
            }
        }
    }
    bool allIsWell = fourTri.isValid();
    if (!allIsWell) {
        std::clog << fourTri.isoSig() << "\n";
        std::clog << fourTri.source(regina::Language::Python) << "\n";
        throw std::logic_error("Something went unexpectedly wrong during the construction");
    }
    if (printDebugInfo) {std::clog << "... success!\n\n";}

    return fourTri;
}

/** 
 * LT: I have subsequently realised this doesn't make much sense topologically, each connected 
 * component of a Kirby diagram defines an ~oriented~ manifold and there is only ever one correct 
 * choice of how to do the connected sums to get the right manifold.
 * 
 * note: does orientation ever matter for one-handles and isolated 0-framed 2-handles?
 * if not, can reduce some unneccesary repetitions
**/
std::vector<regina::Triangulation<4>> katie::katie4AllSums(std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> info, bool orient) {

    if (printDebugInfo) {std::clog << "This link has " << info.size() << " split components\n";}
    int nextID = 0;
    std::vector<graph<4>> gems;
    unsigned long numSums = std::pow(2,info.size()-1);
    if (printDebugInfo) {std::clog << "There are " << numSums << " possible connected sums\n\n";}

    if (printDebugInfo) {std::clog << "Building graph for component 0 ...\n\n";}
    graph<4> firstGem = katie4Graph(std::get<0>(info[0]),std::get<1>(info[0]),std::get<2>(info[0]),nextID);
    int lastID = nextID;
    if (printDebugInfo) {std::clog << "\n... size " << firstGem.nodes().size();}

    if (printDebugInfo) {std::clog << "\n\nAdding to each connected sum choice:\n";}
    for (int j=0; j<numSums; ++j) {
        graph<4> gem;
        gem.fromAdjacencyList(firstGem.adjacencyList());
        gems.push_back(gem);
        if (printDebugInfo) {std::clog << "+" << std::flush;}
    }
    if (printDebugInfo) {std::clog << "\n\n";}

    for (int i=1; i<info.size(); ++i) {

        if (printDebugInfo) {std::clog << "Building graph for component " << i << " ...\n\n";}
        graph<4> nextGem = katie4Graph(std::get<0>(info[i]),std::get<1>(info[i]),std::get<2>(info[i]),nextID);
        if (printDebugInfo) {std::clog << "\n... size " << nextGem.nodes().size();}

        if (printDebugInfo) {std::clog << "\n\nAdding to each connected sum choice:\n";}
        for (int j=0; j<numSums; ++j) {
            std::bitset<32> orientations(j); // note: first element of a bitset is ones place, i.e. last entry in binary
            graph<4> gem;
            gem.fromAdjacencyList(nextGem.adjacencyList());
            gems[j].connected_sum(gem,orientations[i-1],lastID);
            if (printDebugInfo) {std::clog << (orientations[i-1] ? "-" : "+") << std::flush;}
        }
        lastID = nextID;
        if (printDebugInfo) {std::clog << "\n\n";}

    }

    if (printDebugInfo) {std::clog << "\n\nFinished building graphs, final size " << gems[0].nodes().size() << "\n";}

    std::vector<regina::Triangulation<4>> tris;
    for (int j=0; j<numSums; ++j) {
        if (printDebugInfo) {
            std::clog << "Building triangulation " << j << ", connected sum +";
            std::bitset<32> orientations(j);
            for (int i=0; i<info.size()-1; ++i) {
                std::clog << (orientations[i] ? "-" : "+");
            }
            std::clog << " ...\n";
        }
        std::vector<std::tuple<int,int,int>> gemGluingList = gluingList(gems[j]);
        std::vector<int> nodeClass = bipartiteClasses(gems[j]);
        regina::Triangulation<4> fourTri;
        regina::Perm<5> swap01(1,0,2,3,4);
        fourTri.newPentachora(gems[j].size());
        for (const auto& [fromIndx,toIndx,facet] : gemGluingList) {
            if (!orient) {
                fourTri.pentachoron(fromIndx)->join(facet,fourTri.pentachoron(toIndx),regina::Perm<5>());
            } else {
                if (nodeClass[fromIndx] == -1) {
                    fourTri.pentachoron(fromIndx)->join(swap01[facet],fourTri.pentachoron(toIndx),swap01);
                } else {
                    fourTri.pentachoron(fromIndx)->join(facet,fourTri.pentachoron(toIndx),swap01);
                }
            }
        }
        bool allIsWell = fourTri.isValid();
        if (!allIsWell) {
            std::clog << fourTri.isoSig() << "\n";
            std::clog << fourTri.source(regina::Language::Python) << "\n";
            throw std::logic_error("Something went unexpectedly wrong during the construction");
        }
        tris.push_back(fourTri);
        if (printDebugInfo) {std::clog << "... success!\n\n";}
    }

    return tris;
}

void katie::katie3PrintGraph(std::vector<std::tuple<std::vector<int>,std::vector<int>>> info) {

    int nextID = 0;
    if (printDebugInfo) {std::clog << "This link has " << info.size() << " split components\n\n";}
    if (printDebugInfo) {std::clog << "--------------------\n\nBuilding graph for component 0 ...\n\n";}
    graph<3> gem = katie3Graph(std::get<0>(info[0]),std::get<1>(info[0]),nextID);
    int lastID = nextID;
    if (printDebugInfo) {std::clog << "\n... size " << gem.nodes().size() << "\n\n";}

    for (auto [rawPDVect,framingVector] : info | std::views::drop(1)) {
        if (printDebugInfo) {std::clog << "--------------------\n\nBuilding graph for next component ...\n\n";}
        graph<3> newGem = katie3Graph(rawPDVect,framingVector,nextID);
        if (printDebugInfo) {std::clog << "\n... size " << newGem.nodes().size() << "\n\n";} 
        if (printDebugInfo) {std::clog << "Connected sum with existing graph ...\n\n";}
        gem.connected_sum(newGem,false,lastID);
        lastID = nextID;
    }

    if (printDebugInfo) {std::clog << "--------------------\n\nFinished building graph, final size " << gem.nodes().size() << "\nBuilding triangulation ...\n";}
    printGluingList(gem);

}

void katie::katie4PrintGraph(std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<bool>>> info) {

    int nextID = 0;
    if (printDebugInfo) {std::clog << "This link has " << info.size() << " split components\n\n";}
    if (printDebugInfo) {std::clog << "--------------------\n\nBuilding graph for component 0 ...\n\n";}
    graph<4> gem = katie4Graph(std::get<0>(info[0]),std::get<1>(info[0]),std::get<2>(info[0]),nextID);
    int lastID = nextID;
    if (printDebugInfo) {std::clog << "\n... size " << gem.nodes().size() << "\n\n";}
    for (auto [rawPDVect,framingVector,isOneHandleVector] : info | std::views::drop(1)) {
        if (printDebugInfo) {std::clog << "--------------------\n\nBuilding graph for next component ...\n\n";}
        graph<4> newGem = katie4Graph(rawPDVect,framingVector,isOneHandleVector,nextID);
        if (printDebugInfo) {std::clog << "\n... size " << newGem.nodes().size() << "\n\n";} 
        if (printDebugInfo) {std::clog << "Connected sum with existing graph ...\n\n";}
        gem.connected_sum(newGem,false,lastID);
        lastID = nextID;
    }

    if (printDebugInfo) {std::clog << "--------------------\n\nFinished building graph, final size " << gem.nodes().size() << "\nBuilding triangulation ...\n";}
    printGluingList(gem);

}