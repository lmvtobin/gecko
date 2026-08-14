/**
 * Simplifying 3-dimensional triangulations using Markov Chain Monte Carlo random walks
 * 
 * This is a direct C++ adaptation of code by Eduardo Altmann and Jonathan Spreer, with some
 * additional functionality applying this to simplification. The original is available at 
 * 
 * https://github.com/jspreer/MCMCForTriangulations
 * 
 * and is described in the paper
 * 
 * Altmann, E. G., & Spreer, J. (2026).
 * Sampling Triangulations of Manifolds Using Monte Carlo Methods.
 * Experimental Mathematics, 35(1), 260–274. 
 * https://doi.org/10.1080/10586458.2024.2433506
 **/

#ifndef __MCMC_3_H
#define __MCMC_3_H

#include <string>
#include <set>
#include <vector>
#include <random>

#include <triangulation/dim3.h>

namespace mcmc {

std::set<std::string> neighbours(std::string iso, bool up);

std::pair<std::string,std::vector<size_t>> choosemove(std::string iso, std::vector<size_t> f, double gamma);

std::string simplify(std::string iso, int steps=10000, int simplifyFreq=100, bool verbose=false);

}

#endif