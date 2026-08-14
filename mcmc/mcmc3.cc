/**
 * Implementation for mcmc3.h
 * 
 * This is a C++ adaptation of code by Eduardo Altmann and Jonathan Spreer, available at
 * https://github.com/jspreer/MCMCForTriangulations
 **/

#include <cmath>
#include <iterator>

#include <triangulation/isosigtype.h>

#include "mcmc3.h"

using namespace mcmc;

std::set<std::string> mcmc::neighbours(std::string iso, bool up) {
    std::set<std::string> nbrs;
    const regina::Triangulation<3> base(iso);
    if (up) {
        for (int i=0; i<base.countTriangles(); ++i) {
            regina::Triangulation<3> target(base);
            if (target.pachner(target.triangle(i))) {
                std::string tiso = target.isoSig<regina::IsoSigDegrees<3,1>>();
                nbrs.insert(tiso);
            }
        }
    } else {
        for (int i=0; i<base.countEdges(); ++i) {
            regina::Triangulation<3> target(base);
            if (target.pachner(target.edge(i))) {
                std::string tiso = target.isoSig<regina::IsoSigDegrees<3,1>>();
                nbrs.insert(tiso);
            }
        }
    }
    return nbrs;
}

std::pair<std::string,std::vector<size_t>> mcmc::choosemove(std::string iso, std::vector<size_t> f, double gamma) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<> dist(0.0,1.0);

    bool up = (dist(mt) < std::exp((-1)*gamma*f[3]));
    if (!up && f[3] < 3) {
        // there are not enough tetrahedra for a 3-2 move, stay here
        return {iso,f};
    }
    std::set<std::string> nbrs = neighbours(iso,up);

    if (up) {
        if (dist(mt) > (double)nbrs.size()/(double)f[2]) {
            // stay here
            return {iso,f};
        } else {
            // go up (very likely): choose a random string from nbrs
            std::uniform_int_distribution<> choice(0,nbrs.size()-1);
            auto it = nbrs.begin();
            std::advance(it,choice(mt));
            return {*it, {f[0],f[1]+1,f[2]+2,f[3]+1}};
        }
    } else {
        if (dist(mt) > (double)nbrs.size()/(double)(f[2]-2)) {
            // stay here
            return {iso,f};
        } else {
            // go down (unlikely): choose a random string from nbrs
            std::uniform_int_distribution<> choice(0,nbrs.size()-1);
            auto it = nbrs.begin();
            std::advance(it,choice(mt));
            return {*it, {f[0],f[1]-1,f[2]-2,f[3]-1}};
        }
    }
}

std::string mcmc::simplify(std::string iso, int steps, int simplifyFreq, bool verbose) {

    regina::Triangulation<3> start(iso);
    std::vector<size_t> f = start.fVector();
    double gamma = 4.0/double(f[3]+10); // aiming to average roughly initial size + 10 during MCMC

    std::string best = iso;
    int best_size = f[3];

    int simplifyCounter = 1;
    for (int i=1; i<=steps; ++i) {
        std::tie(iso,f) = choosemove(iso,f,gamma);
        if (f[3] < best_size) {
            regina::Triangulation<3> trig(iso);
            if (trig.simplify()) {
                best = trig.isoSig<regina::IsoSigDegrees<3,1>>();
                best_size = trig.countTetrahedra();
            } else {
                best = iso;
                best_size = f[3];
            }
        }
        if (simplifyCounter == simplifyFreq) {
            if (verbose) {
                std::cout << i << " | n = " << f[3] << " | best so far = " << best_size << "\n";
            }
            regina::Triangulation<3> trig(iso);
            if (trig.simplify()) {
                if (verbose) {std::cout << "--> simplified to n = " << trig.countTetrahedra() << "\n";}
                if (trig.countTetrahedra() < best_size) {
                    best = trig.isoSig<regina::IsoSigDegrees<3,1>>();
                    best_size = trig.countTetrahedra();
                }
            }
            simplifyCounter = 0;
        }
        ++simplifyCounter;
    }

    return best;

}