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
 * Interactive demo program which generates a single random 3- or 4-manifold
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

int main() {
    
    std::string temp;

    int n;
    std::cout << "Enter a grid size:\n";
    std::cin >> n;

    bool knot;
    std::string knotString;
    std::cout << "All grids, or only knots? (enter 'link' or 'knot')\n";
    std::cin >> knotString;
    if (knotString == "knot") {
        knot = true;
    } else if (knotString == "link") {
        knot = false;
    } else {
        std::cerr << "Error: invalid input\n";
        return 1;
    }

    std::cout << "\nPress enter to generate a random " << (knot ? "grid knot" : "grid") <<  " of size " << n << "\n";
    std::cin.ignore();
    std::getline(std::cin,temp);
    auto start = std::chrono::high_resolution_clock::now();
    gecko::Grid grid;
    if (knot) {
        grid = gecko::randomKnot(n);
    } else { 
        grid = gecko::randomGrid(n);
    }
    auto now = std::chrono::high_resolution_clock::now();
    std::cout << "(took " << std::chrono::duration_cast<std::chrono::microseconds>(now-start).count() << " mu s)\n";
    std::cout << "X: ";
    gecko::printVector(grid.x);
    std::cout << "O: ";
    gecko::printVector(grid.o);
    gecko::drawGrid(grid);

    int dimInt;
    std::cout << "Dimension 3 or 4?\n";
    std::cin >> dimInt;
    if (!(dimInt==3 || dimInt==4)) {
        std::cerr << "Error: invalid input\n";
        return 1;
    }

    std::cout << "\n";

    if (dimInt==3) {

        std::cout << "Press enter to build the corresponding 3-manifold\n" << std::flush;
        std::cin.ignore();
        std::getline(std::cin,temp);

        std::cout << "Building ..." << std::flush;
        start = std::chrono::high_resolution_clock::now();
        regina::Triangulation<3> trig = gecko::build3Manifold(grid,false);
        now = std::chrono::high_resolution_clock::now();
        std::cout << " done (took " << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n\n";

        start = std::chrono::high_resolution_clock::now();
        std::string sig = trig.sig();
        now = std::chrono::high_resolution_clock::now();
        std::cout << "Isomorphism signature ... " << std::flush << sig 
                  << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n\n" << std::flush;

        start = std::chrono::high_resolution_clock::now();
        std::cout << "Approx hyperbolic volume ... "  << std::flush << regina::SnapPeaTriangulation(trig).volume();
        now = std::chrono::high_resolution_clock::now();
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n" << std::flush;

        start = std::chrono::high_resolution_clock::now();
        std::cout << "Homology ... H1 = "  << std::flush << trig.homology(1).str() 
                  << ", H2 = " << std::flush << trig.homology(2).str();
        now = std::chrono::high_resolution_clock::now();
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n" << std::flush;

        start = std::chrono::high_resolution_clock::now();
        std::cout << "Attempting to identify ... " << std::flush << gecko::identify3Manifold(sig);
        now = std::chrono::high_resolution_clock::now();
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n" << std::flush;

    } else if (dimInt==4) {

        std::cout << "Press enter to build the corresponding 4-manifold, with 2-handles only\n";
        std::cin.ignore();
        std::getline(std::cin,temp);

        std::cout << "Building ..." << std::flush;
        start = std::chrono::high_resolution_clock::now();
        regina::Triangulation<4> trig = gecko::build4ManifoldCanonical(grid,true,false);
        now = std::chrono::high_resolution_clock::now();
        std::cout << " done (took " << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n\n";

        start = std::chrono::high_resolution_clock::now();
        std::string sig = trig.sig();
        now = std::chrono::high_resolution_clock::now();
        std::cout << "Isomorphism signature ... " << std::flush << sig;
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n\n" << std::flush;

        start = std::chrono::high_resolution_clock::now();
        std::cout << "Homology ... H1 = "  << std::flush << trig.homology(1).str() 
                  << ", H2 = " << std::flush << trig.homology(2).str() 
                  << ", H3 = " << std::flush << trig.homology(3).str();
        now = std::chrono::high_resolution_clock::now();
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n" << std::flush;

        start = std::chrono::high_resolution_clock::now();
        std::cout << "Closed? (S3 boundary filled) ... " << std::flush << (trig.isClosed() ? "yes" : "no, ideal boundary");
        now = std::chrono::high_resolution_clock::now();
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n" << std::flush;

        start = std::chrono::high_resolution_clock::now();
        if (trig.isClosed()) {
            regina::IntersectionForm form = trig.intersectionForm();
            std::cout << "Intersection form ... " << std::flush << (form.even() ? "even" : "odd") << " rank " << form.rank() << " signature " << form.signature();
        } else {
            std::cout << "Attempting to identify boundary ... " << std::flush << gecko::identify3Manifold(trig.boundaryComponent(0)->build().sig());
        }
        now = std::chrono::high_resolution_clock::now();
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n" << std::flush;

        std::cout << "Simplifying ... " << std::flush;
        start = std::chrono::high_resolution_clock::now();
        trig.simplify();
        now = std::chrono::high_resolution_clock::now();
        std::cout << trig.sig();
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n\n" << std::flush; 


        std::cout << "Press enter to build the 4-manifold with canonical 1-handles\n";
        // std::cin.ignore();
        std::getline(std::cin,temp);

        std::cout << "Building ..." << std::flush;
        start = std::chrono::high_resolution_clock::now();
        trig = gecko::build4ManifoldCanonical(grid,false,false);
        now = std::chrono::high_resolution_clock::now();
        std::cout << " done (took " << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n\n";

        start = std::chrono::high_resolution_clock::now();
        sig = trig.sig();
        now = std::chrono::high_resolution_clock::now();
        std::cout << "Isomorphism signature ... " << std::flush << sig;
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n\n" << std::flush;

        start = std::chrono::high_resolution_clock::now();
        std::cout << "Homology ... H1 = "  << std::flush << trig.homology(1).str() 
                  << ", H2 = " << std::flush << trig.homology(2).str() 
                  << ", H3 = " << std::flush << trig.homology(3).str();
        now = std::chrono::high_resolution_clock::now();
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n" << std::flush;

        start = std::chrono::high_resolution_clock::now();
        std::cout << "Closed? (S3 boundary filled) ... " << std::flush << (trig.isClosed() ? "yes" : "no, ideal boundary");
        now = std::chrono::high_resolution_clock::now();
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n" << std::flush;

        start = std::chrono::high_resolution_clock::now();
        if (trig.isClosed()) {
            regina::IntersectionForm form = trig.intersectionForm();
            std::cout << "Intersection form ... " << std::flush << (form.even() ? "even" : "odd") << " rank " << form.rank() << " signature " << form.signature();
        } else {
            std::cout << "Attempting to identify boundary ... " << std::flush << gecko::identify3Manifold(trig.boundaryComponent(0)->build().sig());
        }
        now = std::chrono::high_resolution_clock::now();
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n" << std::flush;

        start = std::chrono::high_resolution_clock::now();
        std::cout << "Simplifying ... " << std::flush;
        trig.simplify();
        now = std::chrono::high_resolution_clock::now();
        std::cout << trig.sig();
        std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << " ms)\n" << std::flush;
    }

    return 0;
}