#include "include/solver.hpp"
#include "solver/include/solver.hpp"
#include <iostream>

int main(int argc, char** argv) {
    std::cout << "===Gross Pitaevski Solver===" << std::endl;

    GrossPitaevskiSolver solver;
    solver.solve();
}
