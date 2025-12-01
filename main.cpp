#include "include/solver.hpp"
#include "solver/include/solver.hpp"
#include "include/output.hpp"

int main(int argc, char** argv) {
    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Gross Pitaevski Equation Solver");
    OutputFormatter::printBorderLine();

    GrossPitaevskiSolver solver;
    solver.solve();
}
