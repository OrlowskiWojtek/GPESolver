#include "output.hpp"
#include "manager/sim_manager.hpp"

int main(int argc, char** argv) {
    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Gross Pitaevski Equation Solver");
    OutputFormatter::printBorderLine();

    std::unique_ptr<SimulationManager> manager = std::make_unique<SimulationManager>();
    try{
        manager->initialize();
        manager->run_simulation();
    } catch (const std::exception& e){
        std::cerr << e.what() << std::endl;
    }
}
