#include "output.hpp"
#include "filemanager/file_manager.hpp"
#include "manager/sim_manager.hpp"
#include <cstring>
#include <iostream>

int main(int argc, char** argv) {
    if (argc == 2 && std::strcmp(argv[1], "--print_default_json") == 0) {
        FileManager::print_default_json(std::cout);
        return 0;
    }

    if (argc == 2 && std::strcmp(argv[1], "--print_default_binary_json") == 0) {
        FileManager::print_default_binary_json(std::cout);
        return 0;
    }

    if (argc > 1) {
        std::cerr << "Usage: " << argv[0]
                  << " [--print_default_json|--print_default_binary_json]" << std::endl;
        return 1;
    }

    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Gross Pitaevski Equation Solver");
    OutputFormatter::printBorderLine();

    std::srand(std::time(0));

    std::unique_ptr<SimulationManager> manager = std::make_unique<SimulationManager>();
    try{
        manager->initialize();
        manager->run_simulation();
    } catch (const std::exception& e){
        std::cerr << e.what() << std::endl;
    }
}
