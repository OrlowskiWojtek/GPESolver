#ifndef FILE_MANAGER_HPP
#define FILE_MANAGER_HPP

#include "context/context.hpp"
#include "parameters/parameters.hpp"
#include "manager/sim_mediator.hpp"

class FileManager {
public:
    FileManager(AbstractSimulationMediator *mediator);
    ~FileManager();

    void load_params();

    void save_params();
    void save_energies(energies_container_t &);

private:
    static const char PARAMS_FILENAME[];
    static const char ENERGIES_FILENAME[];
    static const char TEXT_FILE_EXTENSION[];
    static const char BINARY_FILE_EXTENSION[];

    void load_from_text_file(std::string filename);
    void load_from_binary_file(std::string filename);
    void save_to_text_file(wavefunction_t &, std::string filename);
    void save_to_binary_file(wavefunction_t &, std::string filename);

    void check_params();

    wavefunction_t psi_loading_buffer;
    AbstractSimulationMediator *mediator;
    PhysicalParameters *params;
};

#endif
