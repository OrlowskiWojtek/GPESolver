#ifndef FILE_MANAGER_HPP
#define FILE_MANAGER_HPP

#include "context/context.hpp"
#include "parameters/parameters.hpp"
#include "manager/sim_mediator.hpp"
#include <functional>

//! \class FileManager
//!
//! \brief Manages all io to files.
//!
//! \details
//! \todo Auto specify whether load from text or binary file based on extension.
class FileManager {
public:
    //! \brief Constructor.
    //! \param mediator Pointer to AbstractSimulationMediator.
    FileManager(AbstractSimulationMediator *mediator);
    
    //! \brief Destructor.
    ~FileManager();

    //! \brief Loads parameters from file.
    void load_params();

    //! \brief Saves parameters to file.
    void save_params();
    
    //! \brief Saves energies to file.
    //! \param energies Container of energies to save.
    void save_energies(const energies_container_t &energies);

    //! \brief Loads wavefunction from a text file.
    //! \param filename Name of the text file.
    void load_from_text_file(std::string filename);
    
    //! \brief Loads wavefunction from a binary file.
    //! \param filename Name of the binary file.
    void load_from_binary_file(std::string filename);
    
    //! \brief Loads wavefunction from a text file and applies a callback.
    //! \param filename Name of the text file.
    //! \param callback Function to apply to the loaded wavefunction - usually copy to another data buffer.
    void load_from_text_file(std::string filename, std::function<void(wavefunction_t&)> callback);
    
    //! \brief Loads wavefunction from a binary file and applies a callback.
    //! \param filename Name of the binary file.
    //! \param callback Function to apply to the loaded wavefunction - usually copy to another data buffer.
    void load_from_binary_file(std::string filename, std::function<void(wavefunction_t&)> callback);

    //! \brief Saves wavefunction to a text file.
    //! \param psi Wavefunction to save.
    //! \param filename Name of the text file.
    void save_to_text_file(const wavefunction_t &psi, std::string filename);
    
    //! \brief Saves wavefunction to a binary file.
    //! \param psi Wavefunction to save.
    //! \param filename Name of the binary file.
    void save_to_binary_file(const wavefunction_t &psi, std::string filename);

private:
    //! \brief Name of the parameters file.
    static const char PARAMS_FILENAME[];
    //! \brief Name of the energies file.
    static const char ENERGIES_FILENAME[];
    //! \brief Extension for text files.
    static const char TEXT_FILE_EXTENSION[];
    //! \brief Extension for binary files.
    static const char BINARY_FILE_EXTENSION[];

    //! \brief Checks parameters for validity.
    void check_params();

    //! \brief Checks for best number of n for FFTW
    bool is_fft_compatible(int n);

    //! \brief Buffer for loading wavefunctions.
    wavefunction_t psi_loading_buffer;
    //! \brief Pointer to the simulation mediator.
    AbstractSimulationMediator *mediator;
    //! \brief Pointer to the physical parameters.
    PhysicalParameters *params;
};

#endif
