#ifndef ABSTRACT_SIMULATION_MEDIATOR_HPP
#define ABSTRACT_SIMULATION_MEDIATOR_HPP

#include "context/context.hpp"
#include <string>

class AbstractSimulationMediator{
public:
    virtual void on_data_loaded(const wavefunction_t&) = 0;
    virtual void on_data_initialized(const wavefunction_t& ) = 0;
    virtual void on_binary_data_initialized(const BinaryMixtureState&) = 0;
    virtual void on_pote_initialized(const potential_t& ) = 0;
    virtual void on_pote_changed(const potential_t& ) = 0;
    virtual void on_params_loaded() = 0;

    virtual void request_load_from_binary(const std::string&, wavefunction_t&) = 0;
    virtual void request_load_from_text(const std::string&, wavefunction_t&) = 0;
    virtual void request_free_potential() = 0;

    virtual void save_data(const wavefunction_t&) = 0;
    virtual void save_checkpoint(const wavefunction_t&) = 0;
    virtual void save_initial_state(const wavefunction_t&) = 0;
    virtual void save_binary_checkpoint(const BinaryMixtureState&) = 0;
    virtual void save_binary_initial_state(const BinaryMixtureState&) = 0;
    virtual void save_energies(const energies_container_t&) = 0;
};

#endif
