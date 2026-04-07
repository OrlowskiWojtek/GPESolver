#ifndef ABSTRACT_SIMULATION_MEDIATOR_HPP
#define ABSTRACT_SIMULATION_MEDIATOR_HPP

#include "context/context.hpp"

class AbstractSimulationMediator{
public:
    virtual void on_data_loaded(const wavefunction_t&) = 0;
    virtual void on_data_initialized(const wavefunction_t& ) = 0;
    virtual void on_pote_initialized(const potential_t& ) = 0;
    virtual void on_pote_changed(const potential_t& ) = 0;
    virtual void on_params_loaded() = 0;

    virtual void request_load_from_binary(wavefunction_t&) = 0;
    virtual void request_load_from_text(wavefunction_t&) = 0;
    virtual void request_free_potential() = 0;

    virtual void save_data(const wavefunction_t&) = 0;
    virtual void save_checkpoint(const wavefunction_t&) = 0;
    virtual void save_initial_state(const wavefunction_t&) = 0;
    virtual void save_energies(const energies_container_t&) = 0;
};

#endif
