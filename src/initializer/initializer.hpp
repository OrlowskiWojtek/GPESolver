#ifndef INITIALIZER_HPP
#define INITIALIZER_HPP

#include "context/context.hpp"
#include "manager/sim_mediator.hpp"

class DataInitializer {
public:
    DataInitializer(AbstractSimulationMediator*);
    void initialize_wavefunction();

private:
    void init_from_text_file();
    void init_from_binary_file();
    void init_with_gaussian();
    void init_with_multiple_gaussian();
    void init_with_cos();

    wavefunction_t _data;
    PhysicalParameters* params;
    AbstractSimulationMediator* p_mediator;
    SimulationContext* p_sctx;
};

#endif
