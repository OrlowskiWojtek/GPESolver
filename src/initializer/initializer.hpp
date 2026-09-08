#ifndef INITIALIZER_HPP
#define INITIALIZER_HPP

#include "context/context.hpp"
#include "manager/sim_mediator.hpp"
#include <functional>

class DataInitializer {
public:
    DataInitializer(AbstractSimulationMediator *);
    void initialize_wavefunction();
    void initialize_potential();
    void change_potential(std::string pote_key);

private:
    //! Init for psi text file
    void init_from_text_file();
    //! Init for psi binary file
    void init_from_binary_file();
    //! apply already set pote function to _pote container
    void init_pote();
    //! apply already set data function to _data container
    void init_wavefunction();

    wavefunction_t _data;
    potential_t _pote;
    std::function<double(double, double, double)> _pote_func;
    std::function<std::complex<double>(double, double, double)> _data_func;

    PhysicalParameters *params;
    AbstractSimulationMediator *p_mediator;
    SimulationContext *p_sctx;
};

#endif
