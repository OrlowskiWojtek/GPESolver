#ifndef INITIALIZER_HPP
#define INITIALIZER_HPP

#include "context/context.hpp"
#include "manager/sim_mediator.hpp"
#include "parameters/init_parameters.hpp"
#include <functional>

class DataInitializer {
public:
    DataInitializer(AbstractSimulationMediator *);
    void initialize_wavefunction(const WavefunctionInitConfig& conf);
    void initialize_binary_wavefunctions(const WavefunctionInitConfig& conf_a,
                                         const WavefunctionInitConfig& conf_b);
    void initialize_potential(const PotentialInitConfig& conf);

private:
    //! Init for psi text file
    void init_from_text_file(const WavefunctionInitConfig& conf,
                             wavefunction_t& target);
    //! Init for psi binary file
    void init_from_binary_file(const WavefunctionInitConfig& conf,
                               wavefunction_t& target);
    //! apply already set pote function to _pote container
    void init_pote(const PotentialInitConfig& conf);
    //! apply already set data function to _data container
    void init_wavefunction(const WavefunctionInitConfig& conf,
                           wavefunction_t& target);
    //! applies small random noise to the wavefunction, just to break symmetries
    void add_initial_noise(wavefunction_t& target);
    //! applies small random noise to the potential, just to break symmetries
    //! potential is then pote * noise (noise is the same between experiments)
    void add_initial_pote_noise();

    wavefunction_t _data;
    BinaryMixtureState _binary_data;
    potential_t _pote;
    std::function<double(double, double, double)> _pote_func;
    std::function<std::complex<double>(double, double, double,
                                       const WavefunctionInitConfig&)> _data_func;

    PhysicalParameters* params;
    AbstractSimulationMediator *p_mediator;
    SimulationContext *p_sctx;
};

#endif
