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
    //! Init with single normalized gaussian
    void init_with_gaussian();
    //! Init with multiple gaussian system
    void init_with_multiple_gaussian();
    //! Init with spatially separated gaussians in 3D (in format nx | ny | nz)
    void init_with_setup_gaussian();
    //! Init with cosinus
    void init_with_cos();

    //! rotate centers around z axis
    void rotate_centers(std::vector<double> &centers_x, std::vector<double> &centers_y, double degrees);

    //! apply already set pote function to _pote container
    void init_pote();

    wavefunction_t _data;
    potential_t _pote;
    std::function<double(int, int, int)> _pote_func;

    PhysicalParameters *params;
    AbstractSimulationMediator *p_mediator;
    SimulationContext *p_sctx;
};

#endif
