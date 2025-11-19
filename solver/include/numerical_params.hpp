#ifndef NUMERICAL_PARAMS_HPP
#define NUMERICAL_PARAMS_HPP

#include <cstddef>

struct NumericalParameters{
    static constexpr size_t iter_imag_evo = 15000;
    static constexpr size_t iter_real_evo = 100000;

    static constexpr double relax_omega = 0.8;
    static constexpr double imag_time_dt = 1.25e11; //!< \todo might need to fix into some units 
    static constexpr double real_time_dt = 1.00e10;
};

#endif
