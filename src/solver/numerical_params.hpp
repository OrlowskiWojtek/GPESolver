#ifndef NUMERICAL_PARAMS_HPP
#define NUMERICAL_PARAMS_HPP

#include <cstddef>

struct NumericalParameters {
    static constexpr size_t iter_imag_speed_test = 500;
    static constexpr size_t iter_real_speed_test = 500;

    static double imag_time_dt;
    static double real_time_dt;
};

#endif
