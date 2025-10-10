#ifndef BASIC_TYPES_HPP
#define BASIC_TYPES_HPP

#include <array>

class Point3D {
    std::array<double, 3> data;

public:
    double &x() {
        return data[0];
    };

    double &y() {
        return data[1];
    };

    double &z() {
        return data[2];
    };
};

#endif
