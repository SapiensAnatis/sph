/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * kernel.cpp implements the class methods defined in kernel.hpp.
 */

#include <cmath>
#include "kernel.hpp"

double Kernel::kernel(double q) {
    double w;
    double sigma = 2.0/3.0;

    // Third-order Schoenberg B-spline (M_4), as detailed in Price (2010) pg. 3.
    if (q < 1) {
        w = 0.25 * std::pow(2 - q, 3) - std::pow(1 - q, 3);
        return sigma * w;
    } else if (q >= 1 && q < 2) {
        w = 0.25 * std::pow(2 - q, 3);
        return sigma * w;
    } else {
        // > 2
        return 0;
    }
}

double Kernel::d_kernel(double q) {
    double w;
    double sigma = 2.0/3.0;

    // Derivative of kernel(), analytically
    if (q < 1) {
        w = 0.25 * -3 * std::pow(2 - q, 2) + 3 * std::pow(1 - q, 2);
        return sigma * w;
    } else if (q >= 1 && q < 2) {
        w = 0.25 * -3 * std::pow(2 - q, 2);
        return sigma * w;
    } else {
        // > 2
        return 0;
    }
}

// This returns a constant for now, but may later use variable smoothing lengths
double Kernel::smoothing_length(const Config &c) {
    return c.d_unit;
}
