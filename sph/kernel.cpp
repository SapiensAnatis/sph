/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * kernel.cpp implements the functions from kernel.hpp.
 */

#ifndef kernel_hpp // Include guard
#define kernel_hpp

#include <cmath>

double kernel(double q) {
    double w;
    double sigma = 2.0/3.0;

    // Third-order Schoenberg B-spline (M_4), as detailed in Price (2012) pg. 3.
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

double dkernel_dq(double q) {
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

#endif