/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * kernel.cpp implements the functions from kernel.hpp.
 */

#ifndef kernel_hpp // Include guard
#define kernel_hpp

#include <cmath>
#include <stdexcept>
#include <iostream>

double kernel(double q) {
    double w;
    double sigma = 1./24.;

    if (q < 0) {
        // The kernel is symmetric, so we expect absolute values only!
        std::cerr << "[ERROR] Illegal attempt to evaluate kernel at q = " << q << std::endl;
        throw std::invalid_argument("Negative argument given to kernel!");
    }

    // M_5 quartic Schoenberg B-spline, as detailed in Price (2012) pg. 3.
    if (q < 0.5) {
        w = std::pow(5./2. - q, 4) - 5 * std::pow(3./2. - q, 4) + 10 * std::pow(1./2. - q, 4);
        return sigma * w;
    } else if (q >= 0.5 && q < 3./2.) {
        w = std::pow(5./2. - q, 4) - 5 * std::pow(3./2. - q, 4);
        return sigma * w;
    } else if (q >= 3./2. && q < 5./2.) {
        w = std::pow(5./2. - q, 4);
        return sigma * w;
    } else {
        // > 2.5
        return 0;
    }
}

double dkernel_dq(double q) {
    double w;
    double sigma = 2.0/3.0;

    if (q < 0) {
        // The kernel is symmetric, so we expect absolute values only!
        std::cerr << "[ERROR] Illegal attempt to evaluate kernel derivative at q = " << q << std::endl;
        throw std::invalid_argument("Negative argument given to kernel!");
    }

    // Derivative of kernel(), analytically
    if (q < 0.5) {
        w = -4 * std::pow(5./2. - q, 3) + 20 * std::pow(3./2. - q, 3) - 40 * std::pow(1./2. - q, 3);
        return sigma * w;
    } else if (q >= 0.5 && q < 3./2.) {
        w = -4 * std::pow(5./2. - q, 3) + 20 * std::pow(3./2. - q, 3);
        return sigma * w;
    } else if (q >= 3./2. && q < 5./2.) {
        w = -4 * std::pow(5./2. - q, 3);
        return sigma * w;
    } else {
        // > 2.5
        return 0;
    }
}

#endif