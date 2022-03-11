/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * density.cpp implements the class methods defined in density.hpp.
 */

#include "density.hpp"
#include <cmath>

void DensityCalculator::operator()(Particle &p, const ParticleVector &p_vec) {
    double h = smoothing_length();
    double density = 0;

    for (Particle p_i : p_vec) {
        double q = std::abs(p_i.pos - p.pos) / h;
        double w = this->kernel(q);
        
        density += p_i.mass * (w / h);
        
    }

    p.density = density;
}

// This returns a constant for now, but may later use variable smoothing lengths
double DensityCalculator::smoothing_length() {
    return 0.3* this->config.d_unit;
}

double DensityCalculator::kernel(double q) {
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