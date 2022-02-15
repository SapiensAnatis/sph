/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * density.cpp implements the class methods defined in density.hpp.
 */

#include <cmath>

#include "density.hpp"
#include "kernel.hpp"

void DensityCalculator::operator()(Particle &p, const ParticleVector &p_vec) {
    double h = smoothing_length();
    double density = 0;

    for (Particle p_i : p_vec) {
        if (p_i != p) { // i != j
            double q = std::abs(p_i.pos - p.pos) / h;
            double w = Kernel::kernel(q);
            
            density += p_i.mass * (w / h);
        }
    }

    p.density = density;
}