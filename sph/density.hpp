/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * density.hpp defines DensityCalculator, a class which can operate on a Particle to set its
 * density property. DensityCalculator also contains the function `kernel` which is the density
 * kernel that should be used. The default is the third-order Schoenberg B-spline approximation to
 * the Gaussian, as detailed in Price (2010) pg. 3.
 */


#ifndef density_hpp
#define density_hpp

#include "particle.hpp"
#include "setup.hpp"

// Visitor-like design pattern: DensityCalculator operates on a Particle and sets the `density`
// property by using the smoothing kernel and iterating over neighbours.

class DensityCalculator {
    public:
        Config config;

        void operator()(Particle &p, const ParticleVector &p_vec);
        DensityCalculator(Config c) : config(c) {};
    
    private:
        // w(q)
        double kernel(double q);
        double smoothing_length();
};

#endif