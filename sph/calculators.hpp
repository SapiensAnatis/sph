/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * calculators.hpp defines a few classes which set properties on Particles. For the moment, these
 * are DensityCalculator and AccelerationCalculator. These are visitor-like classes which operate
 * on a particle and return void, instead setting a property on the object.
 */


#ifndef calculators_hpp
#define calculators_hpp

#include "calculators.hpp"
#include "setup.hpp"

// Visitor-like design pattern: DensityCalculator operates on a Particle and sets the `density`
// property by using the smoothing kernel and iterating over neighbours.

class DensityCalculator {
    public:
        const Config &config;

        void operator()(Particle &p, const ParticleVector &p_vec);
        DensityCalculator(const Config &c) : config(c) {};
};

class AccelerationCalculator {
    public:
        const Config &config;

        // Equation 2.27 of Bate thesis
        void operator()(Particle &p, const ParticleVector &p_vec);
        AccelerationCalculator(const Config &c) : config(c) {};

    private:
        /*
         * These 'intermediate quantities' could have their own class as well as be properties of
         * Particle. After all, density was -- so why not everything else?
         * 
         * I decided against this because density is probably the only quantity that would have to
         * be analyzed (as per the project brief). Also, having extraneous fields on the Particle
         * would increase the memory footprint for each particle slightly, and this would accumulate
         * quickly for a simulation with many thousands of particles.
         */

        // Calculate pressure using isothermal equation of state (Bate thesis 2.22)
        double pressure_isothermal(const Particle &p);
        // Get sound speed -- just constant value, but disentangled from method to be easily changed
        double sound_speed();
};

#endif