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
        DensityCalculator(Config &c) : config(c) {};
};

class AccelerationCalculator {
    public:
        const Config &config;

        void operator()(Particle &p, const ParticleVector &p_vec);
        AccelerationCalculator(Config &c) : config(c) {};
};

#endif