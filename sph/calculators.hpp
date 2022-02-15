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

        // Artificial viscosity params
        const double alpha = 1;
        const double beta = 2;
        const double eta_coeff = 0.01; // multiplied by h^2 in viscosity

        // Equation 2.27 of Bate thesis
        void operator()(Particle &p_i, const ParticleVector &p_vec);
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

        /* Calculate pressure using isothermal equation of state (Bate thesis 2.22)
         * Parameters:
         *      p: calculate the pressure using the density estimate at this particle
         *      c_s: sound speed
         */
        double pressure_isothermal(const Particle &p, double c_s);
        
        // Get sound speed -- just constant value, but disentangled from method to be easily changed
        double sound_speed();
        
        // 
        
        /*
         * Get artificial viscosity Π_ij between two particles. 
         *
         * Parameters:
         *      p_i: first particle,
         *      p_j: second particle
         *      r_ij: distance between first and second particle
         *      h: smoothing length
         *      c_s: sound speed
         * 
         * r_ij and h, and c_s are parameters, rather than calculated insitu, as when this method is
         * called in operator(), they have already been calculated for the pressure, so we reuse
         * them.
         */
        double artificial_viscosity(
            const Particle &p_i, 
            const Particle &p_j, 
            double r_ij,
            double h,
            double c_s
        );
};

#endif