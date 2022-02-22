/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * calculators.hpp defines the functions which calculate specific quantities such as the weighting
 * function, the density evaluation function, and the acceleration evaluation function.
 */


#ifndef calculators_hpp
#define calculators_hpp

#include "calculators.hpp"
#include "setup.hpp"

class Kernel {
    public:
        // W
        static double kernel(double q);
        // grad(W)
        static double d_kernel(double q);
        // h. Not much point in this until I get around to doing dynamic smoothing lengths
        static double smoothing_length(const Config &c);
};


// Calculators adopt a visitor design pattern. This is so that they can be instantiated and store
// certain information that would otherwise be needed for every function call e.g. particle array
// pointer, Config data, etc. 

// Base type of calculator
class Calculator {
    public:
        // Calculation function
        virtual void operator()(Particle &p) {
            throw new std::logic_error("Attempt to call un-implemented operator() function!");
        }

        Calculator(const Config &c, ParticleArrayPtr p_arr_ptr) 
            : config(c), p_all(p_arr_ptr) {}
    protected:
        const Config config;
        const ParticleArrayPtr p_all;

};


class DensityCalculator : public Calculator {
    public:
        // ctor -- just call base class
        DensityCalculator(const Config &c, ParticleArrayPtr p_arr_ptr) 
            : Calculator(c, p_arr_ptr) {};
        // Equation 2.21 of Bate thesis
        void operator()(Particle &p) override;
};

class AccelerationCalculator : public Calculator {
    public:
        // ctor
        AccelerationCalculator(const Config &c, ParticleArrayPtr p_arr_ptr) 
            : Calculator(c, p_arr_ptr) {};
        // Artificial viscosity params
        const double alpha = 1;
        const double beta = 2;
        const double eta_coeff = 0.01; // multiplied by h^2 in viscosity

        // Equation 2.27 of Bate thesis
        void operator()(Particle &p_i) override;

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
        
        // Get sound speed -- just constant value, but disentangled from method to be modifiable
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