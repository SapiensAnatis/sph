/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * calculators.hpp defines the functions which calculate specific quantities such as the weighting
 * function, the density evaluation function, and the acceleration evaluation function.
 */


#ifndef calculators_hpp
#define calculators_hpp

#include "basictypes.hpp"

// Calculators adopt a visitor design pattern. This is so that they can be instantiated and store
// certain information that would otherwise be needed for every function call e.g. particle array
// pointer, Config data, etc. 

// Base type of calculator. Defines constructor (storing config and particle array) and an
// override-able operator method
class Calculator {
    public:
        // ctor
        Calculator(const Config c, const ParticleArrayPtr p_arr_ptr) 
            : config(c), p_arr(p_arr_ptr) {}
        // Calculation function
        virtual void operator()(Particle &p) {
            throw new std::logic_error("Attempt to call un-implemented operator() function!");
        }
        // Update function: supply new pointer and config (due to ghost particle reinitialization
        // changing config.n_part and possibly reallocating the array)
        void update(Config c, ParticleArrayPtr p_arr_ptr) {
            config = c;
            p_arr.swap(p_arr_ptr);
        }
    protected:
        Config config;
        ParticleArrayPtr p_arr;

        // Calculate the gradient of W between p_i, p_j with respect to the coordinates of p_i.
        // Used in acceleration and energy calculators.
        double grad_W(const Particle &p_i, const Particle &p_j, double h);

};

class DensityCalculator : public Calculator {
    public:
        // ctor -- just call base class
        DensityCalculator(const Config &c, ParticleArrayPtr p_arr_ptr) 
            : Calculator(c, p_arr_ptr) {};
        
        // Calculate the smoothing length for a particle and then the density. This void method sets
        // the properties on p.
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

        // Get sound speed -- seperate method to be disentangled from pressure calculation methods.
        // Public as it's used to set initial velocities in setup.cpp for the adiabatic test.
        // If isothermal pressure calculation is enabled, then this just returns 1. If
        // adiabatic pressure calculation is enabled, it uses sqrt(gamma * pressure / density)
        double sound_speed(Particle p);

    protected:
        /* Calculate pressure using isothermal equation of state (Bate thesis 2.22)
         * Parameters:
         *      p: calculate the pressure using the density estimate at this particle
         *      c_s: sound speed
         */
        double pressure_isothermal(const Particle &p, double c_s);

        /* Calculate pressure using adiabatic equation of state (Bate thesis 2.23)
         * Parameters:
                p: calculate the pressure using the density and internal energy of this particle
         */
        double pressure_adiabatic(const Particle &p);
        
        
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

        // Check that p's density is not less than epsilon, and throw an error if it is.
        void ensure_nonzero_density(const Particle &p);
};

// Inherit from AccelerationCalculator instead of base Calculator, as we require use of artificial
// viscosity and sound speed methods
class EnergyCalculator : public AccelerationCalculator {
    public:
        // ctor -- just call base class
        EnergyCalculator(const Config &c, ParticleArrayPtr p_arr_ptr) 
            : AccelerationCalculator(c, p_arr_ptr) {};
            
        // Calculate du/dt for a particle and set it as a property
        void operator()(Particle &p) override;
};

#endif