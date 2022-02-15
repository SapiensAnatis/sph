/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * kernel.hpp contains the kernel weighting function and spatial derivatives thereof. They are
 * defined here, in their own class, because several different steps use them, so it is best to have
 * them all call into a single compilation unit where they can be modified and viewed together.
 */

#ifndef kernel_hpp
#define kernel_hpp

#include "setup.hpp"

class Kernel {
    public:
        // W
        static double kernel(double q);
        // grad(W)
        static double d_kernel(double q);
        // h
        static double smoothing_length(Config c);
};

#endif