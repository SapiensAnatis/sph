/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * kernel.hpp implements the M_4 kernel that is used as the weighting function. It is in its own
 * separate file so that it is easy to find and modify (as the choice of kernel is based to some
 * extent on personal taste). It also makes the dependency structure a bit less tangled.
 */

#ifndef kernel_hpp // Include guard
#define kernel_hpp

const double KERNEL_RADIUS = 2.5;

double kernel(double q);
double dkernel_dq(double q);

#endif