# Smoothed particle hydrodynamics

## Overview

This program is a 1-dimensional implementation of the smoothed particle hydrodynamics (SPH) method of simulating colliding shock streams. It was originally written for a coursework assignment, and the associated report which was submitted alongside it is included as a PDF file.

SPH is a mesh-free, Lagrangian-based approach to computationally simulating fluid dynamics. It models the motion of point-like particles individually. The term 'smoothed' in the name refers to the fact that physical quantities such as acceleration are calculated by iterating over a particle's neighbours, using a smoothing length (which is adjusted up or down for low/high densities, respectively, to maintain a roughly constant number of neighbours). 

This work is largely based upon the methods described in reviews by Price (2012) and Rosswog (2009), and aims to reproduce analytical results described in chapter 3 of a thesis by Bate (1995). 

## Building

There are two ways of building the program. There is the fancy build system, using Bazel, which also enables unit testing. But there is also a backup Makefile in the sph/ directory that allows the code to be run with `make && ./sph`. Please see the below information about `define.hpp` for an associated warning if you are using `make` and plan to tinker with the code.

When invoked, the program takes one positional argument, which is a path to a config file. If it doesn't find it, it'll just use "./config", which works fine when using `make`, but since Bazel puts the binary in some weird directory, you may need to pass a hardcoded path e.g. `bazel run -- /full/path/to/config.txt`

The program should run fine and doesn't require any particularly esoteric external dependencies or libraries -- the main ones are GNU Scientific Library and a C++17 compiler. Google Test is used for the unit tests, but the Bazel build system automatically downloads that (I think).

## Reflections

Having reached the end of this project, there are two things I wish I would have done differently.

Firstly, I wish I had been more strict with unit testing. The end phases of the project had me debugging a lot of logical errors that were very tricky to track down. I gave up on unit testing much earlier because I didn't really understand how to isolate all the various pieces of the code to unit test while keeping things clustered together and sub-methods private with object-oriented programming. I also wasn't clear about how to get analytical results for some of the more mundane quantities, like the gradient of the weighting function. It would have been far easier to isolate the logical errors that plagued the later stages of this project with a robust unit testing framework, and perhaps even a routine in the integrator that checked momentum/energy conservation.

Secondly, I regret using shared pointers for the array, and I would have benefitted more generally from planning out the project and the requirements for data structures before writing the first lines of code. I thought a vanilla array would work well, and be more efficient, because if I was just creating all the particles once and never deleting them, I surely wouldn't need any of the fancy rearranging/insertion/etc. methods of a vector, right? Wrong! I later realized that when creating ghost particles, you will need to modify the size of the array, which lead to some unnecessarily complicated code in ghost_particles.cpp. Moreover, using an array meant that I had to pass the Config struct almost *everywhere*, because there was no intrinsic size information stored in the pointer, and I was iterating over the array almost everywhere (such is the nature of SPH). Though, in hindsight, I very easily could have made a custom struct that contained the array and some size information (...at that point, having also done resizing for the ghost particles, I'd have pretty much just recreated a vector but with my own worse code).

## Reviewing/editing

A brief overview of what each file contains is as follows:

- config.txt: Sets runtime properties, such as number of particles, timestep, boundary size, adiabatic/isothermal etc.
- basictypes.hpp: Defines the Config and Particle struct, which are types used in almost every other file
- calculators.cpp/hpp: Defines DensityCalculator, AccelerationCalculator, and EnergyCalculator, which are called into by the integrator as well as the setup. This is where the bulk of the maths happens and is where most equations are implemented.
- define.hpp: Defines some compile-time settings and constants for the program such as whether to use variable smoothing lengths, and whether to print root-finding diagnostic messages. WARNING: If any of these settings are changed, and you are using `make`, it is highly advisable to do a clean build afterwards (`make clean && make`) as make will otherwise re-use .o files compiled under old settings.
- ghost_particles.cpp/hpp: Contains the method to set up the ghost particles, which is done on setup and also in the middle of each timestep.
- kernel.cpp/hpp: Contains the SPH smoothing kernel.
- main.cpp: The main entrypoint for the program.
- plot.py: Sample plotting code to visualize the results of the program.
- setup.cpp/hpp: Contains the code that sets up the initial conditions of the simulation and the particle array. Called into by main.cpp.
- smoothing_length.cpp/hpp: Contains the root-finding algorithm that enables variable smoothing lengths, as well as a method to calculate 'omega' parameters (since both require calculating dW/dh).
- sph_simulation.cpp/hpp: Provides the integrator (velocity Verlet) and also file output routines.

## Bibliography

- Price, Daniel J. "Smoothed particle hydrodynamics and magnetohydrodynamics." Journal of Computational Physics 231.3 (2012): 759-794.
- Rosswog, Stephan. "Astrophysical smooth particle hydrodynamics." New Astronomy Reviews 53.4-6 (2009): 78-104.
- Bate, M. "Publication: Ph. D. Thesis Pub Date: 1995." (http://www.astro.ex.ac.uk/people/mbate/Preprints/thesis/thesis.html)