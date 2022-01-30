struct Particle {
    /* 1-dimensional values */
    double pos;
    double vel;
    double density;
};

// There are no strings or other fields that require manual memory allocation
// for Particle, so actually a factory is unneeded
// Particle particle_factory(double pos, double vel);