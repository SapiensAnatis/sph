/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * particle.hpp defines the basic particle unit. These are set up in an array according to n_part.
 */

#ifndef particle_hpp // Include guard
#define particle_hpp

enum ParticleType {
    Alive,
    Ghost
};

// Used to display particle type as a string, rather than number, in dump files
const char* const ParticleTypeNames[2] = {
    "Alive",
    "Ghost"
};

// 'global' particle counter, so creator of Particle doesn't have to keep track
static int _particle_counter = 0;

struct Particle {
    const int id; // Unique numerical identifier
    double mass; // Should be const, but can't be if not given in constructor :(

    double pos;
    double vel;
    double acc;

    double h; // Smoothing length

    double u; // energy
    double density;
    double pressure;

    ParticleType type;

    // Full initializer for unit tests
    Particle(double pos, double vel, double mass)
        : id(_particle_counter), mass(mass), pos(pos), vel(vel), acc(0), u(0), density(0), pressure(0), type(Alive)
    {
        _particle_counter++;
    }

    // Default initializer for creating arrays
    Particle() : id(_particle_counter), type(Alive)
    {
        _particle_counter++;
    }
    
    // Assignment operator
    Particle& operator =(Particle &p) {
        // Don't copy id. This ensures that id is unique between particles. It does however result
        // in the slightly odd issue with them not starting at 0 due to array reallocation, as
        // described in setup.cpp
        mass = p.mass;
        pos = p.pos;
        vel = p.vel;
        acc = p.acc;
        u = p.u;
        density = p.density;
        pressure = p.pressure;
        type = p.type;

        return *this;
    }

    // Equality operator
    bool operator ==(Particle p) {
        // Check id
        return (id == p.id);
    }

    // Inequality operator
    bool operator !=(Particle p) {
        return (id != p.id);
    }
};

#endif