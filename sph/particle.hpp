/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * particle.hpp defines the basic particle unit. These are set up in an array according to n_part.
 */

#ifndef particle_hpp // Include guard
#define particle_hpp


enum ParticleType {
    Alive,
    Dead,
    Ghost
};

// 'global' particle counter, so creator of Particle doesn't have to keep track
static int _particle_counter = 0;

struct Particle {
    const int id; // Unique numerical identifier
    double mass; // Should be const, but can't be if not given in constructor :(

    double pos;
    double vel;
    double acc;

    double u; // energy
    double density;
    double pressure;

    ParticleType type;

    // Full initializer for unit tests
    Particle(double pos, double vel, double mass)
        : id(_particle_counter), mass(mass), pos(pos), vel(vel), acc(0), u(0), density(0), type(Alive)
    {
        _particle_counter++;
    }

    // Default initializer for creating arrays
    Particle() : id(_particle_counter), type(Alive)
    {
        _particle_counter++;
    }

    // Inequality operator
    bool operator !=(Particle p) {
        // Check id
        return (id != p.id);
    }
};

#endif