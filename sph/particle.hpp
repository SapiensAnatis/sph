/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * particle.hpp defines the basic particle unit. These are set up in an array according to n_part.
 */

#ifndef particle_hpp // Include guard
#define particle_hpp

#include <vector>

enum ParticleType {
    Alive,
    Dead,
    Ghost
};

// 'global' particle counter, so creator of Particle doesn't have to keep track
static int _particle_counter = 0;

struct Particle {
    const int id; // Unique numerical identifier

    double pos;
    double vel;
    double acc;

    double mass;
    double u; // energy
    double density;

    ParticleType type;

    Particle(double pos, double vel)
        : id(_particle_counter), pos(pos), vel(vel), acc(0), mass(1), u(0)
    {
        _particle_counter++;
    }

    // Inequality operator
    bool operator !=(Particle p) {
        // Check id
        return (this->id != p.id);
    }
};

typedef std::vector<Particle> ParticleVector;

#endif