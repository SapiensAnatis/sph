/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * basictypes.hpp defines two structures that are used by almost every other file: Particle, Config,
 * and also ParticleArrayPtr
 */

#ifndef basictypes_hpp // Include guard
#define basictypes_hpp

#include <memory>
#include <iostream>
#include <boost/shared_ptr.hpp>

// ===== CONFIG =====

enum PressureCalc {
    Isothermal,
    Adiabatic
};

struct Config {
    int n_part;
    double mass;
    PressureCalc pressure_calc;
    double limit;
    double v_0;
    double h_factor;
    double t_i;
    // Runtime properties; not set from ConfigReader
    int n_ghost; // Number of ghost particles
};

// ===== PARTICLES ===== 

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

    double du_dt; // Thermal energy derivative w.r.t time
    double u; // Thermal energy
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
        // Don't copy id. This ensures that id is unique between particles.
        mass = p.mass;
        pos = p.pos;
        vel = p.vel;
        acc = p.acc;
        h = p.h;
        u = p.u;
        du_dt = p.du_dt;
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

    // ostream operator -- enables 'printing' of particles
    friend std::ostream& operator <<(std::ostream& os, const Particle& p) {
        return os << "<Particle> id: " << p.id 
                  << " type: " << ParticleTypeNames[p.type] 
                  << " pos: " << p.pos 
                  << " vel: " << p.vel
                  << " density: " << p.density
                  << " h: " << p.h;
    }
};

// Particle array
typedef boost::shared_ptr<Particle []> ParticleArrayPtr;

#endif