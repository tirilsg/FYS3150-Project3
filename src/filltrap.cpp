#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "particle_class.h"
#include "penningtrap_class.h"
#include "filltrap.h"

void fillPenningTrapWithRandomParticles(PenningTrap& trap, int number_of_particles) { //function that takes a trap and number of particles, and fills the trap 
    for (int i = 0; i < number_of_particles; ++i) { //loops through number_of_particles instances of particles to be added to trap
        Particle particle; //defines paticle
        double q = 1;  //charge
        double m = 40.078; //mass 
        arma::vec r = vec(3).randn() * 0.1 * trap.d; //randomly generated initial position
        arma::vec v = vec(3).randn() * 0.1 * trap.d; //randomly generated initial velocity
        particle.set(q,m,r,v); //defines particle
        trap.particle_add(particle); //adds particle to trap
    }
}
