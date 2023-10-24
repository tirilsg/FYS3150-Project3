#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "particle_class.h"
#include "penningtrap_class.h"
#include "logging.h"
#include "twoqmotion.h"

void twoq(const Particle& particle1, const Particle& particle2, double max_t, int iterations, const std::string& filename, bool particle_int) {
    PenningTrap trap; //defines the trap
    if(particle_int == true){
        trap.Particle_interactions=true;
    }
    trap.particle_add(particle1); //adds both particles into the trap environment
    trap.particle_add(particle2);
    double dt = max_t / static_cast<double>(iterations); //defines size of timestep
    arma::vec t = arma::linspace(0, max_t, iterations); //defines time-vector
    arma::mat data(iterations, 13); //defines a huge matrix data, that will store our data for both particles
    for (int i = 0; i < iterations; i++) {//looping through each step in time
        trap.evolve_RK4(dt, max_t); //evolves the system with time using runge-kutta
        const auto& r1 = trap.particle_vector[0].r; //constructs information-vectors
        const auto& v1 = trap.particle_vector[0].v;
        const auto& r2 = trap.particle_vector[1].r;
        const auto& v2 = trap.particle_vector[1].v;
        data(i, 0) = t(i); //fills the data-array with time for this timestep
        data(i, span(1, 3)) = r1.t(); //then fills the data array with estimations of position and velcity of particle 1
        data(i, span(4, 6)) = v1.t();
        data(i, span(7, 9)) = r2.t();//then fills the data array with estimations of position and velcity of particle 2
        data(i, span(10, 12)) = v2.t();
    }
    saveDataToTxt(filename, data); //we call a function that simply saves this data in a txt file
}

