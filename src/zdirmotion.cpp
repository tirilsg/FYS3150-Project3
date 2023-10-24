#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "particle_class.h"
#include "penningtrap_class.h"
#include "logging.h"
#include "zdirmotion.h"

void singleqmotion(const Particle& particle, double max_t, int iterations, const std::string& filename, bool useRK4) {//function that simulates a single particle within a penning trap 
    PenningTrap trap; //no interaction, only one particle 
    trap.particle_add(particle); //adds particle to trap
    double dt = max_t / static_cast<double>(iterations); //defines timestep from time interval and timesteps
    arma::vec t = arma::linspace(0, max_t, iterations);
    arma::mat rs(iterations, 3); //defines position array
    arma::mat vs(iterations, 3); //defines velocity array
    rs.row(0) = particle.r.t(); //the first elements of our rs matrix is the initial position of the particle
    vs.row(0) = particle.v.t(); //the first elements of our vs matrix is the initial velocity of the particle
    for (int i = 1; i < iterations; i++) { //if we want to use runge-kutta for our simulation
        if (useRK4) {
            trap.evolve_RK4(dt, max_t);
        }
        else {//if not, we make use of forward-euler
            trap.evolve_forward_Euler(dt, max_t);
        }
        rs.row(i) = trap.particle_vector[0].r.t(); //for each iteration, we estimate a new position and velocity for the particle and fills our arrays
        vs.row(i) = trap.particle_vector[0].v.t();
    }
    arma::mat data(iterations, 7);
    data.col(0) = t; //defines the first column of our text file to contain the time
    for (int i = 0; i < iterations; i++) { //the data is filled into a new array, that compiles both position and velocity data in the same array
        for (int j = 1; j < 4; j++) {
            data(i, j) = rs(i, j - 1);
            data(i, j + 3) = vs(i, j - 1);
        }
    }
    saveDataToTxt(filename, data); //we call a function that simply saves this data in a txt file
}

