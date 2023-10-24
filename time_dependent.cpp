//this program is simulates 100 particles within a penning trap during an interval in time, with and without particle interactions
//due to a huge number of calculations needed to store data, this program will take a long time to run
//to test how it works, adjust the number of particles inserted into the trap, and number of iterations for calculations


#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <cmath>
#include "particle_class.h"
#include "penningtrap_class.h"
#include "zdirmotion.h"
#include "twoqmotion.h"
#include "filltrap.h"

using namespace std;
using namespace arma;

void simulateAndLogData(PenningTrap& trap, const string& filename, const arma::vec& w_v, const arma::vec& fs, double max_t, int iterations) {//function that simulates and logs the data frequency w_v and number of trapped particles as function of ampliture f_amp
    double dt = max_t / (double)iterations; //calculates size of timestep
    arma::vec t = arma::linspace(0, max_t, iterations); //time vector
    arma::mat still_trapped = arma::zeros(w_v.size(), fs.size()); //empty array to be filled with data on amount of particles within system
    for (int k = 0; k < fs.size(); k++) {//loops through all the aplitudes providfed as argument
        for (int j = 0; j < w_v.size(); j++) { //looping through every w_v from the vector function took as argument
            PenningTrap penningtrap = trap; //creates a copy of our penning trap provided as argument
            penningtrap.f_amp = fs(k); //sets the amplitude argument for the trap-copy, constant for the rest of this loop
            penningtrap.w_v = w_v(j); //sets the w_v for the copy of penning trap provided as argument, is set a new for every iteration of loop
            for (int i = 0; i < iterations; i++) {//evolves the system in time, for the amount of iterations decided by argument
                penningtrap.evolve_RK4(dt, t(i));
            }
            still_trapped(j, k) = penningtrap.count_particles(); //counts the amount of particles that are still trapped, and stores the data within our data-array
        }
    }
    arma::mat data = arma::mat(w_v.size(), 4); //we have to create another array, so we can store the w_v data as well
    data.col(0) = w_v; //the first column is to contain w_v data
    for (int j = 0; j < w_v.size(); j++) {
        for (int k = 1; k < 4; k++) {
            data(j, k) = still_trapped(j, k - 1); //stores the data at the right place in our array
        }
    }
    data.save(filename, arma::raw_ascii); //saves the data in file, to be imported in python
}

int main() {
    double max_t = 500.0; 
    int n_particles = 100;
    int iterations = 10000; //a large amount of iterations
    arma::vec fs = arma::vec{0.1, 0.4, 0.7}; //amplitude array
    arma::vec w_v = arma::linspace(0.2, 2.5, 115); //frequency array

    PenningTrap nointer; //defines a penning trap
    nointer.time_dependence = true; //we want time-dependence
    fillPenningTrapWithRandomParticles(nointer, n_particles); //calls the function to fill our trap with n_particles with randomly generated initial conditions
    simulateAndLogData(nointer, "nointertest.txt", w_v, fs, max_t, iterations);

    PenningTrap iter; //defines a penning trap
    iter.Particle_interactions = true; //we want to take into consideration particle interactions
    iter.time_dependence = true;
    fillPenningTrapWithRandomParticles(iter, n_particles); 
    simulateAndLogData(iter, "intertest.txt", w_v, fs, max_t, iterations);

    w_v = arma::linspace(1.0, 1.8, 115); //we zoom into an interesting w_v interval, and then repeats the code to store data as previously
    PenningTrap nointerzoom; 
    nointerzoom.time_dependence = true;
    fillPenningTrapWithRandomParticles(nointerzoom, n_particles);
    simulateAndLogData(nointerzoom, "nointertestzoom.txt", w_v, fs, max_t, iterations);

    PenningTrap iterzoom; 
    iterzoom.Particle_interactions = true;
    iterzoom.time_dependence = true;
    fillPenningTrapWithRandomParticles(iterzoom, n_particles);
    simulateAndLogData(iterzoom, "intertestzoom.txt", w_v, fs, max_t, iterations);
    return 0;
}