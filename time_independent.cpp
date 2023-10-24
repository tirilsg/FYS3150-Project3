// program simulating multiple instances of particles within a penning trap
// all simulationsuse a time-independent instance of electric field potential 

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

int main(){
  //Problem 8 : simulation of particle 1 in penning trap in z-direction
  double q = 1; //charge of the Calium particle, unit e
  double m = 40.078; //mass of the Calium particle, unit u
  arma::vec r = {20.0 , 0.0 , 20.0}; //initial position, units \mu m
  arma::vec v = {0.0 , 25.0 , 0.0}; //initial velocity, units \mu m / \mu s
  Particle particle1; //we declare the particle 1
  particle1.set(q,m,r,v); //we define the particle 1
  double time = 50.0; //time interval for our simulation
  double iterations = 4000; //amount of steps we want to run simulation for
  string filename = "particle1motion.txt"; 
  bool useRK4=true; //define that we wish to use Runge-Kutta to estimate position in our simulation 
  singleqmotion(particle1, time, iterations, filename, useRK4); //we call a function that runs our simulation for the conditions specified, and stores data in the .txt doc

  //Problem 8 : simulation of two particles with particle interactions in the (x,y)-plane, as well as logging of phase space
  bool particle_int=false; //we do not want the particles to interact with oneanother
  r = {25.0 , 25.0 , 0.0}; //initial position, units \mu m
  v = {0.0 , 40.0 , 5.0}; //initial velocity, units \mu m / \mu s
  Particle particle2; //we declare the particle 2
  particle2.set(q,m,r,v); //we define the particle 2
  filename="withoutinteraction.txt";
  twoq(particle1, particle2, time, iterations, filename, particle_int); //we call our function that performs simulation and stores data

  //Problem 8 : simulation of two particles without particle interactions in the (x,y)-plane, as well as logging of phase space
  filename="withinteraction.txt";
  particle_int=true; //we let the particles interact, and influence eachothers movement
  twoq(particle1, particle2, time, iterations, filename, particle_int);//we call our function that performs simulation and stores data

  //Problem 8 : error for particle1 using RK4 and Forward Euler, for steps 4000, 8000, 16000, 32000
  arma::vec iterations_vec = {4000, 8000, 16000, 32000}; //vector containing amount of steps we wish to use
  useRK4=true; //we want to use runge-kutta to estimate movement
  for (const int iter : iterations_vec) {
    filename="qRK4"+to_string(iter)+".txt"; //filename describing exactly what the file contains
    singleqmotion(particle1, time, iter, filename, useRK4); //we run simulation and store data in files with appropriate names
  }
  useRK4=false; //we want to use forward-euler to estimate movement
  for (const int iter : iterations_vec) {
    filename="qeuler"+to_string(iter)+".txt";
    singleqmotion(particle1, time, iter, filename, useRK4);//we run simulation and store data in files with appropriate names
  }
  return 0;
}