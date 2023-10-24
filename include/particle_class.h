#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>

#ifndef __particle_class_h__  
#define __particle_class_h__

using namespace std;
using namespace arma;

class Particle{
  public:
  double q; // particle charge
  double m; // particle mass
  vec r; // position vector
  vec v; // velocity vector
  void set(double q, double m , arma::vec r, arma::vec v){ //a function that lets us assign values to our arguments
      //we call can set the values for a particle class by calling particle.set()
      Particle::q = q;
      Particle::m = m;
      Particle::r = r;
      Particle::v = v;
  }
  bool escape_test(double d){//a test that checks whether a particle has escaped the trap
    bool escaped = false; //the norm is set to false
    double norm = arma::norm(r);
    if (norm > d){//if the distance from origo in xy-plane is larger than d, the particle escapes 
      escaped = true; //and the function returns a value "true" for the test
    }
  return escaped; //we return the test result
  }
};

#endif