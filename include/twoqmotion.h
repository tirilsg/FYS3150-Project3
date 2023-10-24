#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "particle_class.h"
#include "penningtrap_class.h"
#include "logging.h"

#ifndef __twoqmotion_h__  
#define __twoqmotion_h__

using namespace std;
using namespace arma;

void twoq(const Particle& particle1, const Particle& particle2,double max_t, int iterations, const std::string& filename, bool particle_int);

#endif