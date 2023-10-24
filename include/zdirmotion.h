#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "particle_class.h"
#include "penningtrap_class.h"
#include "logging.h"

#ifndef __zdirmotion_h__  
#define __zdirmotion_h__

using namespace std;
using namespace arma;

void singleqmotion(const Particle& particle,double max_t, int iterations, const std::string& filename, bool useRK4);

#endif