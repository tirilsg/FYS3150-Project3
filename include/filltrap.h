#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "particle_class.h"
#include "penningtrap_class.h"

#ifndef __filltrap_h__  
#define __filltrap_h__

void fillPenningTrapWithRandomParticles(PenningTrap& trap, int number_of_particles);

#endif