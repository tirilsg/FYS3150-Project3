# FYS3150-Project3

This repository contains c++ code constructing a model for a penning trap, used to run multiple simulations testing both the accuracy of the model, and the physics at work within the system. The majority of the code is devided into files of their own, with header and source files sorted in maps `src` (source files) and `include` (header files), which is has to be linked together when testing the code. The code is split into multiple files: 

- ### `particle_class.h`:
Contains a class `Particle` that is used to define a particle by mass, charge, a position vector and a velocity vector. This class also contains a function `escape_test()` that checks whether the particle exsists within the bounds of the penning trap. 

- ### `penningtrap_class.h`:
Contains a definition of the class `PenningTrap`, that defines the environment within a penning trap, and uses this to define methods `evolve_forward_Euler(dt,t)` and `evolve_RK4(dt,t)` that estimates particle's movement within the trap, with the method `evolve_RK4(dt,t)` set as the default trajectory-estimation-method. The class contains definitions of booleans that defines whether the system is modelled as time-dependent (which is set to `false` by default), and whether the system estimates particle-trajectories by taking into consideration forces between particles that exsist within the trap (also `false` by default). Also contains a method `particle_add()` that adds a particle into our system, and a function `count_particles()` that conts the amount of particles within the system.

- ### `filltrap.cpp`:
Contains a single function `fillPenningTrapWithRandomParticles(trap, number_of_particles)` that adds a number of Calcium ions into a system PenningTrap trap. 

- ### `logging.cpp`:
Contains a single function `saveDataToTxt(filename, data)` that takes an arbitrary filename, as well as a set of data, and writes the data into the file. 

- ###  `twoqmotion.cpp`:
Contains code for a function `twoq(article1, particle2, max_t, iterations,filename, particle_int)` that adds two particles into a penning trap, and simulates their movement within the trap in an intervall of time, and saves the data into a file by calling  `saveDataToTxt(filename, data)`. 

- ###  `zdirmotion.cpp`:
Contains the implementations of a function `singleqmotion(particle, max_t, iterations, filename, useRK4)` that simulates the movement of a single particle within the environment of a penning trap. This function takes a boolean useRK4 as argument, that dictates whether we use the method `evolve_RK4(dt,t)` to estimate the particle trajectory or not. The data is saved in a file.



This repository contains the model for the penning trap, constructed in c++, used to


`#include <iostream>`
