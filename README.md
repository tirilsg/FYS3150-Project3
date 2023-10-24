# FYS3150-Project3

This repository contains c++ code constructing a model for a penning trap, used to run multiple simulations testing both the accuracy of the model, and the physics at work within the system. The majority of the code is devided into files of their own, with header and source files sorted in maps `src` (source files) and `include` (header files), which is has to be linked together when testing the code. The code is split into multiple files: 

----------------------

### `particle_class.h`:
Contains a class `Particle` that is used to define a particle by mass, charge, a position vector and a velocity vector. This class also contains a function `escape_test()` that checks whether the particle exsists within the bounds of the penning trap. 

### `penningtrap_class.h`:
Contains a definition of the class `PenningTrap`, that defines the environment within a penning trap, and uses this to define methods `evolve_forward_Euler(dt,t)` and `evolve_RK4(dt,t)` that estimates particle's movement within the trap, with the method `evolve_RK4(dt,t)` set as the default trajectory-estimation-method. The class contains definitions of booleans that defines whether the system is modelled as time-dependent (which is set to `false` by default), and whether the system estimates particle-trajectories by taking into consideration forces between particles that exsist within the trap (also `false` by default). Also contains a method `particle_add()` that adds a particle into our system, and a function `count_particles()` that conts the amount of particles within the system.

### `filltrap.cpp`:
Contains a single function `fillPenningTrapWithRandomParticles(trap, number_of_particles)` that adds a number of Calcium ions into a system PenningTrap trap. 

###`logging.cpp`:
Contains a single function `saveDataToTxt(filename, data)` that takes an arbitrary filename, as well as a set of data, and writes the data into the file. 

### `twoqmotion.cpp`:
Contains code for a function `twoq(article1, particle2, max_t, iterations,filename, particle_int)` that adds two particles into a penning trap, and simulates their movement within the trap in an intervall of time, and saves the data into a file by calling  `saveDataToTxt(filename, data)`. 

### `zdirmotion.cpp`:
Contains the implementations of a function `singleqmotion(particle, max_t, iterations, filename, useRK4)` that simulates the movement of a single particle within the environment of a penning trap. This function takes a boolean useRK4 as argument, that dictates whether we use the method `evolve_RK4(dt,t)` to estimate the particle trajectory or not. The data is saved in a file.


------------------------

The simulations themselves is done by calling the functions implemented in the files above, which we simply link to each of our main programs. We decided to split the program code into two seperate files, dependent on time-dependency, since the program `time_dependent.cpp` takes a long time to run.

### `time_independent.cpp`:
All the simulations done in this program is independent of time.
This file contains code that estimates the movement of a single particle 1 in a penning trap, two particles 1 and 2 in a penning trap with and without particle interactions, and simulations for usage of both `evolve_forward_Euler(dt,t)` and `evolve_RK4(dt,t)` for a single particle 1. 


### `time_dependent.cpp`:
All the simulations done in this program is dependent of time.
This program defines a function `simulateAndLogData(trap, filename, w_v, fs, max_t, iterations)` that fills a penning trap with particles, and estimates the trajectories of each of these particles by `evolve_RK4(dt,t)` and estimates how many particles are still trapped, for a vector containing different values for amplitudes and frequencies $w_v$ and stores the data for frequencies, and amount of trapped particles in a file. This function is called in our `main()`, and is used to create multiple simulations for a time-dependent system.

--------------------

### Linking and compiling of our project files:
To run the program `time_independent.cpp`:
```sh
g++ time_independent.cpp src/*.cpp -I include -o time_independent -O2 -llapack -lblas -larmadillo
```
```sh
./time_independent
```
To run the program `time_dependent.cpp`:
```sh
g++ time_dependent.cpp src/*.cpp -I include -o time_dependent -O2 -llapack -lblas -larmadillo
```
```sh
./time_dependent
```

----------------

The simulations ran by the programs `time_independent.cpp` and `time_dependent.cpp` produce sets of data, which can be imported and visualized by running the python programs with the similar names: 

### `time_independent.py`:
By running our program `time_independent.cpp`, we will obtain a number of files, which we can interpret and visualize by running the code in this python file. The program imports these files, and creates plots showing the particles movement in the xyz-plane, as well as relevant phase space plots, and an error analysis. 


### `time_dependent.py`:
By running our program `time_dependent.cpp`, we will obtain a number of files, which we can interpret and visualize by running the code in this python file. The program imports these files, and returns plots visualizing the fraction of trapped particles as a function of the frequency $w_v$. 


-------------

Finally, within our map `doc-construction`, one can find a compilation of all the plots created by `time_independent.py` and `time_dependent.py`, as well as the latex file that imports these plots into our Project report.
