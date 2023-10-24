#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>

#ifndef __penningtrap_class_h__  
#define __penningtrap_class_h__


class PenningTrap{
  private: 
  double ke=1.389e5; //we define ke as a private argument that cannot be edited
  public:
  double B0 = 96.5; //we define arguments needed to 
  double V0 = 2.41*pow(10,6); 
  double d = 500.0;
  std::vector<Particle> particle_vector; //particles within our system
  std::vector<bool> particles_escaped; //vector that will contain information about whether a particle exists within the system or not 
  //arguments that changes how we model the environment within the trap will act:
  bool time_dependence=false; //we set the norm to correspond to no consideration to particle interaction
  bool Particle_interactions = false;//an argumet for particle interaction, that can be changed for an arbitrary penningtrap, set to false by default
  //variables needed to the time-dependent model the system containing many particles:
  int particle_amount = 0; //this is an integer that we change only for simulation of particle_amount Ca particles with randomly generated initial conditions
  double f_amp = 0.0; // we define two variables f_amp and w_v that has to do with modelling our time-dependence
  double w_v = 0.0; //these will be changed in the instance when we wish to model a system with time-dependence

  void particle_add(Particle p_in){//this function adds a particle into our system
    particle_vector.push_back(p_in);
    particle_amount += 1;
    particles_escaped.push_back(p_in.escape_test(d));
  }

  int count_particles() const{//this function counts the amount of particles that exists within our system
    int count = 0;
    for (int i = 0; i < particle_vector.size(); i++){
      if (!particles_escaped[i]){//test that checks if a particle has escaped. if not, it is counted
        const arma::vec& r = particle_vector[i].r;
        double distance = arma::norm(r);
        if (distance <= d){
          count++;
        }
      }
    }
    return count; //returns the amount of particles that exists within the system 
  }

  arma::vec ext_E_field(arma::vec r,double t){//model of the electric field 
    arma::vec E = arma::vec(3);
    double r_norm = arma::norm(r);
    if (time_dependence == true){//checks if we want to model by using time dependence
      double V_def = V0*(1.0+f_amp*cos(w_v*t)); //the definition for V0 for time dependence, formula (8) in report
      if (r_norm <= d) {//definition of electric field, formula (6) in the report
        E(0) = V_def/pow(d,2)*r(0);
        E(1) = V_def/pow(d,2)*r(1);     
        E(2) = -2.0*V_def/pow(d,2)*r(2);
      }
    }
    else{ //if no time dependence
      if (r_norm <= d) {//definition of electric field, formula (6) in the report
        E(0) = V0/pow(d,2)*r(0);
        E(1) = V0/pow(d,2)*r(1);            
        E(2) = -2.0*V0/pow(d,2)*r(2);
      }
    }
    return E; //returns the electric field
  }

  arma::vec ext_B_field(const arma::vec& r){//model of the magnetic field
    arma::vec B(3); // using the definition of B-field from theory
    double r_norm = arma::norm(r);
    if (r_norm <= d) {
      B(0) = 0;
      B(1) = 0;
      B(2) = B0; //only one component in z-direction
    }
    return B; //returns the B-field vector
  }

  arma::vec force_particle(int i, int j) {//calculates the force on a particle, from 
    if (Particle_interactions == false) {//checks if there are muntiple particles within the system
        return arma::zeros(3); //if not, returns a vector containing only zeros
    }
    else{
      const arma::vec r_vec = particle_vector[i].r - particle_vector[j].r;//distance between the particles
      const double abs_r_vec = arma::norm(r_vec); //distance between the particles
      const arma::vec E_ = ke * particle_vector[j].q * (r_vec / (pow(abs_r_vec, 3))); //coloumb's law, formula (9) in the report
      return particle_vector[i].q * E_; //returns the lorentz force vector, affecting our particle
    }
  }

  arma::vec total_force_external(int i, double t) {
    if (particles_escaped[i]) {//checks if the particle has escaped
      return arma::zeros(3); //if escaped, returns a vector containing only zeros
    }
    const arma::vec ri = particle_vector[i].r;
    const arma::vec vi = particle_vector[i].v;
    const double qi = particle_vector[i].q;
    const arma::vec external_force_Efield = qi * ext_E_field(ri, t); //calculation of force from electric field
    const arma::vec external_force_Bfield = {qi*vi(1)*B0,-qi*vi(0)*B0,0.0}; //calculation of force from magnetic field
    return external_force_Bfield + external_force_Efield; // returns a calculation of total force from the environment 
  }
  
  arma::vec total_force_particles(int i){//calculation of all the forces working between the particles within the system
    arma::vec particles_force(3, arma::fill::zeros);
    for (int j=0; j < particle_vector.size(); j++){
      if (particles_escaped[j] == false){
        if (j != i){
            particles_force += force_particle(i, j); //simply sums the forces on particle i from all other particles within the system
        }
      }
    }
  return particles_force; 
  }

  arma::vec total_force(int i, double t) {//calculates the sum of external forces, and forces between particles in the system
    arma::vec external_force = total_force_external(i,t); //call our functions 
    arma::vec particles_force = total_force_particles(i);//call our functions 
    arma::vec total_force = external_force + particles_force; //sums
    return total_force; //returns the total force vector
  }
  
  void evolve_forward_Euler(double dt, double t) {//implementation of forward-euler-method, using the definition presented in 2.4.2
    for (int i = 0; i < particle_vector.size(); i++) {
      if (particles_escaped[i] == false){//checks if the particle has escaped
        arma::vec total_forces = total_force(i,t); // Use the existing total_force function to calculate force
        particle_vector[i].v += (total_forces / particle_vector[i].m) * dt; //uses a=F/m to calculate new v
        particle_vector[i].r += particle_vector[i].v * dt; //uses r=v*dt to calculate new position
        particles_escaped[i] = particle_vector[i].escape_test(d); //checks if the particle escapes
      }
    }
  }

  void evolve_RK4(double dt, double t) {//implementation of runge-kutta defined in 2.4.1 in the report
    std::vector<Particle> particles_ = particle_vector; //creates a copy of the position-vector
    for (int i = 0; i < particle_vector.size(); i++) {//for each dimention, xyz, the operation is done
      arma::vec k_r1 = particle_vector[i].v * dt;
      arma::vec k_v1 = (total_force(i, t) / particle_vector[i].m) * dt; //uses a=F/m to estimate v
      particle_vector[i].r = particles_[i].r + k_r1 / 2.0;
      particle_vector[i].v = particles_[i].v + k_v1 / 2.0;

      arma::vec k_r2 = particle_vector[i].v * dt;
      arma::vec k_v2 = (total_force(i, t + dt / 2.0) / particle_vector[i].m) * dt;
      particle_vector[i].r = particles_[i].r + k_r2 / 2.0;
      particle_vector[i].v = particles_[i].v + k_v2 / 2.0;

      arma::vec k_r3 = particle_vector[i].v * dt;
      arma::vec k_v3 = (total_force(i, t + dt / 2.0) / particle_vector[i].m) * dt;
      particle_vector[i].r = particles_[i].r + k_r3;
      particle_vector[i].v = particles_[i].v + k_v3;

      arma::vec k_r4 = particle_vector[i].v * dt;
      arma::vec k_v4 = (total_force(i, t + dt) / particle_vector[i].m) * dt;
      particle_vector[i].r = particles_[i].r + (k_r1 + 2.0 * k_r2 + 2.0 * k_r3 + k_r4) / 6.0;
      particle_vector[i].v = particles_[i].v + (k_v1 + 2.0 * k_v2 + 2.0 * k_v3 + k_v4) / 6.0;
      particles_escaped[i] = particle_vector[i].escape_test(d); //checks if the particle escapes
    }
  }
};

#endif