import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

#a couple arguments that will be needed
B0 = 9.65e1; V0 = 2.41e6
d = 500; q = 1
m = 40.078
w_0 = q * B0 / m ; w_z = np.sqrt(2 * q * V0 / (m * d ** 2))

 #definition of movement in the z-direction as defined by formula (12)
def z_position(t, z0, vz0, w_z):
    return z0 * np.cos(w_z * t) + (vz0 / w_z) * np.sin(w_z * t)

time = np.linspace(1, 50, 4000) #time array with the same dimensions as the data simulated in c++
z0_particle1 = 20 ; vz0_particle1 = 0 #initial condition for particle 1
z_positions= z_position(time, z0_particle1, vz0_particle1, w_z) #calculates the position along the z-axis for particle 1 as a function of time

data = np.loadtxt("particle1motion.txt") #loading the data for the movement of particle 1 simulated by the c++ program time_independent.cpp
t = data[:, 0]; z= data[:, 3] #extracting data needed

plt.figure(figsize=(12, 6)) #we create plot
plt.plot(t, z, label="Particle 1")
plt.plot(time, z_positions, "--", label="Particle 1, analytical movement estimation")
plt.xlabel(r"t ($\mu$s)")
plt.ylabel(r"z ($\mu$m)")
plt.title("Particle 1 movement in z direction as a function of time")
plt.legend()
plt.grid()
plt.savefig("zdir_particle1motion.pdf")

data_no_interaction = np.loadtxt("withoutinteraction.txt") #imports data for two particles 1 and 2 in a penning system without consideration to particle interaction, simulated in the c++ program time_independent.cpp
no_t = data_no_interaction[:, 0] #we extract data we will need for our plots
no_x1 = data_no_interaction[:, 1]
no_y1 = data_no_interaction[:, 2]
no_z1 = data_no_interaction[:, 3]
no_vx1 = data_no_interaction[:, 4]
no_vy1 = data_no_interaction[:, 5]
no_vz1 = data_no_interaction[:, 6]
no_x2 = data_no_interaction[:, 7]
no_y2 = data_no_interaction[:, 8]
no_z2 = data_no_interaction[:, 9]
no_vx2 = data_no_interaction[:, 10]
no_vy2 = data_no_interaction[:, 11]
no_vz2 = data_no_interaction[:, 12]

data_with_interaction = np.loadtxt("withinteraction.txt") #imports data for two particles 1 and 2 in a penning system with consideration to particle interaction, simulated in the c++ program time_independent.cpp
int_t = data_with_interaction[:, 0] #we extract data we will need for our plots
int_x1 = data_with_interaction[:, 1]
int_y1 = data_with_interaction[:, 2]
int_z1 = data_with_interaction[:, 3]
int_vx1 = data_with_interaction[:, 4]
int_vy1 = data_with_interaction[:, 5]
int_vz1 = data_with_interaction[:, 6]
int_x2 = data_with_interaction[:, 7]
int_y2 = data_with_interaction[:, 8]
int_z2 = data_with_interaction[:, 9]
int_vx2 = data_with_interaction[:, 10]
int_vy2 = data_with_interaction[:, 11]
int_vz2 = data_with_interaction[:, 12]

plt.figure(figsize=(12, 6)) # we make a simple plot of the particles movement in xy-plane with and without consideration to particle interaction
plt.suptitle("Movement of particle 1 and 2 in the (x , y) -plane")
plt.subplot(1, 2, 1)
plt.title("No interaction between particles")
plt.plot(no_x1, no_y1, label="Particle 1")
plt.plot(no_x2, no_y2, label="Particle 2")
plt.xlabel(r"x ($\mu$s)")
plt.ylabel(r"y ($\mu$m)")
plt.legend()
plt.grid()
plt.subplot(1, 2, 2)
plt.title("With interaction between particles")
plt.plot(int_x1, int_y1, label="Particle 1")
plt.plot(int_x2, int_y2, label="Particle 2")
plt.xlabel(r"x ($\mu$s)")
plt.ylabel(r"y ($\mu$m)")
plt.legend()
plt.grid()
plt.savefig("xy_movement.pdf") #saves the files

fig = plt.figure(figsize=(12, 6))
plt.suptitle("Movement of particle 1 and 2 in the (x , y, t) space")# we make a  plot of the particles movement in xy-plane as a function of time, with and without consideration to particle interaction
ax1 = fig.add_subplot(121, projection="3d")
ax1.set_title("No interaction between particles")
ax1.plot(no_x1, no_y1, np.linspace(0, 50, len(no_x1)), label="Particle 1 ")
ax1.plot(no_x2, no_y2, np.linspace(0, 50, len(no_x2)), label="Particle 2 ")
ax1.set_xlabel(r"x ($\mu$m)")
ax1.set_ylabel(r"y ($\mu$m)")
ax1.set_zlabel(r"t ($\mu$s)")
ax1.legend()
ax2 = fig.add_subplot(122, projection="3d")
ax2.set_title("With interaction between particles")
ax2.plot(int_x1, int_y1, np.linspace(0, 50, len(int_x1)), label="Particle 1")
ax2.plot(int_x2, int_y2, np.linspace(0, 50, len(int_x2)), label="Particle 2")
ax2.set_xlabel(r"x ($\mu$m)")
ax2.set_ylabel(r"y ($\mu$m)")
ax2.set_zlabel(r"t ($\mu$s)")
ax2.legend()
plt.tight_layout()
plt.savefig("xyt_movement.pdf") #saves the plot as pdf file

########################TWO PARTICLES IN PENNING TRAP, TRAJECTORIES IN THE (X,VX) AND (Z,VZ) PLANES WITH AND WITHOUT PARTICLE INTERACTIONS
plt.figure(figsize=(12, 6))
plt.suptitle(r"Phase space plots for the particles 1 and 2 (x , $v_x$)") # we make a phase space plot of the particles (x,v_x), with and without consideration to particle interaction
plt.subplot(1,2,1)
plt.title("No interaction between particles")
plt.plot(no_x1, no_vx1, label="Particle 1")
plt.plot(no_x2, no_vx2, label="Particle 2")
plt.xlabel(r"x ($\mu$m)")
plt.ylabel(r"$v_x$ ($\mu$m / $\mu$s)")
plt.grid()
plt.legend()
plt.subplot(1,2,2)
plt.title("With interaction between particles")
plt.plot(int_x1, int_vx1, label="Particle 1")
plt.plot(int_x2, int_vx2, label="Particle 2")
plt.xlabel(r"x ($\mu$m)")
plt.ylabel(r"$v_x$ ($\mu$m / $\mu$s)")
plt.legend()
plt.grid()
plt.savefig("Phasespace_(x,v_x).pdf") #saves plot as pdf file

fig1 = plt.figure(figsize=(12, 6))
fig1.suptitle(r"Phase space plots for the particles 1 and 2 (x, $v_x$ , t)")# we make a phase space plot of the particles (x,v_x,t), with and without consideration to particle interaction
ax1 = fig1.add_subplot(121, projection="3d")
ax1.set_xlabel(r"x ($\mu$m)")
ax1.set_ylabel(r"$v_x$ ($\mu$m / $\mu$s)")
ax1.set_zlabel(r"t ($\mu$s)")
ax1.set_title("No interaction between particles")
ax1.plot(no_x1, no_vx1, np.linspace(0, 50, len(no_x1)), label="Particle 1")
ax1.plot(no_x2, no_vx2, np.linspace(0, 50, len(no_x2)), label="Particle 2")
ax1.legend()
ax2 = fig1.add_subplot(122, projection="3d")
ax2.set_xlabel(r"x ($\mu$m)")
ax2.set_ylabel(r"$v_x$ ($\mu$m / $\mu$s)")
ax2.set_zlabel(r"t ($\mu$s)")
ax2.set_title("With interaction between particles")
ax2.plot(int_x1, int_vx1, np.linspace(0, 50, len(int_x1)), label="Particle 1")
ax2.plot(int_x2, int_vx2, np.linspace(0, 50, len(int_x2)), label="Particle 2")
ax2.legend()
plt.tight_layout()
plt.savefig("Phasespace_(x,v_x,t).pdf")

plt.figure(figsize=(12, 6))
plt.suptitle(r"Phase space plots for the particles 1 and 2 (z , $v_z$)")# we make a phase space plot of the particles (z,v_z), with and without consideration to particle interaction
plt.subplot(1,2,1)
plt.title("No interaction between particles")
plt.plot(no_z1, no_vz1, label="Particle 1")
plt.plot(no_z2, no_vz2, label="Particle 2")
plt.xlabel(r"z ($\mu$m)")
plt.ylabel(r"$v_z$ ($\mu$m / $\mu$s)")
plt.grid()
plt.legend()
plt.subplot(1,2,2)
plt.title("With interaction between particles")
plt.plot(int_z1, int_vz1, label="Particle 1")
plt.plot(int_z2, int_vz2, label="Particle 2")
plt.xlabel(r"z ($\mu$m)")
plt.ylabel(r"$v_z$ ($\mu$m / $\mu$s)")
plt.legend()
plt.grid()
plt.savefig("Phasespace_(z,v_z).pdf")

fig2 = plt.figure(figsize=(12, 6))
fig2.suptitle(r"Phase space plots for the particles 1 and 2 (z, $v_z$ , t)")# we make a phase space plot of the particles (z,v_z,t), with and without consideration to particle interaction
ax3 = fig2.add_subplot(121, projection="3d")
ax3.set_xlabel(r"z ($\mu$m)")
ax3.set_ylabel(r"$v_z$ ($\mu$m / $\mu$s)")
ax3.set_zlabel(r"t ($\mu$s)")
ax3.set_title("No interaction between particles")
ax3.plot(no_z1, no_vz1, np.linspace(0, 50, len(no_z1)), label="Particle 1")
ax3.plot(no_z2, no_vz2, np.linspace(0, 50, len(no_z2)), label="Particle 2")
ax3.legend()
ax4 = fig2.add_subplot(122, projection="3d")
ax4.set_xlabel(r"z ($\mu$m)")
ax4.set_ylabel(r"$v_z$ ($\mu$m / $\mu$s)")
ax4.set_zlabel(r"t ($\mu$s)")
ax4.set_title("With interaction between particles")
ax4.plot(int_z1, int_vz1, np.linspace(0, 50, len(int_z1)), label="Particle 1")
ax4.plot(int_z2, int_vz2, np.linspace(0, 50, len(int_z2)), label="Particle 2")
ax4.legend()
plt.tight_layout()
plt.savefig("Phasespace_(z,v_z,t).pdf")

########################3D PLOT OF TWO PARTICLES IN THE (X,Y,Z) SPACE WITH AND WITHOUT INTERACTIONS
plt.figure(figsize=(12, 6))
ax = plt.axes(projection="3d")
ax.set_title("Movement of the particles 1 and 2 with interactions")# we make a plot of particle movement in the xyz-plane, with and without particle interactions
ax.plot3D(int_x1, int_y1, int_z1, label="Particle 1")
ax.plot3D(int_x2, int_y2, int_z2, label="Particle 2")
ax.set_xlabel(r"x ($\mu$m)")
ax.set_ylabel(r"y ($\mu$m)")
ax.set_zlabel(r"z ($\mu$m)")
ax.legend()
plt.tight_layout()
plt.savefig(r"3Dwithinteractions.pdf")

plt.figure(figsize=(12, 6))
ax = plt.axes(projection="3d")
ax.set_title("Movement of the particles 1 and 2 without interactions")
ax.plot3D(no_x1, no_y1, no_z1, label="Particle 1")
ax.plot3D(no_x2, no_y2, no_z2, label="Particle 2")
ax.set_xlabel(r"x ($\mu$m)")
ax.set_ylabel(r"y ($\mu$m)")
ax.set_zlabel(r"z ($\mu$m)")
ax.legend()
plt.tight_layout()
plt.savefig(r"3Dwithoutinteractions.pdf")

########################PLOTTING SIZE OF RELATIVE ERROR AT EACH TIMESTEP
x0 = 20; y0 = 0; z0 = 20 #initial condition of particle 1
v0 = 25 #this is vy0, initial velocity in the y-direction
w_p = (w_0 + np.sqrt(w_0 ** 2 - 2 * w_z ** 2)) / 2 #definition from equation (16) in the report
w_m = (w_0 - np.sqrt(w_0 ** 2 - 2 * w_z ** 2)) / 2
A_p = (v0 + w_m * x0) / (w_m - w_p) #these are definitions as result of the 2.3 special case of the initial conditions of particle 1
A_m = -(v0 + w_p * x0) / (w_m - w_p)
p_p = 0 ; p_m = 0

def xy_positions_analytical(t): #define the expected movement of particle 1 by formulas (17)
    x = A_p * np.cos(w_p * t + p_p) + A_m * np.cos(w_m * t + p_m)
    y = -A_p * np.sin(w_p * t + p_p) - A_m * np.sin(w_m * t + p_m)
    return x, y

def z_position_with_velocity(t): #redefined function to only apply to particle 1, since v0z=0
    return z0 * np.cos(w_z * t)

rel_rk4 = [] ; rel_eul = [] #empty lists to contain relative errors
abs_rk4 = [] ; abs_euler = [] #empty lists to contain the largest absolute error for each file
step_sizes = [4000, 8000, 16000, 32000] #amount of steps used
file_names = ["qRK4", "qeuler"] #first name of each file, as stored by c++ code
for step_size in step_sizes: #loops through the step_size list
    for name in file_names: #loops through the file_names list
        filename = f"{name}{step_size}.txt" #constructs the name of each file we need to read, loads
        data = np.loadtxt(filename)
        if name == "qRK4":
            data_rk4 = data
            method_name = "RK4"
        else:
            data_euler = data
            method_name = "Euler"
        t = data[:, 0]
        x1 = data[:, 1]
        y1 = data[:, 2]
        z1 = data[:, 3]
        x_analytical, y_analytical = xy_positions_analytical(t) #computes analytical position as function of time
        z_analytical = z_position_with_velocity(t)
        abs_errors = np.sqrt((x_analytical - x1) ** 2 + (y_analytical - y1) ** 2 + (z_analytical - z1) ** 2) #uses formula (18) to compute absolute error
        rel_errors = abs_errors / np.sqrt(x_analytical ** 2 + y_analytical ** 2 + z_analytical ** 2) #uses formula (19) to calculare relative error
        if name == "qRK4": #stores data in correct list, based on what method was used to simulate the data
            rel_rk4.append(rel_errors)
            abs_rk4.append(np.max(abs_errors)) #checks what the largest absolute error was and adds it to the correct list
        else:
            rel_eul.append(rel_errors)
            abs_euler.append(np.max(abs_errors))

con_rk4 = 0; con_eul = 0
for i in range(1, len(abs_rk4)):
    con_rk4 += (1/3) * np.log(abs_rk4[i-1]/abs_rk4[i]) / np.log(step_sizes[i-1]/step_sizes[i]) #formula (20) for error convergence
for i in range(1, len(abs_euler)):
    con_eul += (1/3) * np.log(abs_euler[i-1]/abs_euler[i]) / np.log(step_sizes[i-1]/step_sizes[i]) #formula (20) for error convergence

#we plot the data, and add the convergence of relative error as the title of each subplot
plt.figure(figsize=(12, 6))
plt.suptitle("Relative errors for Particle 1")
plt.subplot(1,2,1)
for i in range(4):
    t=np.linspace(0,50,step_sizes[i])
    plt.plot(t,(rel_rk4[i]),label=f"n={step_sizes[i]}")
plt.xlabel(r"t ($\mu$s)")
plt.ylabel("Relative Error")
plt.title(f"RK4-method, with an error convergence rate {abs(con_rk4):.3f}")
plt.legend()
plt.grid(True)  
plt.subplot(1,2,2)
for i in range(4):
    t=np.linspace(0,50,step_sizes[i])
    plt.plot(t,(rel_eul[i]),label=f"n={step_sizes[i]}")
plt.xlabel(r"t ($\mu$s)")
plt.ylabel("Relative Error")
plt.title(f"Forward-Euler-method, with an error convergence rate {abs(con_eul):.3f}")
plt.legend()
plt.grid(True)
plt.savefig(r"Relativeerror.pdf")

#We create the same plot, with logaritmix y-axis
plt.figure(figsize=(12, 6))
plt.suptitle("Relative errors for Particle 1, scaled")
plt.subplot(1,2,1)
for i, step_size in enumerate(step_sizes):
    t=np.linspace(0,50,step_sizes[i])
    plt.plot(t,(rel_rk4[i]),label=f"n={step_sizes[i]}")
plt.xlabel(r"t ($\mu$s)")
plt.ylabel(r"$log_{10}$(Relative Error)")
plt.yscale("log") 
plt.title(f"RK4-method, with an error convergence rate {abs(con_rk4):.3f}")
plt.legend()
plt.grid(True)  
plt.subplot(1,2,2)
for i, step_size in enumerate(step_sizes):
    t=np.linspace(0,50,step_sizes[i])
    plt.plot(t,(rel_eul[i]),label=f"n={step_sizes[i]}")
plt.xlabel(r"t ($\mu$s)")
plt.ylabel(r"$log_{10}$(Relative Error)")
plt.yscale("log")
plt.title(f"Forward-Euler-method, with an error convergence rate {abs(con_eul):.3f}")
plt.legend()
plt.grid(True)
plt.savefig(r"Relativeerrorformatted.pdf")