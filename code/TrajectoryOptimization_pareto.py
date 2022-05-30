import matplotlib.pyplot as plt
import TrajectoryOptimization
import numpy as np
from tikzplotlib import save as tikz_save

# This code computes the Pareto front for a multiobjective optimization problem

mass = 10.0  # Load mass, [kg]
dt = 0.005  # Sampling frequency, [s]
distance = 0.254  # Move Distance, [m]
final_time = 1.0  # Final time of the move, [s]
vel_max = 0.5  # Maximum velocity, [m/s]
accel_max = 9.81 / 2  # Maximum acceleration, [m/s^2]
# m = 1               # Maximum jerk scaling factor
# j_max = accel_max/(m*dt)
j_max = 98.1
p_max = 100.0

trajectory = TrajectoryOptimization.TrajectoryGenerator(distance, dt, mass, [0, vel_max], [-accel_max, accel_max],
                                                        [-j_max, j_max], [-p_max, p_max], tf=final_time)

weightrange = np.linspace(0, 1, num=51)
E = []
P = []

# Loop over weights
for weight in weightrange:

    # Calculate weighted minimum total energy trajectory and peak power
    ret, solProg = trajectory.GenerateTrajectory(var=['p', 'p'], norm=['abs', 'peak'], weights=[dt*weight, 1-weight])  # Energy and peak power
    d_t, v_t, a_t, j_t, P_t, dt = ret

    E_e = sum(abs(P) * dt for P in P_t)
    Ppk_e = max(abs(P) for P in P_t)
    print('weight_e, energy, peak power:', weight, E_e, Ppk_e)

    E.append(E_e)
    P.append(Ppk_e)

plt.figure(8)
plt.title('Minimum Total Energy-Peak Power Pareto Front')
plt.xlabel('Total Energy [J]')
plt.ylabel('Peak Power [W]')
plt.scatter(E, P)

plt.grid()
# tikz_save('pareto_ep.tex')

# Plot normalized to minimums
bestE = min(E)
bestP = min(P)
Eratio = [e/bestE for e in E]
Pratio = [p/bestP for p in P]

plt.figure(9)
plt.title('Minimum Total Energy-Peak Power Pareto Front')
plt.xlabel('Normalized Total Energy [J]')
plt.ylabel('Normalized Peak Power [W]')
plt.scatter(Eratio, Pratio)

plt.grid()
tikz_save('pareto_ep.tex')

plt.show()
