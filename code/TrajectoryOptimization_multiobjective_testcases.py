import matplotlib.pyplot as plt
import TrajectoryOptimization as TrajectoryOptimization
import numpy as np
# from tikzplotlib import save as tikz_save

# This code calculates a multiobjective trajectory optimization problem and plots the resulting trajectory

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

print(accel_max, j_max)

trajectory = TrajectoryOptimization.TrajectoryGenerator(distance, dt, mass, [0, vel_max], [-accel_max, accel_max],
                                                        [-j_max, j_max], [-p_max, p_max], tf=final_time)

# Calculate minimum total energy trajectory
# ret, solProg = trajectory.GenerateTrajectory(var=['v', 'j'], norm=['abs', 'abs'], weights=[1.0, 0.00000001])
ret, solProg = trajectory.GenerateTrajectory(var=['p', 'p'], norm=['abs', 'peak'], weights=[dt*0.5, 0.5])  # Energy and peak power
d_t, v_t, a_t, j_t, P_t, dt = ret
time_t = [dt * i for i in range(0, len(d_t))]
E_e = sum(abs(P) * dt for P in P_t)
Ppk_e = max(abs(P) for P in P_t)
k = 8
t = [round(i * max(time_t) / k, 3) for i in range(0, k + 1)]
print('Total Energy (min time):', E_e)
print('Peak Power (min time):', Ppk_e)
plt.figure(5)
plt.subplot(511)
plt.plot(time_t, d_t)
plt.title('Minimum Mixed Objective')
plt.ylabel('Distance [m]')
plt.xlim([0, max(time_t)])
plt.xticks(t, [])
plt.grid()
plt.subplot(512)
plt.plot(time_t, v_t)
plt.ylabel('Velocity [m/s]')
plt.xlim([0, max(time_t)])
plt.xticks(t, [])
plt.grid()
plt.subplot(513)
plt.plot(time_t, a_t)
plt.ylabel('Acceleration [m/s^2]')
plt.xlim([0, max(time_t)])
plt.xticks(t, [])
plt.grid()
plt.subplot(514)
plt.plot(time_t, j_t)
plt.xlim([0, max(time_t)])
plt.xticks(t, [])
plt.ylabel('Jerk [m/s^3]')
plt.grid()
plt.subplot(515)
plt.plot(time_t, P_t)
plt.ylabel('Power [W]')
plt.xlabel('Time [s]')
plt.xlim([0, max(time_t)])
plt.xticks(t, t)
plt.grid()
# tikz_save('min_e.tex')

plt.show()
