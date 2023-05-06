#
# Copyright (c) 2021, ATA Engineering, Inc.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import matplotlib.pyplot as plt
# import TrajectoryOptimization_orig as TrajectoryOptimization
import TrajectoryOptimization as TrajectoryOptimization
import numpy as np
# from tikzplotlib import save as tikz_save

# This code solves trajectory optimization problems for a variety of different objective functions
# It produces the figures presented in the trajectory optimization paper


def solve_and_plot_trajectories(var, time, n, title=''):
    ret1, solProg1 = trajectory.GenerateTrajectory(var=var, norm=['abs'])
    d_1, v_1, a_1, j_1, P_1, dt_1 = ret1
    ret2, solProg2 = trajectory.GenerateTrajectory(var=var, norm=['rms'])
    d_2, v_2, a_2, j_2, P_2, dt_2 = ret2
    ret_inf, solProginf = trajectory.GenerateTrajectory(var=var, norm=['peak'])
    d_inf, v_inf, a_inf, j_inf, P_inf, dt_inf = ret_inf

    tick_spacing = 8
    t_ticks = [round(i * max(time) / tick_spacing, 3) for i in range(0, tick_spacing + 1)]
    plt.figure(n)
    plt.subplot(511)
    plt.plot(time, d_1, time, d_2, time, d_inf)
    plt.title(title)
    plt.ylabel('Distance [m]')
    plt.xlim([0, max(time)])
    plt.xticks(t_ticks, [])
    plt.grid()
    plt.legend([var[0] + '_1', var[0] + '_2', var[0] + '_inf'])
    plt.subplot(512)
    plt.plot(time, v_1, time, v_2, time, v_inf)
    plt.ylabel('Velocity [m/s]')
    plt.xlim([0, max(time)])
    plt.xticks(t_ticks, [])
    plt.grid()
    plt.subplot(513)
    plt.plot(time, a_1, time, a_2, time, a_inf)
    plt.ylabel('Acceleration [m/s^2]')
    plt.xlim([0, max(time)])
    plt.xticks(t_ticks, [])
    plt.grid()
    plt.subplot(514)
    plt.plot(time, j_1, time, j_2, time, j_inf)
    plt.xlim([0, max(time)])
    plt.xticks(t_ticks, [])
    plt.ylabel('Jerk [m/s^3]')
    plt.grid()
    plt.subplot(515)
    plt.plot(time, P_1, time, P_2, time, P_inf)
    plt.ylabel('Power [W]')
    plt.xlabel('Time [s]')
    plt.xlim([0, max(time)])
    plt.xticks(t_ticks, t_ticks)
    plt.grid()
    E_1 = sum(abs(P) * dt for P in P_1)
    E_2 = sum(abs(P) * dt for P in P_2)
    E_inf = sum(abs(P) * dt for P in P_inf)
    Ppk_1 = max(abs(P) for P in P_1)
    Ppk_2 = max(abs(P) for P in P_2)
    Ppk_inf = max(abs(P) for P in P_inf)
    print('Total Energy (1-norm, ' + var[0] + '):', E_1)
    print('Total Energy (2-norm, ' + var[0] + '):', E_2)
    print('Total Energy (inf-norm, ' + var[0] + '):', E_inf)
    print('Peak Power (1-norm, ' + var[0] + '):', Ppk_1)
    print('Peak Power (2-norm, ' + var[0] + '):', Ppk_2)
    print('Peak Power (inf-norm, ' + var[0] + '):', Ppk_inf)

    return [E_1, E_2, E_inf], [Ppk_1, Ppk_2, Ppk_inf]


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

time = trajectory.time
E_v, Ppk_v = solve_and_plot_trajectories(['v'], time, 1, 'Minimum Velocity')
# tikz_save('min_v_high_j.tex')
E_a, Ppk_a = solve_and_plot_trajectories(['a'], time, 2, 'Minimum Acceleration')
# tikz_save('min_a.tex')
E_j, Ppk_j = solve_and_plot_trajectories(['j'], time, 3, 'Minimum Jerk')
# tikz_save('min_j.tex')
E_p, Ppk_p = solve_and_plot_trajectories(['p'], time, 4, 'Minimum Power')
# tikz_save('min_p.tex')

# Calculate minimum total energy trajectory
ret, solProg = trajectory.GenerateTrajectory(var=['e'], norm=['abs'])
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
plt.title('Minimum Total Energy')
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

# Calculate minimum time trajectory
ret, solProg = trajectory.GenerateTrajectory(var=['t'], norm=['abs'])
d_t, v_t, a_t, j_t, P_t, dt = ret
time_t = [dt * i for i in range(0, len(d_t))]
E_t = sum(abs(P) * dt for P in P_t)
Ppk_t = max(abs(P) for P in P_t)
k = 8
t = [round(i * max(time_t) / k, 3) for i in range(0, k + 1)]
print('Total Energy (min time):', E_t)
print('Peak Power (min time):', Ppk_t)
plt.figure(6)
plt.subplot(511)
plt.plot(time_t, d_t)
plt.title('Minimum Time')
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
# tikz_save('min_dt.tex')

# Compute and print tables
E_values = np.array([E_v, E_a, E_j, E_p])
E_max = E_values.max()
E_min = E_values.min()
E_max_percent = E_values/E_max
E_min_percent = E_values/E_min
print("Energy tables")
print(E_values.T)
print(E_e)
print(E_t)
print(E_max_percent.T)
print(E_e/E_max)
print(E_t/E_max)
print(E_min_percent.T)
print(E_e/E_min)
print(E_t/E_min)

print("Peak power tables")
Ppk_values = np.array([Ppk_v, Ppk_a, Ppk_j, Ppk_p])
Ppk_max = Ppk_values.max()
Ppk_min = Ppk_values.min()
Ppk_max_percent = Ppk_values/Ppk_max
Ppk_min_percent = Ppk_values/Ppk_min
print(Ppk_values.T)
print(Ppk_e)
print(Ppk_t)
print(Ppk_max_percent.T)
print(Ppk_e/Ppk_max)
print(Ppk_t/Ppk_max)
print(Ppk_min_percent.T)
print(Ppk_e/Ppk_min)
print(Ppk_t/Ppk_min)

plt.show()

# input("Press enter to close")
