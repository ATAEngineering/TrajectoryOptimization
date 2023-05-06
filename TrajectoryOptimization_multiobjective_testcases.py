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
from tikzplotlib import save as tikz_save

import TrajectoryOptimization as TrajectoryOptimization

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

# Calculate minimum total energy trajectory and peak power
ret, solProg = trajectory.GenerateTrajectory(var=['p', 'p'], norm=['abs', 'peak'], weights=[dt*0.5, 0.5])
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
tikz_save('multiobj_e_p.tex')

plt.show()
