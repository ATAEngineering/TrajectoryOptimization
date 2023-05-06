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
import math
from tikzplotlib import save as tikz_save

import TrajectoryOptimization_2jointarm as TrajectoryOptimization


def make_plots(ret):
    dt, theta1, theta2, w1, w2, z1, z2, j1, j2, t1, t2, p1, p2,\
        tip1PosX, tip1PosY, tip1VelX, tip1VelY, tip2PosX, tip2PosY, tip2VelX, tip2VelY = ret
    time_t = [dt * i for i in range(0, len(theta1))]
    k = 8
    t = [round(i * max(time_t) / k, 3) for i in range(0, k + 1)]

    plt.figure(1)
    plt.subplot(611)
    plt.title('Minimum Mixed Objective - Joint 1')
    plt.plot(time_t, theta1)
    plt.ylabel('Theta_1 [rad]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.subplot(612)
    plt.plot(time_t, w1)
    plt.ylabel('w_1 [rad/s]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.subplot(613)
    plt.plot(time_t, z1)
    plt.ylabel('z_1 [rad/s^2]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.subplot(614)
    plt.plot(time_t, j1)
    plt.ylabel('j_1 [rad/s^3]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.subplot(615)
    plt.plot(time_t, t1)
    plt.ylabel('t_1 [N-m]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.subplot(616)
    plt.plot(time_t, p1)
    plt.ylabel('p_1 [W]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.xlabel('Time [s]')
    tikz_save('../paper/figures/2link_trajectory_m1.tex')

    plt.figure(2)
    plt.subplot(611)
    plt.title('Minimum Mixed Objective - Joint 2')
    plt.plot(time_t, theta2)
    plt.ylabel('Theta_2 [rad]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.subplot(612)
    plt.plot(time_t, w2)
    plt.ylabel('w_2 [rad/s]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.subplot(613)
    plt.plot(time_t, z2)
    plt.ylabel('z_2 [rad/s^2]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.subplot(614)
    plt.plot(time_t, j2)
    plt.ylabel('j_2 [rad/s^3]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.subplot(615)
    plt.plot(time_t, t2)
    plt.ylabel('t_2 [N-m]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.subplot(616)
    plt.plot(time_t, p2)
    plt.ylabel('p_2 [W]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.xlabel('Time [s]')
    tikz_save('../paper/figures/2link_trajectory_m2.tex')

    plt.figure(3)
    plt.subplot(211)
    plt.title('Minimum Mixed Objective - Tip')
    plt.plot(time_t, tip1PosX, time_t, tip1PosY)
    plt.plot(time_t, tip2PosX, time_t, tip2PosY)
    plt.ylabel('Tip position [m]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.subplot(212)
    plt.plot(time_t, tip1VelX, time_t, tip1VelY)
    plt.plot(time_t, tip2VelX, time_t, tip2VelY)
    plt.ylabel('Tip speed [m/s]')
    plt.xlim([0, max(time_t)])
    plt.xticks(t, [])
    plt.grid()
    plt.xlabel('Time [s]')
    tikz_save('../paper/figures/2link_tip_position_history.tex')

    plt.figure(4)
    plt.title('X-Y Tip Position')
    plt.plot(tip1PosX, tip1PosY)
    plt.plot(tip2PosX, tip2PosY)
    plt.xlabel('Tip X Postion[m]')
    plt.ylabel('Tip Y Position [m]')
    tikz_save('../paper/figures/2link_tip_position_2d.tex')
    plt.grid()

    plt.show()


# This code calculates a multiobjective trajectory optimization problem and plots the resulting trajectory

dt = 0.005  # Sampling frequency, [s]

# Boundary conditions
final_time = 1.0  # Final time of the move, [s]
theta1i = -3.14/2  # Initial angle of joint 1, [rad]
theta2i = -3.14  # Initial angle of joint 2, [rad]
theta1f = 0.0  # Final angle of joint 1, [rad]
theta2f = -3.14/2  # Final angle of joint 2, [rad]
w1f = 5  # Final speed of joint 1, [rad/s]
w2f = 5  # Final speed of joint 2, [rad/s]
tipVelXf = 0.0  # Final tip X speed [m/s]
tipVelYf = 13.0  # Final tip Y speed [m/s]

# Parameters
tipmass = 0.5  # Tip payload mass, [kg]
arm1mass = 3  # Arm1 mass, [kg]
arm2mass = 2  # Arm2 mass, [kg]
J1 = 1.0  # Arm 1 inertia [kg-m^2]
J2 = 1.0  # Arm 2 inertia [kg-m^2]
L1 = 0.6  # Arm 1 length [m]
L2 = 0.4  # Arm 2 length [m]

# Constraints
theta1_min = -math.pi  # Maximum angle, [rad]
theta1_max = math.pi  # Maximum angle, [rad]
theta2_min = -math.pi  # Maximum angle, [rad]
theta2_max = math.pi  # Maximum angle, [rad]
w1_max = 10.0  # Maximum angular acceleration, [rad/s]
w2_max = 10.0  # Maximum angular acceleration, [rad/s]
z1_max = 100.0  # Maximum angular acceleration, [rad/s^2]
z2_max = 100.0  # Maximum angular acceleration, [rad/s^2]
jerk1_max = 1000.0  # Maximum angular jerk, [rad/s^3]
jerk2_max = 1000.0  # Maximum angular jerk, [rad/s^3]

trajectory = TrajectoryOptimization.TrajectoryGenerator(dt,
                                                        theta1i, theta2i, theta1f, theta2f, w1f, w2f,
                                                        tipVelXf, tipVelYf,
                                                        tipmass, arm1mass, arm2mass, J1, J2,
                                                        L1, L2,
                                                        [theta1_min, theta1_max], [theta2_min, theta2_max],
                                                        [-w1_max, w1_max], [-w2_max, w2_max],
                                                        [-z1_max, z1_max], [-z2_max, z2_max],
                                                        [-jerk1_max, jerk1_max], [-jerk2_max, jerk2_max],
                                                        tf=final_time)

# Calculate minimum total energy trajectory
ret, solProg = trajectory.GenerateTrajectory(var=['p1', 'p2'], norm=['peak', 'peak'], weights=[0.5, 0.5])

make_plots(ret)
