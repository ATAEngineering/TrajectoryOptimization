import matplotlib.pyplot as plt
import TrajectoryOptimization_2jointarm as TrajectoryOptimization
import math
import numpy as np
from tikzplotlib import save as tikz_save

def make_plots(ret):
    dt, theta1, theta2, w1, w2, z1, z2, j1, j2, t1, t2, p1, p2, tip1PosX, tip1PosY, tip1VelX, tip1VelY, tip2PosX, tip2PosY, tip2VelX, tip2VelY = ret
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
                                                        theta1i, theta2i, theta1f, theta2f, w1f, w2f, tipVelXf, tipVelYf,
                                                        tipmass, arm1mass, arm2mass, J1, J2,
                                                        L1, L2,
                                                        [theta1_min, theta1_max], [theta2_min, theta2_max],
                                                        [-w1_max, w1_max], [-w2_max, w2_max],
                                                        [-z1_max, z1_max], [-z2_max, z2_max],
                                                        [-jerk1_max, jerk1_max], [-jerk2_max, jerk2_max], tf=final_time)

# Calculate minimum total energy trajectory
ret, solProg = trajectory.GenerateTrajectory(var=['p1', 'p2'], norm=['peak', 'peak'], weights=[0.5, 0.5])
# ret, solProg = trajectory.GenerateTrajectory(var=['v', 'j'], norm=['abs', 'abs'], weights=[1.0, 0.00000001])
# ret, solProg = trajectory.GenerateTrajectory(var=['p', 'p'], norm=['abs', 'peak'], weights=[dt*0.5, 0.5])  # Energy and peak power

make_plots(ret)


