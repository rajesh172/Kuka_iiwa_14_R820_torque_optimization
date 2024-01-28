# Kuka_iiwa_14_R820_torque_optimization
# Project Overview

This GitHub repository encapsulates a thorough exploration and implementation of motion control and torque optimization for the Kuka LBR iiwa 14 R820 robot, featuring a 7-degree-of-freedom manipulator. Executed in two distinct phases, the project comprises an extensive analytical analysis and practical simulations to assess and enhance the manipulator's operational capabilities.

## Phase 1: Analytical Analysis

### Task Representation
We initiated the project by defining the manipulator's workspace, a crucial step in understanding its reach and capabilities.

### Position Analysis
An in-depth kinematic analysis was conducted to determine the manipulator's precise position, orientation, and workspace limits.

### Velocity and Statics Analysis
Exploration of the manipulator's velocity and statics properties was undertaken, providing critical insights for controlling movements and ensuring stability.

### Redundancy Resolution
Given the seven degrees of freedom, redundancy resolution was addressed to optimize movements and task execution efficiently.

### Stiffness Analysis
The stiffness of the robotic system was assessed, offering valuable insights into its compliance and ability to interact with objects delicately.

We utilized tools such as Matlab and CopeliaSim for these analyses, combining symbolic computation, numerical analysis, and graphical representation to understand the Kuka robot better.

## Phase 2: Motion Control and Torque Optimization

Project 2 focused on motion control and torque optimization for the joints of the Kuka LBR iiwa 14 R820 robot. The specific challenge was to simulate a practical scenario where the robot carries a glass of water, applying a constant force in the -z direction due to the weight of the glass. The objective was to optimize the joint space, allowing the robot to follow a predefined trajectory with a constant force/torque at the tip while minimizing the net joint torque.

### Circular Trajectory Demonstration
To illustrate the manipulator's ability to hold a glass of water securely, we conducted a circular trajectory demonstration in a plane, simulating real-world scenarios.

### Glass Transportation Task
The robot was tasked with transporting a glass from one point to another, requiring the end-effector to follow a specific trajectory while maintaining a constant force at its tip. This involved supplying joint torques focused on minimizing energy consumption and execution difficulty.

### Null Space Contribution
To address the challenge of joint torque optimization, we leveraged the Null space contribution of the robot. The formulation of the Null-space projection vector was strategically designed to optimize joint torques throughout the circular trajectory, considering both joint configuration and force/moment at the end-effector tip.

### Manipulability Force Analysis
A manipulability force analysis was conducted, constructing ellipsoids at various points along the trajectory. This simulation, orchestrated through command inputs from Matlab, highlighted the versatility of the robotic manipulator in handling dynamic tasks with a constant force application at the end-effector.

## Video Demonstrations

1. [Analytical Analysis](https://iitgnacin-my.sharepoint.com/:v:/g/personal/20110161_iitgn_ac_in/EZT4lCU1CwNCgOqY3S821_0BdEbVt_KYmOSsnq725n6yHQ?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=uN0liq) - 3D infinity trajectory tracking by the end-effector of kuka-iiwa-14-R820 with different analysis.
2. [Motion Control and Torque Optimization](https://iitgnacin-my.sharepoint.com/:f:/g/personal/20110161_iitgn_ac_in/Ek1POail0GdIgq5dxkJgwJ4BRTzaSccDG7H5HPmveywIAQ?e=ZsJFt9) - Demonstration of the robot's performance in motion control and torque optimization during glass transportation.

You can explore the repository for detailed code implementations, simulations, and documentation related to each project phase. Your feedback and contributions are highly appreciated!
