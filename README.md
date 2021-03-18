A hardware-in-the-loop simulator for physical human-aerial manipulator cooperation
===============
This software allows a simulated unmanned aerial manipulator (UAM) to interact with the real world through an hardware interface. The UAM can perform both autonomous and human/UAM collaborative tasks.

The simulator provides the following elements:
 - A UAM model already implemented inside the Gazebo simulator. The UAV odometry is obtained through the [RotorS](https://github.com/ethz-asl/rotors_simulator) odometry plugin;
 - A bilateral communication interface connecting the hardware with the simulated model;
 - Some example controllers for the UAM and the hardware interface to test the proposed architecture.

The UAM model consists of a quadrotor and a 6DoFs manipulator attached under it.

The hardware interface can be any device as long as its position in the space can be commanded (e.g., manipulators and/or haptic interfaces). The hardware also needs a force sensor, to provide the real force feedback at the the UAM.

Installation Instructions - Ubuntu 18.04 with ROS Melodic
---------------------------------------------------------
 1. Install and initialize ROS Melodic desktop full:
 ```
 $ sudo sh -c 'echo "deb http://packages.ros.org/ros/ubuntu $(lsb_release -sc) main" > /etc/apt/sources.list.d/ros-latest.list'
 $ sudo apt-key adv --keyserver 'hkp://keyserver.ubuntu.com:80' --recv-key C1CF6E31E6BADE8868B172B4F42ED6FBAB17C654
 $ sudo apt-get update
 $ sudo apt-get install ros-melodic-desktop-full
 $ sudo rosdep init
 $ rosdep update
 $ source /opt/ros/melodic/setup.bash
 ```
 2. If you don't have ROS a workspace yet, you can create one by
 ```
 $ mkdir -p ~/catkin_ws/src
 $ cd ~/catkin_ws/src
 $ catkin_init_workspace  # initialize your catkin workspace
 ```
 3. Get and install the [RotorS](https://github.com/ethz-asl/rotors_simulator) simulator.
 4. Clone this repository into your workspace and compile it with
 ```
 $ cd ~/catkin_ws/src
 $ git clone https://github.com/prisma-lab/HIL_airmanip.git
 $ catkin_make
```

Folders
---------------------------------------------------------
**aerialmanip_control**

The *aerialmanip_control* folder contains the Gazebo simulation world file, the quadrotor xacro description and the UAM controllers.

**arm_description**

The *arm_description* folder contains the UAM arm xacro description and its joint position controllers.

**connection_plugin**

The *connection_plugin* folder contains the Gazebo plugin used to connect the simulated and the real world. This plugin reads the hardware measured force and applies it to the UAM.

**hardware_controller**

The *hardware_controller* folder contains the hardware control system and can be modified accordingly to the user chosen hardware interface. The provided hardware controller considers the KUKA LBR IIWA manipulator.
