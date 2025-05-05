# Distributed Lyapunov Function
Codes for the distributed construction of a global Lyapunov function using sum-of-squares (SOS) optimization. Our computational method is designed to maximize the size of the region of attraction estimated by the constructed Lyapunov function.

See the reference below for more details.

# Usage

- `main_CMflock.m` : Compares the real-time optimization of the "centroid model" of flocking (Eq. 1 of Ref. 1) for heterogeneous and homogeneous systems. In this model, the agents are tasked to converge to a pre-specified formation while tracking a virtual target moving across space. The communication network is defined by an all-to-all network, weighted according to the relative distance between agents. Data packages are exchanged periodically, resulting in a piecewise-constant model.

- `main_TDconsensus.m` : Compares the stability optimization of the time-delay consensus model (Eq. 8 of Ref. 1) for heterogeneous and homogeneous systems. In this model, the agents must achieve consensus in position and velocity, converging to the same state. If the system is unstable, the flock of agents fragments into isolated groups.

- `main_OSflock.m` : Compares the real-time optimization of the "Olfati-Saber model" of flocking (Eq. 9 of Ref. 1) for heterogeneous and homogeneous systems. In this model, the agents are tasked to form a cohesive lattice structure while tracking a virtual target. The communication range is limited, resulting in a sparse interaction network. The lattice structure is emergent since the final configuration depends on the agents' initial conditions. The simulation can be performed for agents moving in free space or maneuvering around obstacles.

- `main_distributedopt.m` : Distributed optimization of the largest Lyapunov exponent in the "centroid model" of flocking (Eq. 1 of Ref. 1) for heterogeneous systems. In the distributed optimization, agents only have access to local information of the state of agents within some specified distance radius. The optimization is thus solved in parallel for each agent `i`, returning the corresponding optimal feedback gains. 

# Pseudocodes


For further details on the formulation and implementation of SOS programming (e.g., using the function `sosprogram` from MATLAB's SOSTOOLS), we refer the user to the tutorial article in Ref. [2].

![image](https://github.com/user-attachments/assets/4f0b8f9c-8654-47b6-917f-0f9db6f3d704)

![image](https://github.com/user-attachments/assets/4a0e8d92-05a7-4c00-ad44-07872d6b8c59)

![image](https://github.com/user-attachments/assets/a77ba12a-90ee-4c9c-94df-47b0a71f50b7)



# Dependences

- `EigOptimization` : This folder contains optimization routines (`beta_opt.m`) used to minimize the Lyapunov exponent (calculated in `opteigreal.m`) associated with the corresponding flocking models (`CM`, `TD`, `OS`) for heterogeneous (`het`) and homogeneous (`hom`) systems.
  
- `FlockODEs` : This folder contains the ODEs describing the CM (`CMflock_piecewise.m`) and OS (`OSflock.m`) flocking models. The functions `agent_coord.m` compute the performance metrics (tracking error, lattice deviation energy, relative connectivty) and the XY coordinates of the agents/target for each time instant in the corresponding models.

- The  time-delay consensus model (`main_TDconsensus.m`) requires installation of the DDE-BIFTOOL toolbox. These codes were tested on version 3.0. See Ref. 3 and their website (https://sourceforge.net/projects/ddebiftool/) for more details.

All codes were tested and run in MATLAB 2023a. To run the codes, download all files in this repository to a folder of your choice and run one of the `main` scripts of your choice. All codes generate/include the required data to run the simulations and optimization; simulations can take a few minutes on a standard laptop.

# License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".

# References
1.  AN Montanari, AED Barioni, C Duan, AE Motter. Optimal flock formation induced by heterogeneity. (2024)
2.  A Packard, U Topcu, PJ Seiler Jr, G Balas. Help on SOS. *IEEE Control Systems Magazine*, 30(4), 18-23 (2010).


