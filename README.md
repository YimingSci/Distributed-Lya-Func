# Distributed Lyapunov Function
Codes for the distributed construction of a global Lyapunov function using sum-of-squares (SOS) optimization. Our computational method is designed to maximize the size of the region of attraction estimated by the constructed Lyapunov function.

See the reference below for more details.

# Usage

- `VDP\_dis\_lya\_func.m` : This code adopts the method from **ref\[1]** to construct a Lyapunov function for the Van der Pol oscillator.

# Pseudocodes


For further details on the formulation and implementation of SOS programming (e.g., using the function `sosprogram` from MATLAB's SOSTOOLS), we refer the user to the tutorial article in Ref. [2].

![image](https://github.com/user-attachments/assets/4f0b8f9c-8654-47b6-917f-0f9db6f3d704)

![image](https://github.com/user-attachments/assets/4a0e8d92-05a7-4c00-ad44-07872d6b8c59)

![image](https://github.com/user-attachments/assets/a77ba12a-90ee-4c9c-94df-47b0a71f50b7)



# Dependences

- `SOS Toolbox`: This code requires the SOS (Sum of Squares) Toolbox to be pre-installed. Link: https://github.com/oxfordcontrol/SOSTOOLS
- `WattsStrogatz Toolbox`: This optional toolbox is used to generate basic coupled oscillator networks. Link: https://www.mathworks.com/help/matlab/math/build-watts-strogatz-small-world-graph-model.html

All codes were tested and run in MATLAB 2024a. To run the codes, download all files in this repository to a folder of your choice and run one of the `main` scripts of your choice. All codes generate/include the required data to run the simulations and optimization; simulations can take a few minutes on a standard laptop.

# License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".

# References
1.  AN Montanari, AED Barioni, C Duan, AE Motter. Optimal flock formation induced by heterogeneity. (2024)
2.  A Packard, U Topcu, PJ Seiler Jr, G Balas. Help on SOS. *IEEE Control Systems Magazine*, 30(4), 18-23 (2010).


