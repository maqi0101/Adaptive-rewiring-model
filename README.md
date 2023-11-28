# Adaptive Rewiring Model

The code presented in this repository constitutes the core for simulating the paper titled "Adaptive rewiring shapes structure and stability of a three-guild herbivore-plant-pollinator network."

## Code Description

- **Draw_Ternary.py**: Python code for drawing ternary diagrams.
- **LV_PHM.m**: Matlab function code describing Community Dynamics (differential equations Eqn. 4a-4c).
- **cal_structure.m**: Matlab code for analyzing network Modularity and Nestedness.
- **get_jacmat.m**: Matlab code for obtaining the Jacobian matrix.
- **overlap.m**: Matlab code for calculating the niche overlap between two species.
- **rewiring.m**: Matlab code describing the adaptive interaction rewiring process.
- **wire_back.m**: Matlab code for restoring the original interaction relationship if adaptive rewiring is not successful.

## Core Programs

- **Model_Core_Line.m**: Matlab core program for generating data and figures for line figures in the main text (e.g., Fig.2c-2e; Fig.3b,3e).
- **Model_Core_Ternary.m**: Matlab main program for generating data for ternary diagrams (requires saving corresponding network metric data as a txt file, e.g., “resilience.txt”; e.g., Fig.2a,2b; Fig.3c,3f; Fig.4).
- **Draw_Ternary.py**: Python program used to read txt file data and draw it on ternary diagrams (note that "Draw_Ternary.py" and "resilience.txt" must be placed in the same folder).

## External Library

The Matlab library "BiMat" is utilized for analyzing network structures. It is publicly available [here](https://bimat.github.io/).
