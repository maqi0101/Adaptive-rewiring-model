# Adaptive-rewiring-model

The codes presented here are the main core for simulating the paper "Adaptive interaction rewiring alters the structure and stability of a three-guild herbivore-plant-pollinator network". 

The code "Draw_ Ternary.py" was written using Python, and the other code files were written using Matlab. “LV_PHM.m” is the function code used to describe Community Dynamics (differential equations Eqn. 4a-4c). “cal_structure.m” is used to analyze network Modularity and Nestedness. “get_jacmat.m” is used to get the Jacobian matrix . “overlap.m” is used to calculate the niche overlap between two species. “rewiring.m” is used to describe the adaptive interaction rewiring process. “wire_back.m” is used to restore the original interaction relationship if adaptive rewiring is not successful. 

"Model_Core_Line.m" is the core program used to generate data and figures for line figures in the main text (e.g., Fig.2c-2e; Fig.3b,3e). "Model_Core_Ternary" is the main program used to generate data for ternary diagrams (note that the corresponding network metric data needs to be saved as a txt file, such as “resilience.txt”; e.g., Fig.2a,2b; Fig.3c,3f; Fig.4). "Draw_ Ternary.py" is a Python program used to read txt file data and draw it on ternary diagrams (note that "Draw_ Ternary.py" and "resilience.txt"  must be placed in the same folder).

The Matlab library "BiMat" is used for analyzing the network structures, which is publicly available via the link: https://bimat.github.io/.
