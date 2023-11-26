# Adaptive-rewiring-model

The codes presented here are the main core for simulating the paper "Adaptive interaction rewiring alters the structure and stability of a three-guild herbivore-plant-pollinator network". The code "Draw_Ternary.py" was written using Python, and the other code files were written using Matlab. "Model_Core_Line.m" is the core program used to generate data and figures for line drawings in the main text. "Model_Core_Ternary" is the main program used to generate data for ternary diagrams (note that the corresponding network metric data needs to be saved as a txt file, such as “resilience.txt”). "Draw_Ternary.py" is a Python program used to read txt file data and draw it on ternary diagrams (note that "Draw_Ternary.py" and "resilience.txt"  must be placed in the same folder).

The Matlab library "BiMat" is used for analyzing the network structures, which is publicly available via the link: https://bimat.github.io/.
