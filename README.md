# Adaptive-rewiring-model

The codes presented here are the main core for simulating the adaptive interaction rewiring model between species of different guilds. Except for “Draw_Ternary.py”, which is written for Python, the rest of the code is written for Matlab. “Model_Core_Line.m” is the core program used to generate data and figures for line drawings in the main text. “Model_Core_Ternary” is the main program used to generate data for ternary diagrams (note that the corresponding network metric data needs to be saved as a txt file, such as “resilience.txt”). “Draw_Ternary.py” is a Python program used to read txt file data and draw it on ternary diagrams (note that “Draw_Ternary.py” and “resilience.txt” must be placed in the same folder).

The Matlab library "BiMat" is used for analyzing the network structures, which is publicly available via the link: https://bimat.github.io/.
