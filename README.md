# WEC-Sim-Python
WEC-Sim-Python is the python version of the WEC-Sim provided from DOE. 
WEC-Sim-Python modified some methods and file structures for faster computation.
Original method has been tested and left in a comment. 

Major changes:
1. WEC-Sim-Python does not uses Simulink.
2. Can use Unity 3D Engin, Open AI Gym, and Paraview to make simulations.
3. A seperate Paraview class is made in the objects directory.
4. Can import both mat file and txt file which allows users without MATLAB to have full compatibility.

# Procedure
1. Choose a BEM solver(NEMOH,etc) and run BEMIO.py which is placed in source/objects/functions directory.
2. Run wecsim-python.py in source directory. 
3. If you want to check out procedure to run each class refer to {link} <- I'll put link later

Note: 
1. Due to the  differences in significant digits or methods of calculatioins, all variables have been tested up to 10^-7 decimal places.

References:
View WEC-Sim from: https://github.com/WEC-Sim/WEC-Sim.git
