# WEC-Sim-Python
**WEC-Sim-Python** is the Python version of the [WEC-Sim (Wave Energy Converter SIMulator)](https://github.com/WEC-Sim/WEC-Sim.git) which was originally developed in MATLAB/SIMULINK. 

## Goal of WEC-Sim-Python
**WEC-Sim-Python** aims to help researchers, start up companies, and enthusiasts without access to MATLAB to use the open source code provided from NREL and Sandia lab. Also with growing research in the machine learning, **WEC-Sim-Python** could provid easier transition for those who learnd machine learning in Python.

## Build Status
**WARNING: WEC-Sim-Python has not been completed.** Estimate date of completion: **February, 2021**
### Currently supports
- [x] Wave class
- [x] Body class: Testing nlHydro
### Future update
- [ ] Other source class and BEMIO 
- [ ] Open AI Gym for simulation

## Major changes from original WEC-Sim
1. Does not support Simulink.
3. A seperate Paraview class is made in the objects directory.
4. Modified some methods and file structures for faster computation but the original method has been tested and left in the comment.

## Precedure
To be updated once upon completion.

### Note
1. Due to the  differences in significant digits or methods of calculatioins, all variables have been tested up to 10^-7 decimal places.
