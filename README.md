# WEC-Sim-Python
**WEC-Sim-Python** is the Python version of the [WEC-Sim (Wave Energy Converter SIMulator)](https://github.com/WEC-Sim/WEC-Sim.git), which was originally developed in MATLAB/SIMULINK. 

## Goal of WEC-Sim-Python
**WEC-Sim-Python** aims to help researchers, start-up companies, and enthusiasts without access to MATLAB in order to use the open-source code provided by NREL and Sandia lab. Also, with growing research in the field of machine learning, **WEC-Sim-Python** could be more convenient for those who develop machine learning projects utilizing Python.

## Build Status
**Note: WEC-Sim-Python has not been completed.** Estimated date of completion: **February, 2021**
### Currently supports
- [x] Simulation class: getWecSimPythonVer will be tested upon pre-release
- [x] Wave class
- [x] Body class
### Future update
- [ ] Other source class and BEMIO 
- [ ] Open AI Gym for simulation

## Major changes from original WEC-Sim
1. Do not support Simulink.
2. A seperate Paraview class is incorporated in the objects directory.
3. Modified some methods and file structures for faster computation, and the original method has been tested and incorporated in the comment.

## Procedure
To be updated upon the completion.

### Note
1. Due to the difference in the significant digits or methods of the calculations, all variables have been tested up to 10^-7 decimal places.
