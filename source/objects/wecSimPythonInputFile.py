import os

import simulationClass
import waveClass
import bodyClass
import constraintClass
import ptoClass

# class wecSimPythonInputFileClass:
#     # This class contains WEC-Sim simulation parameters and settings
def setUpSimu():
    # Simulation Data
    simu = simulationClass.SimulationClass()               # Initialize simulationClass
    simu.startTime = 0                     # Simulation Start Time [s]
    simu.rampTime = 100                   	# Wave Ramp Time [s]
    simu.endTime=400                       # Simulation End Time [s]
    simu.dt = 0.1 							# Simulation time-step [s]
    return(simu)

def setUpWaves():    
    ## Wave Information 
    # Regular Waves  
    waves = waveClass.WaveClass('regular')              # Initialize waveClass
    waves.T = 2.5                          # Wave Period [s]
    waves.H = 8                          # Wave Height [m]
    return(waves)

def setUpBody():
    ## Body Data
    # Float 
    body_0 = bodyClass.BodyClass('rm3.h5')          # Initialize bodyClass for Float    
    body_0.geometryFile = 'float.stl'    # Geomtry File
    body_0.mass = 'equilibrium'                   # Mass [kg]
    body_0.momOfInertia = [20907301, 21306090.66, 37085481.11]  # Moment of Inertia [kg*m^2]     
    
    # Spar/Plate
    body_1 = bodyClass.BodyClass('rm3.h5')          # Initialize bodyClass for Spar/Plate
    body_1.geometryFile = 'plate.stl'    # Geometry File   
    body_1.mass = 'equilibrium'                   # Mass [kg]
    body_1.momOfInertia = [94419614.57, 94407091.24, 28542224.82]   # Moment of Inertia [kg*m^2]     
    return(body_0,body_1)

def setUpConstraint():
    ## PTO and Constraint Parameters
    # Floating (3DOF) Joint
    constraint_1 = constraintClass.ConstraintClass('Constraint1') # Initialize constraintClass for Constraint1
    constraint_1.loc = [0, 0, 0]                    # Constraint Location [m]
    return(constraint_1)

def setUpPto():
    # Translational PTO
    pto_1 = ptoClass.PtoClass('PTO1')                      # Initialize ptoClass for PTO1
    pto_1.k = 0                                 # PTO Stiffness [N/m]
    pto_1.c = 1200000                                 # PTO Damping [N/(m/s)]
    pto_1.loc = [0, 0, 0]                           # PTO Location [m]
    return(pto_1)