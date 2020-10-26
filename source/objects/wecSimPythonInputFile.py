
import simulationClass
import waveClass
import bodyClass
import constraintClass
import ptoClass
## Simulation Data
simu = simulationClass.SimulationClass()               # Initialize simulationClass
simu.startTime = 0                     # Simulation Start Time [s]
simu.rampTime = 100                   	# Wave Ramp Time [s]
simu.endTime=400                       # Simulation End Time [s]
simu.dt = 0.1 							# Simulation time-step [s]

## Wave Information 
# Regular Waves  
waves = waveClass.WaveClass('type')              # Initialize waveClass
waves.T = 999                          # Wave Period [s]
waves.H = 999                          # Wave Height [m]

## Body Data
# Float
body_1 = bodyClass.BodyClass('rm3.h5')          # Initialize bodyClass for Float      
body_1.geometryFile = 'float.stl'    # Geomtry File
body_1.mass = 'equilibrium'                   # Mass [kg]
body_1.momOfInertia = [20907301, 21306090.66, 37085481.11]  # Moment of Inertia [kg*m^2]     

# Spar/Plate
body_2 = bodyClass.BodyClass('rm3.h5')          # Initialize bodyClass for Spar/Plate
body_2.geometryFile = 'plate.stl'    # Geometry File   
body_2.mass = 'equilibrium'                   # Mass [kg]
body_2.momOfInertia = [94419614.57, 94407091.24, 28542224.82]   # Moment of Inertia [kg*m^2]     

## PTO and Constraint Parameters
# Floating (3DOF) Joint
constraint_1 = constraintClass.ConstraintClass('Constraint1') # Initialize constraintClass for Constraint1
constraint_1.loc = [0, 0, 0]                    # Constraint Location [m]

# Translational PTO
pto_1 = ptoClass.PtoClass('PTO1')                      # Initialize ptoClass for PTO1
pto_1.k = 999                                 # PTO Stiffness [N/m]
pto_1.c = 999                                 # PTO Damping [N/(m/s)]
pto_1.loc = [0, 0, 0]                           # PTO Location [m]