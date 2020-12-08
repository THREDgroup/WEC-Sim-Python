# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 07:52:07 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

Note: All the struct() has been changed to dictionary

Note: trimesh is used to obtain bodyGeometry properties

"""
import h5py
import numpy as np
import numpy.matlib 
import warnings
import trimesh

from numpy.linalg import inv
from scipy import interpolate
from copy import copy

class BodyClass:
    def hdf5FileProperties(self):#hdf5 file properties
        """
        Hydrodynamic data from BEM or user defined.
        
        """
        self.hydroData={'simulation_parameters':{'scaled':[],
                                                 'wave_dir':[],
                                                 'water_depth':[],
                                                 'w':[],
                                                 'T':[]},
                        'properties':{'name':[],
                                      'body_number':[],
                                      'cg':[],
                                      'cb':[],
                                      'disp_vol':[],
                                      'dof_start':[],
                                      'dof_end':[]},
                        'hydro_coeffs':{'linear_restoring_stiffness':[],
                                        'excitation':{'re':[],
                                                      'im':[],
                                                      'impulse_response_fun':{'f':[],
                                                                              't':[]}},
                                        'added_mass':{'all':[],
                                                      'inf_freq':[]},
                                        'radiation_damping':{'all':[],
                                                             'impulse_response_fun':{'K':[],
                                                                                     't':[]},
                                                             'state_space':{'it':[],
                                                                            'A':{'all':[]},
                                                                            'B':{'all':[]},
                                                                            'C':{'all':[]},
                                                                            'D':{'all':[]}}},
                                        'mean_drift':[]},
                        'gbm':{'mass':[],
                               'stiffness':[],
                               'damping':[]}}                                                  

    def inputFileProperties(self):#input file properties
        self.name              = []                                            # Body name. For WEC bodies this is given in the h5 file.
        self.mass              = []                                            # Mass in kg or specify 'equilibrium' to have mass= dis vol * density
        self.momOfInertia      = []                                            # Moment of inertia [Ixx, Iyy, Izz] in kg*m^2
        self.cg                = []                                            # Center of gravity [x, y, z] in meters. For WEC bodies this is given in the h5 file.
        self.cb                = []                                            # Center of buoyancy [x, y, z] in meters. For WEC bodies this is given in the h5 file.
        self.dispVol           = []                                            # Displaced volume at equilibrium position in meters cubed. For WEC bodies this is given in the h5 file.
        self.dof               = []                                            # Number of DOFs. For WEC bodies this is given in the h5 file. IF not, default is 6
        self.dof_gbm           = []                                            # Number of DOFs for GBM.
        self.dof_start         = []                                            # Index of DOF starts. For WEC bodies this is given in the h5 file. IF not, default is (bodyNumber-1)*6
        self.dof_end           = []                                            # Index of DOF ends. For WEC bodies this is given in the h5 file. IF not, default is (bodyNumber-1)*6+5
        self.geometryFile      = 'NONE'                                        # Location of geomtry stl files
        self.viscDrag          = {                                             # Structure defining the viscous (quadratic) drag
                                  'Drag':                 [0,0,0,0,0,0],       # Viscous (quadratic) drag, matrix 6x6
                                  'cd':                   [0,0,0,0,0,0],       # Viscous (quadratic) drag cd coefficient, vector length 6
                                  'characteristicArea':   [0,0,0,0,0,0]}       # Characteristic area for viscous drag, vector length 6
        self.initDisp          = {                                             # Structure defining the initial displacement
                                  'initLinDisp':          [0,0,0],             # Initial displacement of center fo gravity - used for decay tests (format: [displacment in m], default = [0, 0, 0])
                                  'initAngularDispAxis':  [0,1,0],             # Initial displacement of cog - axis of rotation - used for decay tests (format: [x y z], default = [1, 0, 0])
                                  'initAngularDispAngle': 0}                   # Initial displacement of cog - Angle of rotation - used for decay tests (format: [radians], default = 0)
        self.hydroStiffness    = [0,0,0,0,0,0]                                 # Hydrostatic stiffness matrix overrides BEMIO definition, matrix 6x6
        self.linearDamping     = [0,0,0,0,0,0]                                 # Linear damping coefficient, matrix size of 6x6
        self.userDefinedExcIRF = []                                            # Excitation IRF from BEMIO used for User-Defined Time-Series
        self.viz               = {                                             # Structure defining visualization properties
                                  'color': [1,1,0],                            # Visualization color for either SimMechanics Explorer or Paraview.
                                  'opacity': 1}                                # Visualization opacity for either SimMechanics Explorer or Paraview.
        self.morisonElement   = {                                              # Structure defining the Morrison Elements
                                 'cd':                 [0,0,0],                # Viscous (quadratic) drag cd, vector length 3
                                 'ca':                 [0,0,0],                # Added mass coefficent for Morrison Element (format [Ca_x, Ca_y, Ca_z], default = [0 0 0])
                                 'characteristicArea': [0,0,0],                # Characteristic area for Morrison Elements calculations (format [Area_x, Area_y, Area_z], default = [0 0 0])
                                 'VME':                 0     ,                # Characteristic volume for Morrison Element (default = 0)
                                 'rgME':               [0,0,0]}                # Vector from center of gravity to point of application for Morrison Element (format [X, Y, Z], default = [0 0 0]).
        self.nhBody            = 0                                             # Flag for non-hydro body.
        self.flexHydroBody     = 0                                             # Flag for flexible body. 
        self.meanDriftForce    = 0                                             # Flag for mean drift force. 0: No 1: from control surface 2: from momentum conservation.
        
    def bodyGeometryFileProperties(self):                                      # body geometry stl file properties
        self.bodyGeometry = {                                                  # Structure defining body's mesh
                             'numFace': [],                                    # Number of faces
                             'numVertex': [],                                  # Number of vertices
                             'vertex': [],                                     # List of vertices
                             'face': [],                                       # List of faces
                             'norm': [],                                       # List of normal vectors
                             'area': [],                                       # List of cell areas
                             'center': []}                                     # List of cell centers
        self.meshFile    = []                                                  # body geometry stl file from Trimesh api
    
    def internalProperties(self,filename):                                     # internal properties 
        self.hydroForce        = {'linearHydroRestCoef':[],                    # Hydrodynamic forces and coefficients used during simulation.
                                  'visDrag':[],
                                  'linearDamping':[],
                                  'userDefinedFe':[],
                                  'fExt':{'re':[],
                                          'im':[],
                                          'md':[]},
                                  'fAddedMass':[],
                                  'fDamping':[],
                                  'irkb':[],
                                  'ssRadf':{'A':[],
                                            'B':[],
                                            'C':[],
                                            'D':[]},
                                  'storage':{'mass':[],
                                             'momOfInertia':[],
                                             'fAddedMass':[],
                                             'output_forceAddedMass':[],
                                             'output_forceTotal':[]}}           
        self.tmp               ={'fadm':[],                                    # temporary file
                                 'adjmass':[],
                                 'mass':[],
                                 'momOfInertia':[],
                                 'hydroForce_fAddedMass':[]}                    
        self.h5File            = filename                                      # hdf5 file containing the hydrodynamic data
        self.hydroDataBodyNum  = []                                            # Body number within the hdf5 file.
        self.massCalcMethod    = []                                            # Method used to obtain mass: 'user', 'fixed', 'equilibrium'
        self.bodyNumber        = []                                            # bodyNumber in WEC-Sim as defined in the input file. Can be different from the BEM body number.
        self.bodyTotal         = 0                                             # Total number of WEC-Sim bodies (body block iterations)
        self.lenJ              = []                                            # Matrices length. 6 for no body-to-body interactions. 6*numBodies if body-to-body interactions.

    def __init__(self,filename):
        """
        Initialize Body Class
        Takes string parameter called filename
        filename: name of h5 file

        """
        self.bodyGeometryFileProperties()
        self.internalProperties(filename)
        self.inputFileProperties()
        self.hdf5FileProperties()
        self.meanDriftForce = 0

    def readH5file(self):
        """
        Read and recond properties of h5 file in self.hydroData

        """
        f = h5py.File(self.h5File, 'r')
        name = '/body' + str(self.bodyNumber)
        self.cg = np.transpose(np.array(f.get(name + '/properties/cg')))
        self.cb = np.transpose(np.array(f.get(name + '/properties/cb')))
        self.dispVol = np.array(f.get(name + '/properties/disp_vol'))
        self.name = np.string_(np.array(f.get(name + '/properties/name'))).decode("utf-8")
        self.hydroData['simulation_parameters']['scaled'] = np.array(f.get('/simulation_parameters/scaled'))
        self.hydroData['simulation_parameters']['wave_dir'] = np.transpose(np.array(f.get('/simulation_parameters/wave_dir')))
        if np.array(f.get('/simulation_parameters/water_depth')).dtype == float:
            self.hydroData['simulation_parameters']['water_depth'] = np.array(f.get('/simulation_parameters/water_depth'))
        else:            
            self.hydroData['simulation_parameters']['water_depth'] = np.string_(np.array(f.get('/simulation_parameters/water_depth'))).decode("utf-8")
        self.hydroData['simulation_parameters']['w'] = np.transpose(np.array(f.get('/simulation_parameters/w')))
        self.hydroData['simulation_parameters']['T'] = np.transpose(np.array(f.get('/simulation_parameters/T')))
        self.hydroData['properties']['name'] = np.string_(np.array(f.get(name + '/properties/name'))).decode("utf-8")
        self.hydroData['properties']['body_number'] = np.array(f.get(name + '/properties/body_number'))
        self.hydroData['properties']['cg'] = np.transpose(np.array(f.get(name + '/properties/cg')))
        self.hydroData['properties']['cb'] = np.transpose(np.array(f.get(name + '/properties/cb')))
        self.hydroData['properties']['disp_vol'] = np.array(f.get(name + '/properties/disp_vol'))
        if np.array(f.get(name +'/properties/dof')).all() != None:
            self.hydroData['properties']['dof'] = np.array(f.get(name +'/properties/dof')) 
        else:
            self.hydroData['properties']['dof'] = np.array(6)
        if np.array(f.get(name + '/properties/dof_start')).all() != None:
            self.hydroData['properties']['dof_start'] = np.array(f.get(name + '/properties/dof_start'))
        else:
            self.hydroData['properties']['dof_start'] = np.array((self.bodyNumber-1)*6+1)
        if np.array(f.get(name + '/properties/dof_end')).all() != None:
            self.hydroData['properties']['dof_end'] = np.array(f.get(name + '/properties/dof_end'))
        else:
            self.hydroData['properties']['dof_end'] = np.array((self.bodyNumber-1)*6+6)
        self.dof       = self.hydroData['properties']['dof']
        self.dof_start = self.hydroData['properties']['dof_start']
        self.dof_end   = self.hydroData['properties']['dof_end']
        self.dof_gbm   = self.dof-6
        self.hydroData['hydro_coeffs']['linear_restoring_stiffness'] = np.transpose(np.array(f.get(name + '/hydro_coeffs/linear_restoring_stiffness')))
        self.hydroData['hydro_coeffs']['excitation']['re'] = np.array(f.get(name +  '/hydro_coeffs/excitation/re'))
        self.hydroData['hydro_coeffs']['excitation']['im'] = np.array(f.get(name + '/hydro_coeffs/excitation/im'))
        if np.array(f.get(name + '/hydro_coeffs/excitation/impulse_response_fun/f')).all() != None:
            self.hydroData['hydro_coeffs']['excitation']['impulse_response_fun']['f'] = np.array(f.get(name + '/hydro_coeffs/excitation/impulse_response_fun/f'))
        if np.array(f.get(name + '/hydro_coeffs/excitation/impulse_response_fun/t')).all() != None:
            self.hydroData['hydro_coeffs']['excitation']['impulse_response_fun']['t'] = np.array(f.get(name + '/hydro_coeffs/excitation/impulse_response_fun/t'))
        self.hydroData['hydro_coeffs']['added_mass']['all'] = np.array(f.get(name + '/hydro_coeffs/added_mass/all'))
        self.hydroData['hydro_coeffs']['added_mass']['inf_freq'] = np.array(f.get(name + '/hydro_coeffs/added_mass/inf_freq'))
        self.hydroData['hydro_coeffs']['radiation_damping']['all'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/all'))
        if np.array(f.get(name + '/hydro_coeffs/radiation_damping/impulse_response_fun/K')).all() != None:
            self.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['K'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/impulse_response_fun/K'))
        if np.array(f.get(name +'/hydro_coeffs/radiation_damping/impulse_response_fun/t')).all() != None:
            self.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['t'] = np.transpose(np.array(f.get(name +'/hydro_coeffs/radiation_damping/impulse_response_fun/t')))
        if np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/it')).all() != None:
            self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/it'))
        if np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/A/all')).all() != None:
            self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/A/all'))
        if np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/B/all')).all() != None:
            self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/B/all'))
        if np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/C/all')).all() != None:
            self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/C/all'))
        if np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/D/all')).all() != None:
            self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['D']['all'] = np.array(f.get(name +'/hydro_coeffs/radiation_damping/state_space/D/all'))
        if np.array(f.get(name + '/properties/mass')).all() != None:
            tmp = np.array(f.get(name + '/properties/mass'))
            self.hydroData['gbm']['mass'] = [tmp[0].arange(self.dof_start+5,self.dof_end+1),tmp[1].arange(self.dof_start+5,self.dof_end+1)]
        if np.array(f.get(name + '/properties/stiffness')).all() != None:
            tmp = np.array(f.get(name + '/properties/stiffness'))
            self.hydroData['gbm']['stiffness'] = [tmp[0].arange(self.dof_start+5,self.dof_end+1),tmp[1].arange(self.dof_start+5,self.dof_end+1)]
        if np.array(f.get(name + '/properties/damping')).all() != None:
            tmp = np.array(f.get(name + '/properties/damping'))
            self.hydroData['gbm']['damping'] = [tmp[0].arange(self.dof_start+5,self.dof_end+1),tmp[1].arange(self.dof_start+5,self.dof_end+1)]
        if self.meanDriftForce == 0:
            self.hydroData['hydro_coeffs']['mean_drift'] = 0.*self.hydroData['hydro_coeffs']['excitation']['re']
        elif self.meanDriftForce == 1:
            self.hydroData['hydro_coeffs']['mean_drift'] = np.array(f.get(name + '/hydro_coeffs/mean_drift/control_surface/val'))
        elif self.meanDriftForce == 2:
            self.hydroData['hydro_coeffs']['mean_drift'] = np.array(f.get(name + '/hydro_coeffs/mean_drift/momentum_conservation/val'))
        else:
            warnings.warn("Wrong flag for mean drift force.",DeprecationWarning)
        f.close()
        
    def loadHydroData(self, hydroData):
        """
        Loads user defiend hydroData structure as alternative
        to reading the h5 file. Used in wecSimMCR

        """
        self.hydroData = hydroData
        self.cg        = hydroData['properties']['cg']
        self.cb        = hydroData['properties']['cb']
        self.dispVol   = hydroData['properties']['disp_vol']
        self.name      = hydroData['properties']['name']
        self.dof       = self.hydroData['properties']['dof']
        self.dof_start = self.hydroData['properties']['dof_start']
        self.dof_end   = self.hydroData['properties']['dof_end']
        self.dof_gbm   = self.dof-6
    
    def hydroForcePre(self,w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B):
        """
        This is the most important method of this class
        HydroForce Pre-processing calculations
        1. Set the linear hydrodynamic restoring coefficient, viscous
            drag, and linear damping matrices
        2. Set the wave excitation force    

        """
        self.setMassMatrix(rho,nlHydro)
        if self.dof_gbm > 0:
            #self.linearDamping = [self.linearDamping np.zeros(1,self.dof-np.size(self.linearDamping))] # to use this you need to define self.linearDamping
            tmp0 = self.linearDamping
            tmp1 = np.size(self.linearDamping)
            self.linearDamping = np.zeros(self.dof[0])                
            self.linearDamping[0][:tmp1[0]] = tmp0[0]
            self.linearDamping[1][:tmp1[1]] = tmp0[1]

            tmp0 = self.viscDrag['Drag']
            tmp1 = np.size(self.viscDrag['Drag'])
            self.viscDrag['Drag'] = np.zeros(self.dof)                
            self.viscDrag['Drag'][0][:tmp1[0]] = tmp0[0]
            self.viscDrag['Drag'][1][:tmp1[1]] = tmp0[1]
            
            self.viscDrag['cd']   = np.append(self.viscDrag['cd'], np.zeros(self.dof[0]-np.size(self.viscDrag['cd'])))
            self.viscDrag['characteristicArea'] = np.append(self.viscDrag['characteristicArea'],np.zeros(1,self.dof-np.size(self.viscDrag['characteristicArea'])))

        if self.hydroStiffness.any() == 1:  #check if self.hydroStiffness is defined
            self.hydroForce['linearHydroRestCoef'] = self.hydroStiffness
        else:
            k = self.hydroData['hydro_coeffs']['linear_restoring_stiffness']#(:,self.dof_start:self.dof_end)
            self.hydroForce['linearHydroRestCoef'] = k*rho*g

        if  self.viscDrag['Drag'].any() == 1:  #check if self.viscDrag['Drag'] is defined
            self.hydroForce['visDrag'] = self.viscDrag['Drag']
        else:
            self.hydroForce['visDrag'] = np.diag(0.5*rho*self.viscDrag['cd']*self.viscDrag['characteristicArea'])

        self.hydroForce['linearDamping'] = self.linearDamping
        self.hydroForce['userDefinedFe'] = np.zeros((len(waveAmpTime[1]),int(self.dof[0])))  #initializing userDefinedFe for non imported wave cases
        if waveType == 'noWave':
            self.noExcitation()
            self.constAddedMassAndDamping(w,CIkt,rho,B2B)
        elif waveType == 'noWaveCIC':
            self.noExcitation()
            self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        elif waveType == 'regular':
            self.regExcitation(w,waveDir,rho,g)
            self.constAddedMassAndDamping(w,CIkt,rho,B2B)
        elif waveType == 'regularCIC':
            self.regExcitation(w,waveDir,rho,g)
            self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        elif waveType == 'irregular' or waveType == 'spectrumImport':
            self.irrExcitation(w,numFreq,waveDir,rho,g)
            self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        elif waveType == 'etaImport':
            self.userDefinedExcitation(waveAmpTime,dt,waveDir,rho,g)
            self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)

        gbmDOF = self.dof_gbm
        if gbmDOF>0:
            self.hydroForce['gbm']['stiffness']=self.hydroData['gbm']['stiffness']
            self.hydroForce['gbm']['damping']=self.hydroData['gbm']['damping']
            self.hydroForce['gbm']['mass_ff']=[self.hydroForce['fAddedMass'].arange(7,self.dof+1)[self.hydroForce['fAddedMass'].arange(self.dof_start+6,self.dof_end+1)]]+self.hydroData['gbm']['mass']   # need scaling for hydro part
            self.hydroForce['fAddedMass'][7:self.dof+1] = np.zeros(len(np.arange(7,self.dof+1)))
            self.hydroForce['fAddedMass'][(self.dof_start[0]+6):(self.dof_end[0]+1)] = np.zeros(len(np.arange(self.dof_start+6,self.dof_end+1)))
            self.hydroForce['gbm']['mass_ff_inv']=inv(self.hydroForce['gbm']['mass_ff'])
            
            # state-space formulation for solving the GBM
            self.hydroForce['gbm']['state_space']['A'] = [np.zeros((gbmDOF,gbmDOF)), np.eye(gbmDOF,gbmDOF)-inv(self.hydroForce['gbm']['mass_ff'])*self.hydroForce['gbm']['stiffness'],-inv(self.hydroForce['gbm']['mass_ff'])*self.hydroForce['gbm']['damping']]    # move to ... hydroForce sector with scaling .         # or create a new fun for all flex parameters
            self.hydroForce['gbm']['state_space']['B'] = np.eye(2*gbmDOF,2*gbmDOF)
            self.hydroForce['gbm']['state_space']['C'] = np.eye(2*gbmDOF,2*gbmDOF)
            self.hydroForce['gbm']['state_space']['D'] = np.zeros((2*gbmDOF,2*gbmDOF))
            self.flexHydroBody = 1
            self.nhBody=0
            
    def adjustMassMatrix(self,adjMassWeightFun,B2B):
        """
        Merge diagonal term of added mass matrix to the mass matrix
        1. Store the original mass and added-mass properties
        2. Add diagonal added-mass inertia to moment of inertia
        3. Add the maximum diagonal traslational added-mass to body
        mass - this is not the correct description

        """
        iBod = self.bodyNumber
        self.hydroForce['storage']['mass'] = copy(self.mass) # use copy method to set it as seperate variable
        self.hydroForce['storage']['momOfInertia'] = copy(self.momOfInertia)
        self.hydroForce['storage']['fAddedMass'] = copy(self.hydroForce['fAddedMass'])
        if B2B == 1:
            self.tmp['fadm'] = np.diag(self.hydroForce['fAddedMass'],k =(iBod-1)*6)
            self.tmp['adjmass'] = sum(self.tmp['fadm'][0:3]*adjMassWeightFun)
            self.mass = self.mass + self.tmp['adjmass']
            self.momOfInertia = self.momOfInertia+self.tmp['fadm'][3:6]
            self.hydroForce['fAddedMass'][0,(iBod-1)*6] = self.hydroForce['fAddedMass'][0,(iBod-1)*6] - self.tmp['adjmass']
            self.hydroForce['fAddedMass'][1,1+(iBod-1)*6] = self.hydroForce['fAddedMass'][1,1+(iBod-1)*6] - self.tmp['adjmass']
            self.hydroForce['fAddedMass'][2,2+(iBod-1)*6] = self.hydroForce['fAddedMass'][2,2+(iBod-1)*6] - self.tmp['adjmass']
            self.hydroForce['fAddedMass'][3,3+(iBod-1)*6] = 0
            self.hydroForce['fAddedMass'][4,4+(iBod-1)*6] = 0
            self.hydroForce['fAddedMass'][5,5+(iBod-1)*6] = 0
        else:
            self.tmp['fadm'] = np.diag(self.hydroForce['fAddedMass'])
            self.tmp['adjmass'] = sum(self.tmp['fadm'][0:3])*adjMassWeightFun # scalar value for adjMassWeighFun
            self.mass = self.mass + self.tmp['adjmass']
            self.momOfInertia = self.momOfInertia + self.tmp['fadm'][3:6]
            self.hydroForce['fAddedMass'][0,0] = self.hydroForce['fAddedMass'][0,0] - self.tmp['adjmass']
            self.hydroForce['fAddedMass'][1,1] = self.hydroForce['fAddedMass'][1,1] - self.tmp['adjmass']
            self.hydroForce['fAddedMass'][2,2] = self.hydroForce['fAddedMass'][2,2] - self.tmp['adjmass']
            self.hydroForce['fAddedMass'][3,3] = 0
            self.hydroForce['fAddedMass'][4,4] = 0
            self.hydroForce['fAddedMass'][5,5] = 0

    def restoreMassMatrix(self):
        """
        Restore the mass and added-mass matrix back to the original value
        Used copy method to set as new variables
        """
        self.tmp['mass'] = copy(self.mass)
        self.tmp['momOfInertia'] = copy(self.momOfInertia)
        self.tmp['hydroForce_fAddedMass'] = copy(self.hydroForce['fAddedMass'])
        self.mass = copy(self.hydroForce['storage']['mass'])
        self.momOfInertia = copy(self.hydroForce['storage']['momOfInertia'])
        self.hydroForce['fAddedMass'] = copy(self.hydroForce['storage']['fAddedMass'])
        self.hydroForce['storage']['mass'] = copy(self.tmp['mass'])
        self.hydroForce['storage']['momOfInertia'] = copy(self.tmp['momOfInertia'])
        self.hydroForce['storage']['fAddedMass'] = copy(self.tmp['hydroForce_fAddedMass'])
             
    def storeForceAddedMass(self,am_mod,ft_mod):
        """
        Store the modified added mass and total forces history (inputs)
        """
        self.hydroForce['storage']['output_forceAddedMass'] = am_mod
        self.hydroForce['storage']['output_forceTotal'] = ft_mod
 
    def setInitDisp(self, x_rot, ax_rot, ang_rot, addLinDisp):
        """
        Function to set the initial displacement when having initial rotation
        x_rot: rotation point
        ax_rot: axis about which to rotate (must be a normal vector)
        ang_rot: rotation angle in radians
        addLinDisp: initial linear displacement (in addition to the displacement caused by rotation)
        """
        cg = self.cg
        relCoord = cg - x_rot
        rotatedRelCoord = self.rotateXYZ(relCoord,ax_rot,ang_rot)
        newCoord = rotatedRelCoord + x_rot
        linDisp = newCoord-cg
        self.initDisp = {'initLinDisp':(linDisp + addLinDisp), 
                         'initAngularDispAxis': ax_rot, 
                         'initAngularDispAngle': ang_rot}     
    
    def listInfo(self):
        """
        List body info{:f}\n'.format(self.T)
        """
        print('\n\t***** Body Number ', self.hydroData['properties']['body_number'][0][0] ,', Name: ' , self.hydroData['properties']['name'] , ' *****\n')
        print('\tBody CG                          (m) = [',self.hydroData['properties']['cg'][0],']\n')
        print('\tBody Mass                       (kg) = ',self.mass[0][0],' \n')
        print('\tBody Diagonal MOI              (kgm2)= [',self.momOfInertia,']\n')

    def bodyGeo(self,fname):
        """
        Reads mesh file and calculates areas and centroids
        Use trimesh to get geometry values besides using the method from WEC-Sim
        Takes in string parameter called fname. fname is name of stl file

        """
        your_mesh = trimesh.load_mesh(fname)
        self.meshFile = your_mesh
        v = your_mesh.triangles.reshape(int(np.size(your_mesh.triangles)/3),3)
        [vertex,faces] = np.unique(v,return_inverse=True,axis=0)
        face = faces.reshape(int(np.size(faces)/3),3)+1
        self.bodyGeometry['vertex'] = vertex
        self.bodyGeometry['numVertex'] = np.size(vertex,0)
        self.bodyGeometry['face'] = face
        self.bodyGeometry['numFace'] = np.size(face,0)
        self.bodyGeometry['norm'] = your_mesh.face_normals
        self.bodyGeometry['center'] = your_mesh.triangles_center
        self.bodyGeometry['area'] = np.transpose([your_mesh.area_faces])        
        
    def plotStl(self):
        """
        Plots the body's mesh
        networkx is required if you want to use this function

        """
        # Plots the body's mesh
        your_mesh = self.meshFile
        your_mesh.show()
                    
    def checkinputs(self):
        """
        Checks the user inputs

        """
        # hydro data file
        if self.h5File == None and self.nhBody == 0:
            warnings.warn("The hdf5 file does not exist")
        # geometry file
        if self.geometryFile == None:
            warnings.warn("Could not locate and open geometry file")
        

    def noExcitation(self):
        """
        Set exciation force for no excitation case
        Gets used in hydroForcePre for noWave condition
        """
        nDOF = int(self.dof[0])
        self.hydroForce['fExt']['re'] = np.zeros(nDOF)
        self.hydroForce['fExt']['im'] = np.zeros(nDOF)
    
    
    def regExcitation(self,w,waveDir,rho,g):
        """
        Regular wave excitation force
        Used by hydroForcePre       

        """
        nDOF = int(self.dof[0])
        re = self.hydroData['hydro_coeffs']['excitation']['re']*rho*g
        im = self.hydroData['hydro_coeffs']['excitation']['im']*rho*g
        md = self.hydroData['hydro_coeffs']['mean_drift']*rho*g
        self.hydroForce['fExt']['re'] = np.zeros(nDOF)
        self.hydroForce['fExt']['im'] = np.zeros(nDOF)
        self.hydroForce['fExt']['md'] = np.zeros(nDOF)
        for ii in range(nDOF):
            if np.size(self.hydroData['simulation_parameters']['wave_dir']) > 1:
                x = self.hydroData['simulation_parameters']['w'][0]
                y = self.hydroData['simulation_parameters']['wave_dir'][0]
                s1 = interpolate.interp2d(x,y, np.squeeze(re[ii]))# interpolate using interp2d to get interpolation of 2d space
                self.hydroForce['fExt']['re'][ii] = s1(w[0],waveDir)
                s2 = interpolate.interp2d(x,y, np.squeeze(im[ii]))
                self.hydroForce['fExt']['im'][ii] = s2(w[0],waveDir)
                s3 = interpolate.interp2d(x,y, np.squeeze(md[ii]))
                self.hydroForce['fExt']['md'][ii] = s3(w[0],waveDir)
            elif self.hydroData['simulation_parameters']['wave_dir'] == waveDir:
                x = self.hydroData['simulation_parameters']['w'][0]
                s1 = interpolate.CubicSpline(x, np.squeeze(re[ii][0]))# interpolate using CubicSline to get interpolation of spline 3d space
                s2 = interpolate.CubicSpline(x, np.squeeze(im[ii][0]))
                s3 = interpolate.CubicSpline(x, np.squeeze(md[ii][0]))
                self.hydroForce['fExt']['re'][ii] = s1(w)
                self.hydroForce['fExt']['im'][ii] = s2(w)
                self.hydroForce['fExt']['md'][ii] = s3(w)

    def irrExcitation(self,wv,numFreq,waveDir,rho,g):
        """
        Irregular wave excitation force
        Used by hydroForcePre

        """
        nDOF = int(self.dof[0])
        re = self.hydroData['hydro_coeffs']['excitation']['re']*rho*g
        im = self.hydroData['hydro_coeffs']['excitation']['im']*rho*g
        md = self.hydroData['hydro_coeffs']['mean_drift']*rho*g
        self.hydroForce['fExt']['re'] = np.zeros((np.size(waveDir),numFreq,nDOF))
        self.hydroForce['fExt']['im'] = np.zeros((np.size(waveDir),numFreq,nDOF))
        self.hydroForce['fExt']['md'] = np.zeros((np.size(waveDir),numFreq,nDOF))
        for ii in range(nDOF):
            if np.size(self.hydroData['simulation_parameters']['wave_dir']) > 1:
                x = self.hydroData['simulation_parameters']['w'][0]
                y = self.hydroData['simulation_parameters']['wave_dir'][0]
                s1 = interpolate.interp2d(x,y, np.squeeze(re[ii]))# interpolate using interp2d to get interpolation of 2d space
                self.hydroForce['fExt']['re'][:,:,ii] = s1(wv[0],waveDir)
                s2 = interpolate.interp2d(x,y, np.squeeze(im[ii]))
                self.hydroForce['fExt']['im'][:,:,ii] = s2(wv[0],waveDir)
                s3 = interpolate.interp2d(x,y, np.squeeze(md[ii]))
                self.hydroForce['fExt']['md'][:,:,ii] = s3(wv[0],waveDir)
            elif self.hydroData['simulation_parameters']['wave_dir'] == waveDir:
                x = self.hydroData['simulation_parameters']['w'][0]
                s1 = interpolate.CubicSpline(x, np.squeeze(re[ii][0]))# interpolate using CubicSline to get interpolation of spline 3d space
                s2 = interpolate.CubicSpline(x, np.squeeze(im[ii][0]))
                s3 = interpolate.CubicSpline(x, np.squeeze(md[ii][0]))
                self.hydroForce['fExt']['re'][:,:,ii] = s1(wv)
                self.hydroForce['fExt']['im'][:,:,ii] = s2(wv)
                self.hydroForce['fExt']['md'][:,:,ii] = s3(wv)
    
    def userDefinedExcitation(self,waveAmpTime,dt,waveDir,rho,g):
        """
        Calculated User-Defined wave excitation force with non-causal convolution
        Used by hydroForcePre
        
        """
        nDOF = int(self.dof[0])
        kf = self.hydroData['hydro_coeffs']['excitation']['impulse_response_fun']['f']*rho*g
        kt = self.hydroData['hydro_coeffs']['excitation']['impulse_response_fun']['t'][0]
        t =  arange_MATLAB(np.min(kt),np.max(kt)+dt,dt)
        for ii in range(nDOF):
            if np.size(self.hydroData['simulation_parameters']['wave_dir']) > 1:
                y = self.hydroData['simulation_parameters']['wave_dir'][0] 
                s1 = interpolate.interp2d(kt,y, np.squeeze(kf[ii])) # interpolate using interp2d to get interpolation of 2d space
                self.userDefinedExcIRF = s1(t,waveDir)
            elif self.hydroData['simulation_parameters']['wave_dir'] == waveDir:
                s1 = interpolate.CubicSpline(kt, np.squeeze(kf[ii][0])) # interpolate using CubicSline to get interpolation of spline 3d space
                self.userDefinedExcIRF = s1(t)
            else:
                warnings.warn("Default wave direction different from hydro database value. Wave direction (waves.waveDir) should be specified on input file.",DeprecationWarning)
                
            self.hydroForce['userDefinedFe'][:,ii] = np.convolve(waveAmpTime[1],self.userDefinedExcIRF,'valid')*dt
        
        self.hydroForce['fExt']['re'] = np.zeros(nDOF)
        self.hydroForce['fExt']['im'] = np.zeros(nDOF)
        self.hydroForce['fExt']['md'] = np.zeros(nDOF)
    
    
    def constAddedMassAndDamping(self,w,CIkt,rho,B2B):
        """
        Set added mass and damping for a specific frequency
        Used by hydroForcePre
        
        """
        am = self.hydroData['hydro_coeffs']['added_mass']['all']*rho
        rd = self.hydroData['hydro_coeffs']['radiation_damping']['all']*rho
        for i in range(len(self.hydroData['simulation_parameters']['w'][0])):
            rd[:,:,i] *= self.hydroData['simulation_parameters']['w'][0][i]
        # Change matrix size: B2B [6x6n], noB2B [6x6]
        if B2B == 1:
            lenJ = 6*int(self.bodyTotal[0])
            self.hydroForce['fAddedMass'] = np.zeros((6,lenJ))
            self.hydroForce['fDamping'] = np.zeros((6,lenJ))
            self.hydroForce['totDOF']  = np.zeros((6,lenJ))
            for ii in range(6):
                for jj in range(lenJ):
                    s1 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(am[ii,jj,:]))
                    self.hydroForce['fAddedMass'][ii,jj] = s1(w)
                    s2 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(rd[ii,jj,:]))
                    self.hydroForce['fDamping'][ii,jj] = s2(w)
        else: #B2B =2
            nDOF = int(self.dof[0])
            self.hydroForce['fAddedMass'] = np.zeros((nDOF,nDOF))
            self.hydroForce['fDamping'] = np.zeros((nDOF,nDOF))
            self.hydroForce['totDOF']  = np.zeros((nDOF,nDOF))
            for ii in range(nDOF):
                for jj in range(nDOF):
                    jjj = int(self.dof_start[0])-1+jj
                    s1 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(am[ii,jjj,:]))
                    self.hydroForce['fAddedMass'][ii,jj] = s1(w)
                    s2 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(rd[ii,jjj,:]))
                    self.hydroForce['fDamping'][ii,jj] = s2(w)
    
    
    def irfInfAddedMassAndDamping(self,CIkt,CTTime,ssCalc,rho,B2B):
        """
        Set radiation force properties using impulse response function\
        Added mass at infinite frequency
        Convolution integral raditation dampingiBod
        State space formulation
        Used by hydroForcePre
        
        """
        nDOF = int(self.dof[0])
        if B2B == 1:
            LDOF = int(self.bodyTotal[0])*6
        else:
            LDOF = int(self.dof[0])
        
        # Convolution integral formulation
        if B2B == 1:
            self.hydroForce['fAddedMass'] = self.hydroData['hydro_coeffs']['added_mass']['inf_freq']*rho
        else:
            self.hydroForce['fAddedMass'] = self.hydroData['hydro_coeffs']['added_mass']['inf_freq'][:,int(self.dof_start[0])-1:int(self.dof_end[0])]*rho
        
        # Radition IRF
        self.hydroForce['fDamping'] = np.zeros((nDOF,LDOF))
        irfk = self.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['K']*rho
        irft = self.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['t'][0]
        self.hydroForce['irkb'] = np.zeros((len(CTTime),nDOF,LDOF))
        if B2B == 1:
            for ii in range(nDOF):
                for jj in range(LDOF):
                    s1 = interpolate.CubicSpline(irft, np.squeeze(irfk[ii,jj,:]))
                    self.hydroForce['irkb'][:,ii,jj] = s1(CTTime)
        else:
            for ii in range(nDOF):
                for jj in range(LDOF):
                    jjj = int(self.dof_start[0])-1+jj
                    s1 = interpolate.CubicSpline(irft, np.squeeze(irfk[ii,jjj,:]))
                    self.hydroForce['irkb'][:,ii,jj] = s1(CTTime)
                
        # State Space Formulation
        if ssCalc == 1:
            if B2B == 1:
                for ii in range(nDOF):
                    for jj in range(LDOF):
                        arraySize = int(self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'][ii,jj])
                        if ii == 0 and jj == 0: # Begin construction of combined state, input, and output matrices
                            Af = np.zeros((arraySize,arraySize))
                            Bf = np.zeros((arraySize,LDOF))
                            Cf = np.zeros((nDOF,arraySize))
                            Af[:arraySize,:arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            Bf[:arraySize,jj]         = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cf[ii,:arraySize]         = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                        else:
                            Ay = np.size(Af,0)
                            Ax = np.size(Af,1)
                            Af = np.pad(Af, ((0,arraySize),(0,arraySize)), mode='constant', constant_values=0)
                            Af[Ay:Ay+arraySize,Ax:Ax+arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            By = np.size(Bf,0)
                            Bf = np.pad(Bf, ((0,arraySize),(0,0)), mode='constant', constant_values=0)
                            Bf[By:By+arraySize,jj] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cx = np.size(Cf,1)
                            Cf = np.pad(Cf, ((0,0),(0,arraySize)), mode='constant', constant_values=0)
                            Cf[ii,Cx:Cx+arraySize]  = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                self.hydroForce['ssRadf']['D'] = np.zeros((nDOF,LDOF))
            else:
                for ii in range(nDOF):
                    for jj in range(int(self.dof_start[0])-1,int(self.dof_end[0])):
                        jInd = jj-int(self.dof_start[0])+1
                        arraySize = int(self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'][ii,jj])
                        if ii == 0 and jInd == 0: # Begin construction of combined state, input, and output matrices
                            Af = np.zeros((arraySize,arraySize))
                            Bf = np.zeros((arraySize,LDOF))
                            Cf = np.zeros((nDOF,arraySize))
                            Af[:arraySize,:arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            Bf[:arraySize,jInd]       = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cf[ii,:arraySize]         = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                        else:
                            Ay = np.size(Af,0)
                            Ax = np.size(Af,1)
                            Af = np.pad(Af, ((0,arraySize),(0,arraySize)), mode='constant', constant_values=0)
                            Af[Ay:Ay+arraySize,Ax:Ax+arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            By = np.size(Bf,0)
                            Bf = np.pad(Bf, ((0,arraySize),(0,0)), mode='constant', constant_values=0)
                            Bf[By:By+arraySize,jInd] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cx = np.size(Cf,1)
                            Cf = np.pad(Cf, ((0,0),(0,arraySize)), mode='constant', constant_values=0)
                            Cf[ii,Cx:Cx+arraySize]  = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                self.hydroForce['ssRadf']['D'] = np.zeros((nDOF,nDOF))
            self.hydroForce['ssRadf']['A'] = Af
            self.hydroForce['ssRadf']['B'] = Bf
            self.hydroForce['ssRadf']['C'] = Cf*rho
    
    def setMassMatrix(self, rho, nlHydro):
        """
        Sets mass for the special cases of body at equilibrium or fixed
        Used by hydroForcePre
        
        """
        if self.mass == 'equilibrium':
            self.massCalcMethod = self.mass
            if nlHydro == 0:
                self.mass = self.hydroData['properties']['disp_vol'] * rho
            else:
                cg_tmp = self.hydroData['properties']['cg'][0]
                z = np.conj(np.transpose(np.array(self.bodyGeometry['center'])))[2] + cg_tmp[2]
                zr = [0 if x > 0 else x for x in z]
                area = np.conj(np.transpose(self.bodyGeometry['area']))[0]
                av = np.array([area,area,area])*-1*np.conj(np.transpose(self.bodyGeometry['norm']))
                tmp = rho*np.array([zr, zr, zr])*-1*av
                self.mass = sum(tmp[2])
        elif self.mass == 'fixed':
            self.massCalcMethod = self.mass
            self.mass = 999
            self.momOfInertia = [999, 999, 999]
        else:
            self.massCalcMethod = 'user'
        
    def forceAddedMass(self,acc,B2B):
        """
        Calculates and outputs the real added mass force time history
        
        """
        iBod = self.bodyNumber
        fam = np.zeros(np.shape(acc))
        for i in range(6):
            tmp = np.zeros(np.size(acc[:,i]))
            for j in range(6):
                if B2B == 1:
                    jj = (iBod-1)*6+j
                else:
                    jj = j
                iam = self.hydroForce['fAddedMass'][i,jj]
                tmp = tmp + acc[:,j]* iam
            fam[:,i] = tmp
        return(fam)
    
    def rotateXYZ(self,x,ax,t):
        """
        Function to rotate a point about an arbitrary axis
        x: 3-componenet coordiantes
        ax: axis about which to rotate (must be a normal vector)
        t: rotation angle
        xn: new coordinates after rotation
        
        """
        rotMat = np.zeros((3,3))
        rotMat[0,0] = ax[0]*ax[0]*(1-np.cos(t))    + np.cos(t)
        rotMat[0,1] = ax[1]*ax[0]*(1-np.cos(t))    + ax[2]*np.sin(t)
        rotMat[0,2] = ax[2]*ax[0]*(1-np.cos(t))    - ax[1]*np.sin(t)
        rotMat[1,0] = ax[0]*ax[1]*(1-np.cos(t))    - ax[2]*np.sin(t)
        rotMat[1,1] = ax[1]*ax[1]*(1-np.cos(t))    + np.cos(t)
        rotMat[1,2] = ax[2]*ax[1]*(1-np.cos(t))    + ax[0]*np.sin(t)
        rotMat[2,0] = ax[0]*ax[2]*(1-np.cos(t))    + ax[1]*np.sin(t)
        rotMat[2,1] = ax[1]*ax[2]*(1-np.cos(t))    - ax[0]*np.sin(t)
        rotMat[2,2] = ax[2]*ax[2]*(1-np.cos(t))    + np.cos(t)
        xn = np.dot(x,rotMat)
        return(xn)
    
    
    def offsetXYZ(self,verts,x):
        """
        Function to move the position vertices
        
        """
        verts_out = np.zeros(3)#for now change when being tested later
        verts_out[0] = verts[0] + x[0]
        verts_out[1] = verts[1] + x[1]
        verts_out[2] = verts[2] + x[2]
        return verts_out

def arange_MATLAB(start, end, step):
    """
    Change np.arange to have same sequence as MATLAB when step is float
    """
    return step*np.arange(start/step, np.floor(end/step))
    