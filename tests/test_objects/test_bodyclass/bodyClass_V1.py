# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 07:52:07 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
import h5py
import numpy as np
import numpy.matlib 
import warnings

from numpy.linalg import inv
from scipy import interpolate

class BodyClass:
    def hdf5FileProperties(self):#hdf5 file properties
        """
        Hydrodynamic data from BEM or user defined.
        hydroData in WEC-Sim-Python is dictionary whereas in WEC-Sim it is struct()
        
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
        self.momOfInertia      = []                                            # Moment of inertia [Ixx Iyy Izz] in kg*m^2
        self.cg                = []                                            # Center of gravity [x y z] in meters. For WEC bodies this is given in the h5 file.
        self.cb                = []                                            # Center of buoyancy [x y z] in meters. For WEC bodies this is given in the h5 file.
        self.dispVol           = []                                            # Displaced volume at equilibrium position in meters cubed. For WEC bodies this is given in the h5 file.
        self.dof               = []                                            # Number of DOFs. For WEC bodies this is given in the h5 file. IF not, default is 6
        self.dof_gbm           = []                                            # Number of DOFs for GBM.
        self.dof_start         = []                                            # Index of DOF starts. For WEC bodies this is given in the h5 file. IF not, default is (bodyNumber-1)*6+1
        self.dof_end           = []                                            # Index of DOF ends. For WEC bodies this is given in the h5 file. IF not, default is (bodyNumber-1)*6+6
        self.geometryFile      = 'NONE'                                        # Location of geomtry stl files
        self.viscDrag          = {                                             # Structure defining the viscous (quadratic) drag
                                  'Drag':                 [0,0,0,0,0,0],       # Viscous (quadratic) drag, matrix 6x6
                                  'cd':                   [0,0,0,0,0,0],       # Viscous (quadratic) drag cd coefficient, vector length 6
                                  'characteristicArea':   [0,0,0,0,0,0]}       # Characteristic area for viscous drag, vector length 6
        self.initDisp          = {                                             # Structure defining the initial displacement
                                  'initLinDisp':          [0,0,0],             # Initial displacement of center fo gravity - used for decay tests (format: [displacment in m], default = [0 0 0])
                                  'initAngularDispAxis':  [0,1,0],             # Initial displacement of cog - axis of rotation - used for decay tests (format: [x y z], default = [1 0 0])
                                  'initAngularDispAngle': 0}                   # Initial displacement of cog - Angle of rotation - used for decay tests (format: [radians], default = 0)
        self.hydroStiffness    = [0,0,0,0,0,0]                                 # Hydrostatic stiffness matrix overrides BEMIO definition, matrix 6x6
        self.linearDamping     = [0,0,0,0,0,0]                                 # Linear damping coefficient, matrix size of 6x6
        self.userDefinedExcIRF = []                                            # Excitation IRF from BEMIO used for User-Defined Time-Series
        self.viz               = {                                             # Structure defining visualization properties
                                  'color': [1,1,0],                            # Visualization color for either SimMechanics Explorer or Paraview.
                                  'opacity': 1}                                # Visualization opacity for either SimMechanics Explorer or Paraview.
        self.morisonElement   = {                                              # Structure defining the Morrison Elements
                                 'cd':                 [0,0,0],                # Viscous (quadratic) drag cd, vector length 3
                                 'ca':                 [0,0,0],                # Added mass coefficent for Morrison Element (format [Ca_x Ca_y Ca_z], default = [0 0 0])
                                 'characteristicArea': [0,0,0],                # Characteristic area for Morrison Elements calculations (format [Area_x Area_y Area_z], default = [0 0 0])
                                 'VME':                 0     ,                # Characteristic volume for Morrison Element (default = 0)
                                 'rgME':               [0,0,0]}                # Vector from center of gravity to point of application for Morrison Element (format [X Y Z], default = [0 0 0]).
        self.nhBody            = 0                                             # Flag for non-hydro body.
        self.flexHydroBody     = 0                                             # Flag for flexible body. 
        self.meanDriftForce    = 0                                             # Flag for mean drift force. 0: No 1: from control surface 2: from momentum conservation.
        
    def bodyGeometryFileProperties(self):#body geometry stl file properties
        self.bodyGeometry = {                                                  # Structure defining body's mesh
                             'numFace': [],                                    # Number of faces
                             'numVertex': [],                                  # Number of vertices
                             'vertex': [],                                     # List of vertices
                             'face': [],                                       # List of faces
                             'norm': [],                                       # List of normal vectors
                             'area': [],                                       # List of cell areas
                             'center': []}                                     # List of cell centers
    
    def internalProperties(self,filename): #internal properties
        self.hydroForce        = {'linearHydroRestCoef':[],
                                  'visDrag':[],
                                  'linearDamping':[],
                                  'userDefinedFe':[],
                                  'fExt':{'re':[],
                                          'im':[],
                                          'md':[]},
                                  'fAddedMass':[],
                                  'fDamping':[],
                                  'irkb':[]}                                    # Hydrodynamic forces and coefficients used during simulation.
        self.h5File            = filename                                       # hdf5 file containing the hydrodynamic data
        self.hydroDataBodyNum  = []                                             # Body number within the hdf5 file.
        self.massCalcMethod    = []                                             # Method used to obtain mass: 'user', 'fixed', 'equilibrium'
        self.bodyNumber        = []                                             # bodyNumber in WEC-Sim as defined in the input file. Can be different from the BEM body number.
        self.bodyTotal         = []                                             # Total number of WEC-Sim bodies (body block iterations)
        self.lenJ              = []                                             # Matrices length. 6 for no body-to-body interactions. 6*numBodies if body-to-body interactions.

    def __init__(self,filename):
        self.internalProperties(filename)
        self.hdf5FileProperties()
        self.meanDriftForce = 0

    def readH5file(self):
        f = h5py.File(self.h5File, 'r')
        name = '/body' + str(self.bodyNumber)
        self.cg = np.transpose(np.array(f.get(name + '/properties/cg')))
        self.cb = np.transpose(np.array(f.get(name + '/properties/cb')))
        self.dispVol = np.array(f.get(name + '/properties/disp_vol'))
        self.name = np.string_(np.array(f.get(name + '/properties/name'))).decode("utf-8")
        self.hydroData['simulation_parameters']['scaled'] = np.array(f.get('/simulation_parameters/scaled'))
        self.hydroData['simulation_parameters']['wave_dir'] = np.array(f.get('/simulation_parameters/wave_dir'))
        self.hydroData['simulation_parameters']['water_depth'] = np.string_(np.array(f.get('/simulation_parameters/water_depth'))).decode("utf-8")
        if self.hydroData['simulation_parameters']['water_depth'] != "infinite":
            self.hydroData['simulation_parameters']['water_depth'] = np.array(f.get('/simulation_parameters/water_depth'))
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
        Loads hydroData structure from matlab variable as alternative
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
        HydroForce Pre-processing calculations
        1. Set the linear hydrodynamic restoring coefficient, viscous
            drag, and linear damping matrices
        2. Set the wave excitation force    

        """
        if self.dof_gbm > 0:
            #self.linearDamping = [self.linearDamping zeros(1,self.dof-length(self.linearDamping))]
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

        if  self.viscDrag['Drag'].any() == 1:  #check if self.viscDrag.Drag is defined
            self.hydroForce['visDrag'] = self.viscDrag['Drag']
        else:
            self.hydroForce['visDrag'] = np.diag(0.5*rho*self.viscDrag['cd']*self.viscDrag['characteristicArea'])

        self.hydroForce['linearDamping'] = self.linearDamping
        self.hydroForce['userDefinedFe'] = np.zeros((len(waveAmpTime[1]),int(self.dof[0])))  #initializing userDefinedFe for non imported wave cases
        if waveType == 'noWave':
            self.noExcitation()
            self.constAddedMassAndDamping(w,CIkt,rho,B2B)
        # elif waveType == 'noWaveCIC':
        #     self.noExcitation()
        #     self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        # elif waveType == 'regular':
        #     self.regExcitation(w,waveDir,rho,g)
        #     self.constAddedMassAndDamping(w,CIkt,rho,B2B)
        elif waveType == 'regularCIC':
            self.regExcitation(w,waveDir,rho,g)
            self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        # elif waveType == 'irregular' or waveType == 'spectrumImport':
        #     self.irrExcitation(w,numFreq,waveDir,rho,g)
        #     self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        # elif waveType == 'etaImport':
        #     self.userDefinedExcitation(waveAmpTime,dt,waveDir,rho,g)
        #     self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)

        gbmDOF = self.dof_gbm
        if gbmDOF>0:
            self.hydroForce['gbm']['stiffness']=self.hydroData['gbm']['stiffness']
            self.hydroForce['gbm']['damping']=self.hydroData['gbm']['damping']
            self.hydroForce['gbm']['mass_ff']=[self.hydroForce['fAddedMass'].arange(7,self.dof+1)[self.hydroForce['fAddedMass'].arange(self.dof_start+6,self.dof_end+1)]]+self.hydroData['gbm']['mass']   # need scaling for hydro part
            self.hydroForce['fAddedMass'][7:self.dof+1] = np.zeros(len(np.arange(7,self.dof+1)))
            self.hydroForce['fAddedMass'][(self.dof_start[0]+6):(self.dof_end[0]+1)] = np.zeros(len(np.arange(self.dof_start+6,self.dof_end+1)))
            self.hydroForce['gbm']['mass_ff_inv']=inv(self.hydroForce['gbm']['mass_ff'])
            
            # state-space formulation for solving the GBM
            self.hydroForce['gbm']['state_space']['A'] = [np.zeros(gbmDOF,gbmDOF), np.eye(gbmDOF,gbmDOF)-inv(self.hydroForce['gbm']['mass_ff'])*self.hydroForce['gbm']['stiffness'],-inv(self.hydroForce['gbm']['mass_ff'])*self.hydroForce['gbm']['damping']]    # move to ... hydroForce sector with scaling .         # or create a new fun for all flex parameters
            self.hydroForce['gbm']['state_space']['B'] = np.eye(2*gbmDOF,2*gbmDOF)
            self.hydroForce['gbm']['state_space']['C'] = np.eye(2*gbmDOF,2*gbmDOF)
            self.hydroForce['gbm']['state_space']['D'] = np.zeros(2*gbmDOF,2*gbmDOF)
            self.flexHydroBody = 1
            self.nhBody=0

    # def adjustMassMatrix(obj,adjMassWeightFun,B2B):
    #     """
    #     Merge diagonal term of added mass matrix to the mass matrix
    #     1. Store the original mass and added-mass properties
    #     2. Add diagonal added-mass inertia to moment of inertia
    #     3. Add the maximum diagonal traslational added-mass to body
    #     mass - this is not the correct description

    #     """
    #     iBod = self.bodyNumber
    #     self.hydroForce['storage']['mass'] = self.mass
    #     self.hydroForce['storage']['momOfInertia'] = self.momOfInertia
    #     self.hydroForce['storage']['fAddedMass'] = self.hydroForce['fAddedMass']
    #     if B2B == 1:
    #         tmp['fadm']=diag(self.hydroForce.fAddedMass(:,1+(iBod-1)*6:6+(iBod-1)*6))
    #         tmp['adjmass'] = sum(tmp.fadm(1:3))*adjMassWeightFun
    #         self.mass = self.mass + tmp.adjmass
    #         self.momOfInertia = self.momOfInertia+tmp.fadm(4:6)'
    #         self.hydroForce.fAddedMass(1,1+(iBod-1)*6) = self.hydroForce.fAddedMass(1,1+(iBod-1)*6) - tmp.adjmass
    #         self.hydroForce.fAddedMass(2,2+(iBod-1)*6) = self.hydroForce.fAddedMass(2,2+(iBod-1)*6) - tmp.adjmass
    #         self.hydroForce.fAddedMass(3,3+(iBod-1)*6) = self.hydroForce.fAddedMass(3,3+(iBod-1)*6) - tmp.adjmass
    #         self.hydroForce.fAddedMass(4,4+(iBod-1)*6) = 0
    #         self.hydroForce.fAddedMass(5,5+(iBod-1)*6) = 0
    #         self.hydroForce.fAddedMass(6,6+(iBod-1)*6) = 0
    #     else
    #         tmp.fadm=diag(self.hydroForce.fAddedMass)
    #         tmp.adjmass = sum(tmp.fadm(1:3))*adjMassWeightFun
    #         self.mass = self.mass + tmp.adjmass
    #         self.momOfInertia = self.momOfInertia+tmp.fadm(4:6)
    #         self.hydroForce.fAddedMass(1,1) = self.hydroForce.fAddedMass(1,1) - tmp.adjmass
    #         self.hydroForce.fAddedMass(2,2) = self.hydroForce.fAddedMass(2,2) - tmp.adjmass
    #         self.hydroForce.fAddedMass(3,3) = self.hydroForce.fAddedMass(3,3) - tmp.adjmass
    #         self.hydroForce.fAddedMass(4,4) = 0
    #         self.hydroForce.fAddedMass(5,5) = 0
    #         self.hydroForce.fAddedMass(6,6) = 0
    
    # def restoreMassMatrix(obj):
    #     """
    #     Restore the mass and added-mass matrix back to the original value
    #     """
    #     tmp = struct
    #     tmp.mass = self.mass
    #     tmp.momOfInertia = self.momOfInertia
    #     tmp.hydroForce_fAddedMass = self.hydroForce.fAddedMass
    #     self.mass = self.hydroForce.storage.mass
    #     self.momOfInertia = self.hydroForce.storage.momOfInertia
    #     self.hydroForce.fAddedMass = self.hydroForce.storage.fAddedMass
    #     self.hydroForce.storage = tmp 
        
    # def storeForceAddedMass(obj,am_mod,ft_mod):
    #     """
    #     Store the modified added mass and total forces history (inputs)
    #     """
    #     self.hydroForce.storage.output_forceAddedMass = am_mod
    #     self.hydroForce.storage.output_forceTotal = ft_mod
    
    # def setInitDisp(obj, x_rot, ax_rot, ang_rot, addLinDisp):
    #     """
    #     Function to set the initial displacement when having initial rotation
    #     x_rot: rotation point
    #     ax_rot: axis about which to rotate (must be a normal vector)
    #     ang_rot: rotation angle in radians
    #     addLinDisp: initial linear displacement (in addition to the displacement caused by rotation)
    #     """
    #     cg = self.cg
    #     relCoord = cg - x_rot
    #     rotatedRelCoord = self.rotateXYZ(relCoord,ax_rot,ang_rot)
    #     newCoord = rotatedRelCoord + x_rot
    #     linDisp = newCoord-cg
    #     self.initDisp.initLinDisp= linDisp + addLinDisp
    #     self.initDisp.initAngularDispAxis = ax_rot
    #     self.initDisp.initAngularDispAngle = ang_rot
    
    # def listInfo(obj):
    #     """
    #     List body info
    #     """
    #     fprintf('\n\t***** Body Number #G, Name: #s *****\n',self.hydroData.properties.body_number,self.hydroData.properties.name)
    #     fprintf('\tBody CG                          (m) = [#G,#G,#G]\n',self.hydroData.properties.cg)
    #     fprintf('\tBody Mass                       (kg) = #G \n',self.mass)
    #     fprintf('\tBody Diagonal MOI              (kgm2)= [#G,#G,#G]\n',self.momOfInertia)
     
    # def bodyGeo(obj,fname):
    #     # Reads mesh file and calculates areas and centroids
    #     try
    #         [self.bodyGeometry.vertex, self.bodyGeometry.face, self.bodyGeometry.norm] = import_stl_fast(fname,1,1)
    #     catch
    #         [self.bodyGeometry.vertex, self.bodyGeometry.face, self.bodyGeometry.norm] = import_stl_fast(fname,1,2)
        
    #     self.bodyGeometry.numFace = length(self.bodyGeometry.face)
    #     self.bodyGeometry.numVertex = length(self.bodyGeometry.vertex)
    #     self.checkStl()
    #     self.triArea()
    #     self.triCenter()
    
    
    # def triArea(obj):
    #     # Function to calculate the area of a triangle
    #     points = self.bodyGeometry.vertex
    #     faces = self.bodyGeometry.face
    #     v1 = points(faces(:,3),:)-points(faces(:,1),:)
    #     v2 = points(faces(:,2),:)-points(faces(:,1),:)
    #     av_tmp =  1/2.*(cross(v1,v2))
    #     area_mag = sqrt(av_tmp(:,1).^2 + av_tmp(:,2).^2 + av_tmp(:,3).^2)
    #     self.bodyGeometry.area = area_mag
    
    
    # def checkStl(obj):
    #     # Function to check STL file
    #     tnorm = self.bodyGeometry.norm
    #     #av = zeros(length(area_mag),3)
    #     #av(:,1) = area_mag.*tnorm(:,1)
    #     #av(:,2) = area_mag.*tnorm(:,2)
    #     #av(:,3) = area_mag.*tnorm(:,3)
    #     #if sum(sum(sign(av_tmp))) ~= sum(sum(sign(av)))
    #     #    warning(['The order of triangle vertices in ' self.geometryFile ' do not follow the right hand rule. ' ...
    #     #        'This will causes visualization errors in the SimMechanics Explorer'])
    #     #
    #     norm_mag = sqrt(tnorm(:,1).^2 + tnorm(:,2).^2 + tnorm(:,3).^2)
    #     check = sum(norm_mag)/length(norm_mag)
    #     if check>1.01 || check<0.99
    #         error(['length of normal vectors in ' self.geometryFile ' is not equal to one.'])
        
    
    
    # def triCenter(obj):
    #     #Function to caculate the center coordinate of a triangle
    #     points = self.bodyGeometry.vertex
    #     faces = self.bodyGeometry.face
    #     c = zeros(length(faces),3)
    #     c(:,1) = (points(faces(:,1),1)+points(faces(:,2),1)+points(faces(:,3),1))./3
    #     c(:,2) = (points(faces(:,1),2)+points(faces(:,2),2)+points(faces(:,3),2))./3
    #     c(:,3) = (points(faces(:,1),3)+points(faces(:,2),3)+points(faces(:,3),3))./3
    #     self.bodyGeometry.center = c
    
    
    # def plotStl(obj):
    #     # Plots the body's mesh and normal vectors
    #     c = self.bodyGeometry.center
    #     tri = self.bodyGeometry.face
    #     p = self.bodyGeometry.vertex
    #     n = self.bodyGeometry.norm
    #     figure()
    #     hold on
    #     trimesh(tri,p(:,1),p(:,2),p(:,3))
    #     quiver3(c(:,1),c(:,2),c(:,3),n(:,1),n(:,2),n(:,3))
    
    
    # def checkinputs(obj):
    #     # Checks the user inputs
    #     # hydro data file
    #     if exist(self.h5File,'file')==0 && self.nhBody==0
    #         error('The hdf5 file #s does not exist',self.h5File)
        
    #     # geometry file
    #     if exist(self.geometryFile,'file') == 0
    #         error('Could not locate and open geometry file #s',self.geometryFile)
        

    # def noExcitation(self):
    #     # Set exciation force for no excitation case
    #     nDOF = self.dof
    #     self.hydroForce['fExt']['re']=np.zeros(1,nDOF)
    #     self.hydroForce.fExt.im=zeros(1,nDOF)
    
    
    def regExcitation(self,w,waveDir,rho,g):
        # Regular wave excitation force
        # Used by hydroForcePre
        nDOF = int(self.dof[0])
        re = self.hydroData['hydro_coeffs']['excitation']['re']*rho*g
        im = self.hydroData['hydro_coeffs']['excitation']['im']*rho*g
        md = self.hydroData['hydro_coeffs']['mean_drift']*rho*g
        self.hydroForce['fExt']['re'] = np.zeros(nDOF)
        self.hydroForce['fExt']['im'] = np.zeros(nDOF)
        self.hydroForce['fExt']['md'] = np.zeros(nDOF)
        for ii in range(nDOF):
            if np.size(self.hydroData['simulation_parameters']['wave_dir']) > 1:
                pass
                # [X,Y] = np.meshgrid(self.hydroData['simulation_parameters']['w'], self.hydroData['simulation_parameters']['wave_dir'])
                # self.hydroForce['fExt']['re'][ii] = interp2(X, Y, np.squeeze(re[ii]), w, waveDir)
                # self.hydroForce['fExt']['im'][ii] = interp2(X, Y, np.squeeze(im[ii]), w, waveDir)
                # self.hydroForce['fExt']['md'][ii] = interp2(X, Y, np.squeeze(md[ii]), w, waveDir)
            elif self.hydroData['simulation_parameters']['wave_dir'] == waveDir:
                s1 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(re[ii][0]))
                s2 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(im[ii][0]))
                s3 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(md[ii][0]))
                self.hydroForce['fExt']['re'][ii] = s1(w)
                self.hydroForce['fExt']['im'][ii] = s2(w)
                self.hydroForce['fExt']['md'][ii] = s3(w)
                # all the possible method of inter1
                # tck = interpolate.splrep(x, y)
                # a[ii] = interpolate.splev(w, tck)

                # s =interpolate.make_interp_spline(x, y)
                # a[ii] = s(w) 

                # s = interpolate.interp1d(x, y, "cubic")
                # a[ii] = s(w)         
                # l, r = [(2, 0)], [(2, 0)]  # natural spline boundary conditions
                # s = interpolate.make_interp_spline(x, y, k=3, bc_type=(l, r))
                # a[ii] = s(w)

            
        
    
    
    # def irrExcitation(obj,wv,numFreq,waveDir,rho,g)
    #     # Irregular wave excitation force
    #     # Used by hydroForcePre
    #     nDOF = self.dof
    #     re = self.hydroData.hydro_coeffs.excitation.re(:,:,:) .*rho.*g
    #     im = self.hydroData.hydro_coeffs.excitation.im(:,:,:) .*rho.*g
    #     md = self.hydroData.hydro_coeffs.mean_drift(:,:,:)    .*rho.*g
    #     self.hydroForce.fExt.re=zeros(length(waveDir),numFreq,nDOF)
    #     self.hydroForce.fExt.im=zeros(length(waveDir),numFreq,nDOF)
    #     self.hydroForce.fExt.md=zeros(length(waveDir),numFreq,nDOF)
    #     for ii=1:nDOF
    #         if length(self.hydroData.simulation_parameters.wave_dir) > 1
    #             [X,Y] = meshgrid(self.hydroData.simulation_parameters.w, self.hydroData.simulation_parameters.wave_dir)
    #             self.hydroForce.fExt.re(:,:,ii) = interp2(X, Y, squeeze(re(ii,:,:)), wv, waveDir)
    #             self.hydroForce.fExt.im(:,:,ii) = interp2(X, Y, squeeze(im(ii,:,:)), wv, waveDir)
    #             self.hydroForce.fExt.md(:,:,ii) = interp2(X, Y, squeeze(md(ii,:,:)), wv, waveDir)
    #         elif self.hydroData.simulation_parameters.wave_dir == waveDir
    #             self.hydroForce.fExt.re(:,:,ii) = interp1(self.hydroData.simulation_parameters.w,squeeze(re(ii,1,:)),wv,'spline')
    #             self.hydroForce.fExt.im(:,:,ii) = interp1(self.hydroData.simulation_parameters.w,squeeze(im(ii,1,:)),wv,'spline')
    #             self.hydroForce.fExt.md(:,:,ii) = interp1(self.hydroData.simulation_parameters.w,squeeze(md(ii,1,:)),wv,'spline')
            
        
    
    
    # def userDefinedExcitation(obj,waveAmpTime,dt,waveDir,rho,g)
    #     # Calculated User-Defined wave excitation force with non-causal convolution
    #     # Used by hydroForcePre
    #     nDOF = self.dof
    #     kf = self.hydroData.hydro_coeffs.excitation.impulse_response_fun.f .*rho .*g
    #     kt = self.hydroData.hydro_coeffs.excitation.impulse_response_fun.t
    #     t =  min(kt):dt:max(kt)
    #     for ii = 1:nDOF
    #         if length(self.hydroData.simulation_parameters.wave_dir) > 1
    #             [X,Y] = meshgrid(kt, self.hydroData.simulation_parameters.wave_dir)
    #             kernel = squeeze(kf(ii,:,:))
    #             self.userDefinedExcIRF = interp2(X, Y, kernel, t, waveDir)
    #         elif self.hydroData.simulation_parameters.wave_dir == waveDir
    #             kernel = squeeze(kf(ii,1,:))
    #             self.userDefinedExcIRF = interp1(kt,kernel,min(kt):dt:max(kt))
    #         else
    #             error('Default wave direction different from hydro database value. Wave direction (waves.waveDir) should be specified on input file.')
            
    #         self.hydroForce.userDefinedFe(:,ii) = conv(waveAmpTime(:,2),self.userDefinedExcIRF,'same')*dt
        
    #     self.hydroForce.fExt.re=zeros(1,nDOF)
    #     self.hydroForce.fExt.im=zeros(1,nDOF)
    #     self.hydroForce.fExt.md=zeros(1,nDOF)
    
    
    def constAddedMassAndDamping(self,w,CIkt,rho,B2B):
        # Set added mass and damping for a specific frequency
        # Used by hydroForcePre
        am = self.hydroData['hydro_coeffs']['added_mass']['all']*rho
        rd = self.hydroData['hydro_coeffs']['radiation_damping']['all']*rho
        for i in range(len(self.hydroData['simulation_parameters']['w'][0])):
            rd[:,:,i] *= self.hydroData['simulation_parameters']['w'][0][i]
        # Change matrix size: B2B [6x6n], noB2B [6x6]
        if B2B == 1:
            lenJ = 6*self.bodyTotal
            self.hydroForce['fAddedMass'] = np.zeros((lenJ,6))
            self.hydroForce['fDamping'] = np.zeros((lenJ,6))
            self.hydroForce['totDOF']  = np.zeros((lenJ,6))
            for ii in range(6):
                for jj in range(lenJ):
                    s1 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(am[ii,jj,:]))
                    self.hydroForce['fAddedMass'][ii,jj] = s1(w)
                    s2 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(rd[ii,jj,:]))
                    self.hydroForce['fDamping'][ii,jj] = s2(w)
        else:
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
        # Set radiation force properties using impulse response function
        # Used by hydroForcePre
        # Added mass at infinite frequency
        # Convolution integral raditation dampingiBod
        # State space formulation
        nDOF = int(self.dof[0])
        if B2B == 1:
            LDOF = int(self.bodyTotal[0])*6
        else:
            LDOF = int(self.dof[0])
        
        # Convolution integral formulation
        if B2B == 1:
            self.hydroForce['fAddedMass'] = self.hydroData['hydro_coeffs']['added_mass']['inf_freq']*rho
        else:
            self.hydroForce['fAddedMass'] = self.hydroData['hydro_coeffs']['added_mass']['inf_freq'][:,int(self.dof_start[0]):int(self.dof_end[0]+1)]*rho
        
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
            Af = np.zeros((nDOF,LDOF))
            Bf = np.zeros((nDOF,LDOF))
            Cf = np.zeros((nDOF,LDOF))
            if B2B == 1:
                for ii in range(nDOF):
                    for jj in range(LDOF):
                        arraySize = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'][ii,jj]
                        if ii == 0 and jj == 0: # Begin construction of combined state, input, and output matrices
                            Af[:arraySize,:arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            Bf[:arraySize,jj]         = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cf[ii,:arraySize]         = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                        else:
                            Af[np.size(Af,0)+1:np.size(Af,0)+arraySize,np.size(Af,1)+1:np.size(Af,1)+arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            Bf[np.size(Bf,0)+1:np.size(Bf,0)+arraySize,jj] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cf[ii,np.size(Cf,1)+1:np.size(Cf,1)+arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                self.hydroForce['ssRadf']['D'] = np.zeros((nDOF,LDOF))
            else:
                for ii in range(nDOF):
                    for jj in range(int(self.dof_start[0]),int(self.dof_end[0])):
                        jInd = jj-int(self.dof_start[0])+1
                        arraySize = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'][ii,jj]
                        if ii == 0 and jInd == 0: # Begin construction of combined state, input, and output matrices
                            Af[:arraySize,:arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            Bf[:arraySize,jInd]       = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cf[ii,:arraySize]         = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                        else:
                            Af[np.size(Af,0)+1:np.size(Af,0)+arraySize,np.size(Af,1)+1:np.size(Af,1)+arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            Bf[np.size(Bf,0)+1:np.size(Bf,0)+arraySize,jj] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cf[ii,np.size(Cf,1)+1:np.size(Cf,1)+arraySize]  = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                self.hydroForce['ssRadf']['D'] = np.zeros((nDOF,nDOF))
            self.hydroForce['ssRadf']['A'] = Af
            self.hydroForce['ssRadf']['B'] = Bf
            self.hydroForce['ssRadf']['C'] = Cf*rho
    
    # def setMassMatrix(obj, rho, nlHydro)
    #     # Sets mass for the special cases of body at equilibrium or fixed
    #     # Used by hydroForcePre
    #     if strcmp(self.mass, 'equilibrium')
    #         self.massCalcMethod = self.mass
    #         if nlHydro == 0
    #             self.mass = self.hydroData.properties.disp_vol * rho
    #         else
    #             cg_tmp = self.hydroData.properties.cg
    #             z = self.bodyGeometry.center(:,3) + cg_tmp(3)
    #             z(z>0) = 0
    #             area = self.bodyGeometry.area
    #             av = [area area area] .* -self.bodyGeometry.norm
    #             tmp = rho*[z z z].*-av
    #             self.mass = sum(tmp(:,3))
            
    #     elif strcmp(self.mass, 'fixed')
    #         self.massCalcMethod = self.mass
    #         self.mass = 999
    #         self.momOfInertia = [999 999 999]
    #     else
    #         self.massCalcMethod = 'user'
        
    
    # def fam = forceAddedMass(obj,acc,B2B)
    #     # Calculates and outputs the real added mass force time history
    #     iBod = self.bodyNumber
    #     fam = zeros(size(acc))
    #     for i =1:6
    #         tmp = zeros(length(acc(:,i)),1)
    #         for j =1:6
    #             if B2B == 1
    #                 jj = (iBod-1)*6+j
    #             else
    #                 jj = j
                
    #             iam = self.hydroForce.fAddedMass(i,jj)
    #             tmp = tmp + acc(:,j) .* iam
            
    #         fam(:,i) = tmp
    
    # def xn = rotateXYZ(obj,x,ax,t)
    #     # Function to rotate a point about an arbitrary axis
    #     # x: 3-componenet coordiantes
    #     # ax: axis about which to rotate (must be a normal vector)
    #     # t: rotation angle
    #     # xn: new coordinates after rotation
    #     rotMat = zeros(3)
    #     rotMat(1,1) = ax(1)*ax(1)*(1-cos(t))    + cos(t)
    #     rotMat(1,2) = ax(2)*ax(1)*(1-cos(t))    + ax(3)*sin(t)
    #     rotMat(1,3) = ax(3)*ax(1)*(1-cos(t))    - ax(2)*sin(t)
    #     rotMat(2,1) = ax(1)*ax(2)*(1-cos(t))    - ax(3)*sin(t)
    #     rotMat(2,2) = ax(2)*ax(2)*(1-cos(t))    + cos(t)
    #     rotMat(2,3) = ax(3)*ax(2)*(1-cos(t))    + ax(1)*sin(t)
    #     rotMat(3,1) = ax(1)*ax(3)*(1-cos(t))    + ax(2)*sin(t)
    #     rotMat(3,2) = ax(2)*ax(3)*(1-cos(t))    - ax(1)*sin(t)
    #     rotMat(3,3) = ax(3)*ax(3)*(1-cos(t))    + cos(t)
    #     xn = x*rotMat
    
    
    # def verts_out = offsetXYZ(obj,verts,x)
    #     # Function to move the position vertices
    #     verts_out(:,1) = verts(:,1) + x(1)
    #     verts_out(:,2) = verts(:,2) + x(2)
    #     verts_out(:,3) = verts(:,3) + x(3)
    
"""
    def write_paraview_vtp(obj, t, pos_all, bodyname, model, simdate, hspressure,wavenonlinearpressure,wavelinearpressure)
        # Writes vtp files for visualization with ParaView
        numVertex = self.bodyGeometry.numVertex
        numFace = self.bodyGeometry.numFace
        vertex = self.bodyGeometry.vertex
        face = self.bodyGeometry.face
        cellareas = self.bodyGeometry.area
        for it = 1:length(t)
            # calculate new position
            pos = pos_all(it,:)
            vertex_mod = self.rotateXYZ(vertex,[1 0 0],pos(4))
            vertex_mod = self.rotateXYZ(vertex_mod,[0 1 0],pos(5))
            vertex_mod = self.rotateXYZ(vertex_mod,[0 0 1],pos(6))
            vertex_mod = self.offsetXYZ(vertex_mod,pos(1:3))
            # open file
            filename = ['vtk' filesep 'body' num2str(self.bodyNumber) '_' bodyname filesep bodyname '_' num2str(it) '.vtp']
            fid = fopen(filename, 'w')
            # write header
            fprintf(fid, '<?xml version="1.0"?>\n')
            fprintf(fid, ['<!-- WEC-Sim Visualization using ParaView -->\n'])
            fprintf(fid, ['<!--   model: ' model ' - ran on ' simdate ' -->\n'])
            fprintf(fid, ['<!--   body:  ' bodyname ' -->\n'])
            fprintf(fid, ['<!--   time:  ' num2str(t(it)) ' -->\n'])
            fprintf(fid, '<VTKFile type="PolyData" version="0.1">\n')
            fprintf(fid, '  <PolyData>\n')
            # write body info
            fprintf(fid,['    <Piece NumberOfPoints="' num2str(numVertex) '" NumberOfPolys="' num2str(numFace) '">\n'])
            # write points
            fprintf(fid,'      <Points>\n')
            fprintf(fid,'        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
            for ii = 1:numVertex
                fprintf(fid, '          #5.5f #5.5f #5.5f\n', vertex_mod(ii,:))
            
            clear vertex_mod
            fprintf(fid,'        </DataArray>\n')
            fprintf(fid,'      </Points>\n')
            # write tirangles connectivity
            fprintf(fid,'      <Polys>\n')
            fprintf(fid,'        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
            for ii = 1:numFace
                fprintf(fid, '          #i #i #i\n', face(ii,:)-1)
            
            fprintf(fid,'        </DataArray>\n')
            fprintf(fid,'        <DataArray type="Int32" Name="offsets" format="ascii">\n')
            fprintf(fid, '         ')
            for ii = 1:numFace
                n = ii * 3
                fprintf(fid, ' #i', n)
            
            fprintf(fid, '\n')
            fprintf(fid,'        </DataArray>\n')
            fprintf(fid, '      </Polys>\n')
            # write cell data
            fprintf(fid,'      <CellData>\n')
            # Cell Areas
            fprintf(fid,'        <DataArray type="Float32" Name="Cell Area" NumberOfComponents="1" format="ascii">\n')
            for ii = 1:numFace
                fprintf(fid, '          #i', cellareas(ii))
            
            fprintf(fid, '\n')
            fprintf(fid,'        </DataArray>\n')
            # Hydrostatic Pressure
            if ~isempty(hspressure)
                fprintf(fid,'        <DataArray type="Float32" Name="Hydrostatic Pressure" NumberOfComponents="1" format="ascii">\n')
                for ii = 1:numFace
                    fprintf(fid, '          #i', hspressure.signals.values(it,ii))
                
                fprintf(fid, '\n')
                fprintf(fid,'        </DataArray>\n')
            
            # Non-Linear Froude-Krylov Wave Pressure
            if ~isempty(wavenonlinearpressure)
                fprintf(fid,'        <DataArray type="Float32" Name="Wave Pressure NonLinear" NumberOfComponents="1" format="ascii">\n')
                for ii = 1:numFace
                    fprintf(fid, '          #i', wavenonlinearpressure.signals.values(it,ii))
                
                fprintf(fid, '\n')
                fprintf(fid,'        </DataArray>\n')
            
            # Linear Froude-Krylov Wave Pressure
            if ~isempty(wavelinearpressure)
                fprintf(fid,'        <DataArray type="Float32" Name="Wave Pressure Linear" NumberOfComponents="1" format="ascii">\n')
                for ii = 1:numFace
                    fprintf(fid, '          #i', wavelinearpressure.signals.values(it,ii))
                
                fprintf(fid, '\n')
                fprintf(fid,'        </DataArray>\n')
            
            fprintf(fid,'      </CellData>\n')
            # end file
            fprintf(fid, '    </Piece>\n')
            fprintf(fid, '  </PolyData>\n')
            fprintf(fid, '</VTKFile>')
            # close file
            fclose(fid)
"""