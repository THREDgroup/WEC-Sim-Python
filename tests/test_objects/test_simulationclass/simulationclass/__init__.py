#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 08:01:43 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
import numpy as np
import datetime
import os
import warnings

class SimulationClass:
    # This class contains WEC-Sim simulation parameters and settings
    def inputProperties(self):
        self.startTime           = 0                                            # (`float`) Simulation start time. Default = ``0`` s
        self.rampTime            = 100                                          # (`float`) Ramp time for wave forcing. Default = ``100`` s
        self.endTime             = []                                           # (`float`) Simulation end time. Default = ``'NOT DEFINED'``
        self.dt                  = 0.1                                          # (`float`) Simulation time step. Default = ``0.1`` s
        self.dtOut               = []                                           # (`float`) Output sampling time. Default = ``dt``
        self.dtNL                = []                                           # (`float`) Sample time to calculate nonlinear forces. Default = ``dt``
        self.dtCITime            = []                                           # (`float`) Sample time to calculate Convolution Integral. Default = ``dt``
        self.dtME                = []                                           # (`float`) Sample time to calculate Morison Element forces. Default = ``dt``
        self.CITime              = 60                                           # (`float`) Convolution integral time. Default = ``60`` s
        self.domainSize          = 200                                          # (`float`) Size of free surface and seabed. This variable is only used for visualization. Default = ``200`` m
        self.ssCalc              = 0                                            # (`integer`) Option for convolution integral or state-space calculation: convolution integral->0, state-space->1. Default = ``0``
        self.mode                = 'normal'                                     # (`string`) Simulation execution mode, 'normal', 'accelerator', 'rapid-accelerator'. Default = ``'normal'``
        self.solver              = 'ode4'                                       # (`string`) PDE solver used by the Simulink/SimMechanics simulation, 'ode4, 'ode45'. Default = ``'ode4'``
        self.numIntMidTimeSteps  = 5                                            # (`integer`) Number of intermediate time steps. Default = ``5`` for ode4 method
        self.autoRateTranBlk     = 'on'                                         # (`string`) Automatically handle rate transition for data transfer, 'on', 'off'. Default = ``'on'``
        self.zeroCrossCont       = 'DisableAll'                                 # (`string`) Disable zero cross control. Default = ``'DisableAll'``
        self.explorer            = 'on'                                         # (`string`) SimMechanics Explorer 'on' or 'off'. Default = ``'on'``
        self.rho                 = 1000                                         # (`float`) Density of water. Default = ``1000`` kg/m^3
        self.g                   = 9.81                                         # (`float`) Acceleration due to gravity. Default = ``9.81`` m/s
        self.nlHydro             = 0                                            # (`integer`) Option for nonlinear hydrohanamics calculation: linear->0, nonlinear->1. Default = ``0``
        self.yawNonLin           = 0                                            # (`integer`) Option for nonlinear yaw calculation linear->0, nonlinear->1 for nonlinear. Default = ``0`` 
        self.yawThresh           = 1                                            # (`float`) Yaw position threshold (in degrees) above which excitation coefficients will be interpolated in non-linear yaw. Default = ``1`` dg
        self.b2b                 = 0                                            # (`integer`) Option for body2body interactions: off->0, on->1. Default = ``0``
        self.paraview            = 0                                            # (`integer`) Option for writing vtp files for paraview visualization, off->0, on->1. Default = ``0``
        self.StartTimeParaview   = 0;                                           # (`float`) Start time for the vtk file of Paraview. Default = ``0``                                    
        self.EndTimeParaview     = 100                                          # (`float`) End time for the vtk file of Paraview. Default = ``0``                                      
        self.dtParaview          = 0.1                                          # (`float`) Timestep for Paraview. Default = ``0.1``         
        self.pathParaviewVideo = 'vtk'                                          # (`string`) Path of the folder for Paraview vtk files. Default = ``'vtk'``     
        self.adjMassWeightFun    = 2                                            # (`integer`) Weighting function for adjusting added mass term in the translational direction. Default = ``2``
        self.mcrCaseFile         = []                                           # (`string`) mat file that contain a list of the multiple conditions runs with given conditions. Default = ``'NOT DEFINED'``  
        self.morisonElement     = 0                                             # (`integer`) Option for Morrison Element calculation: off->0, on->1. Default = ``0``
        self.outputtxt           = 0                                            # (`integer`) Option to save results as ASCII files off->0, on->1. Default = ``0``
        self.reloadH5Data        = 0                                            # (`integer`) Option to re-load hydro data from hf5 file between runs: off->0, on->1. Default = ``0``
        self.pressureDis         = 0                                            # (`integer`) Option to save pressure distribution: off->0, on->1. Default = ``0``

    def internalProperties(self):
        self.version             = '4.0'                                        # (`string`) WEC-Sim-Python version
        self.simulationDate      = datetime.datetime.now()                                    # (`string`) Simulation date and time
        self.outputDir           = 'output'                                     # (`string`) Data output directory name. Default = ``'output'``
        self.time                = 0                                            # (`float`) Simulation time [s]. Default = ``0`` s
        self.inputFile           = 'wecSimInputFile'                            # (`string`) Name of WEC-Sim input file. Default = ``'wecSimInputFile'``
        self.logFile             = []                                           # (`string`) File with run information summary. Default = ``'log'``
        self.caseFile            = []                                           # (`string`) .mat file with all simulation information. Default = dependent
        self.caseDir             = []                                           # (`string`) WEC-Sim case directory. Default = dependent
        self.CIkt                = []                                           # (`integer`) Number of timesteps in the convolution integral length. Default = dependent
        self.maxIt               = []                                           # (`integer`) Total number of simulation time steps. Default = dependent
        self.CTTime              = []                                           # (`float vector`) Convolution integral time series. Default = dependent
        self.numWecBodies        = []                                           # (`integer`) Number of hydrodynamic bodies that comprise the WEC device. Default = ``'NOT DEFINED'``
        self.numPtos             = []                                           # (`integer`) Number of power take-off elements in the model. Default = ``'NOT DEFINED'``
        self.numConstraints      = []                                           # (`integer`) Number of contraints in the wec model. Default = ``'NOT DEFINED'``
        self.numMoorings         = []                                           # (`integer`) Number of moorings in the wec model. Default = ``'NOT DEFINED'``
    
    def __init__(self):
        """
        Initialize Simulation Class
        Returns
        -------
        None.

        """
        self.internalProperties()
        print('WEC-Sim: An open-source code for simulating wave energy converters\n')
        print('Version: ', self.version,'\n\n')
        print('Initializing the Simulation Class...\n')
        global cwd # set current directory as global variabl
        self.caseDir = os.getcwd()
        print('\tCase Dir: ',self.caseDir,' \n')
        self.outputDir = ['.' ,'/', self.outputDir]

    def setupSim(self):
        """
        Set up simulation property for WEC-Sim-Python
        """
        # Sets simulation properties based on values specified in input file
        self.time = np.arange(self.startTime,self.endTime+1,self.dt)
        self.maxIt = np.floor((self.endTime - self.startTime) / self.dt)
        # Set dtOut if it was not specificed in input file
        if self.dtOut is None or self.dtOut < self.dt:
            self.dtOut = self.dt
        
        # Set dtNL if it was not specificed in input file
        if self.dtNL is None or self.dtNL < self.dt:
            self.dtNL = self.dt
        
        # Set dtCITime if it was not specificed in input file
        if self.dtCITime is None or self.dtCITime < self.dt:
            self.dtCITime = self.dt
        
        # Set dtME if it was not specificed in input file
        if self.dtME is None or self.dtME < self.dt:
            self.dtME = self.dt
                    
        self.CTTime = np.arange(0,self.CITime+1,self.dtCITime)            
        self.CIkt = np.size(self.CTTime)
        #self.caseFile = [obj.caseDir filesep 'output' filesep obj.simMechanicsFile(1:end-4) '_matlabWorkspace.mat'];
        #self.logFile = [obj.caseDir filesep 'output' filesep obj.simMechanicsFile(1:end-4) '_simulationLog.txt'];
        os.mkdir(self.outputDir)
        self.getWecSimPythonVer()
    

    def checkinputs(self):
        """
        Checks user input to ensure that ``simu.endTime`` is specified
        """        
        if self.endTime is None:
            warnings.warn('simu.endTime, the simulation end time must be specified in the wecSimInputFile')

    def rhoDensitySetup(self,rho,g):
        """
        Assigns density and gravity values
        
        Parameters
        ------------
          rho : float
              density of the fluid medium (kg/m^3)
          g : float
              gravitational acceleration constant (m/s^2)
        """
        self.rho = rho
        self.g   = g
    

    def listInfo(self,waveTypeNum):
        """
        List simulation info

        """
        print("\nWEC-Sim Simulation Settings:\n");
        # print("\tTime Marching Solver                 = Fourth-Order Runge-Kutta Formula \n")
        print("\tStart Time                     (sec) = ",self.startTime,"\n")
        print("\tEnd Time                       (sec) = ",self.endTime,"\n")
        print("\tTime Step Size                 (sec) = ",self.dt,"\n")
        print("\tRamp Function Time             (sec) = ",self.rampTime,"\n")
        if waveTypeNum > 10:
            
            print("\tConvolution Integral Interval  (sec) = ",self.CITime,"\n")
        
        print("\tTotal Number of Time Steps           = ",self.maxIt,"\n")
    

    def getWecSimPythonVer(self):
        """
        Determines WEC-Sim version used
        """
        
        # ws_exe = which('wecSim');
        # ws_dir = fileparts(ws_exe);
        # git_ver_file = [ws_dir '/../.git/refs/heads/master'];
        # self.version = textread(git_ver_file,'#s')
        
        self.version = 'No git version available'