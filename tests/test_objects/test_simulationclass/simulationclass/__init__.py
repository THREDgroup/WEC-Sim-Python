#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 08:01:43 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
import datetime
import os

class SimulationClass:
    # This class contains WEC-Sim simulation parameters and settings
    def inputProperties(self):
        self.simMechanicsFile    = 'NOT DEFINED'                                # (`string`) Simulink/SimMecahnics model file. Default = ``'NOT DEFINED'``
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
        self.saveMat             = 1                                            # (`integer`) Option to save .mat file for each run: off->0, on->1. Default = ``1``
        self.pressureDis         = 0                                            # (`integer`) Option to save pressure distribution: off->0, on->1. Default = ``0``

    def internalProperties(self):
        self.version             = '4.0'                                        # (`string`) WEC-Sim version
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

    def simulationClass(self):
        # This method initializes the ``simulationClass``.
        print('WEC-Sim: An open-source code for simulating wave energy converters\n')
        print('Version: ', self.version,'\n\n')
        print('Initializing the Simulation Class...\n')
        global cwd # set current directory as global variabl
        self.caseDir = os.getcwd()
        print('\tCase Dir: ',self.caseDir,' \n')
        self.outputDir = ['.' ,'/', self.outputDir]
    
"""
    def obj = loadSimMechModel(obj,fName)
        # This method loads the simulink model and sets parameters
        # 
        # Parameters
        # ------------
        #   fname : string
        #       the name of the SimMechanics ``.slx`` file
        #
        
        load_system(fName);
             obj.simMechanicsFile = fName;
             [~,modelName,~] = fileparts(obj.simMechanicsFile);
             set_param(modelName,'Solver',obj.solver,...
             'StopTime',num2str(obj.endTime),...
             'SimulationMode',obj.mode,...
             'StartTime',num2str(obj.startTime),...
             'FixedStep',num2str(obj.dt),...
             'MaxStep',num2str(obj.dt),...
             'AutoInsertRateTranBlk',obj.autoRateTranBlk,...
             'ZeroCrossControl',obj.zeroCrossCont,...
             'SimCompilerOptimization','on',...            
             'ReturnWorkspaceOutputs','off',...
             'SimMechanicsOpenEditorOnUpdate',obj.explorer);
    

    def setupSim(obj)
        # Sets simulation properties based on values specified in input file
        obj.time = obj.startTime:obj.dt:obj.endTime;
        obj.maxIt = floor((obj.endTime - obj.startTime) / obj.dt);
        # Set dtOut if it was not specificed in input file
        if isempty(obj.dtOut) || obj.dtOut < obj.dt
            obj.dtOut = obj.dt;
        
        # Set dtNL if it was not specificed in input file
        if isempty(obj.dtNL) || obj.dtNL < obj.dt
            obj.dtNL = obj.dt;
        
        # Set dtCITime if it was not specificed in input file
        if isempty(obj.dtCITime) || obj.dtCITime < obj.dt
            obj.dtCITime = obj.dt;
        
        # Set dtME if it was not specificed in input file
        if isempty(obj.dtME) || obj.dtME < obj.dt
            obj.dtME = obj.dt;
                    
        obj.CTTime = 0:obj.dtCITime:obj.CITime;            
        obj.CIkt = length(obj.CTTime);
        obj.caseFile = [obj.caseDir filesep 'output' filesep obj.simMechanicsFile(1:end-4) '_matlabWorkspace.mat'];
        obj.logFile = [obj.caseDir filesep 'output' filesep obj.simMechanicsFile(1:end-4) '_simulationLog.txt'];
        mkdir(obj.outputDir)
        obj.getWecSimVer;
    

    def checkinputs(obj)
        # Checks user input to ensure that ``simu.endTime`` is specified and that the SimMechanics model exists
        
        if isempty(obj.endTime)
            error('simu.endTime, the simulation end time must be specified in the wecSimInputFile')
                    
        # Check simMechanics file exists
        if exist(obj.simMechanicsFile,'file') ~= 4
            error('The simMechanics file, #s, does not exist in the case directory',value)
        
        # Remove existing output folder
        if exist(obj.outputDir,'dir') ~= 0
            try
                rmdir(obj.outputDir,'s')
            catch
                warning('The output directory could not be removed. Please close any files in the output directory and try running WEC-Sim again')
            
        
    

    def rhoDensitySetup(obj,rho,g)
        # Assigns density and gravity values
        #
        # Parameters
        # ------------
        #   rho : float
        #       density of the fluid medium (kg/m^3)
        #   g : float
        #       gravitational acceleration constant (m/s^2)
        #
        obj.rho = rho;
        obj.g   = g;
    

    def listInfo(obj,waveTypeNum)
        # Lists simulation info
        fprintf('\nWEC-Sim Simulation Settings:\n');
        #fprintf('\tTime Marching Solver                 = Fourth-Order Runge-Kutta Formula \n')
        fprintf('\tStart Time                     (sec) = #G\n',obj.startTime)
        fprintf('\tEnd Time                       (sec) = #G\n',obj.endTime)
        fprintf('\tTime Step Size                 (sec) = #G\n',obj.dt)
        fprintf('\tRamp Function Time             (sec) = #G\n',obj.rampTime)
        if waveTypeNum > 10
            fprintf('\tConvolution Integral Interval  (sec) = #G\n',obj.CITime)
        
        fprintf('\tTotal Number of Time Steps           = #u \n',obj.maxIt)
    

    def getWecSimVer(obj)
        # Determines WEC-Sim version used
        try
            ws_exe = which('wecSim');
            ws_dir = fileparts(ws_exe);
            git_ver_file = [ws_dir '/../.git/refs/heads/master'];
            obj.version = textread(git_ver_file,'#s');
        catch
            obj.version = 'No git version available';
        
        
    


"""

