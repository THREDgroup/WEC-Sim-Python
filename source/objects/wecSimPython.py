#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 08:01:43 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
### WEC-Sim-Python run file
##
## Following Classes are required to be defined in the WEC-Sim input file:
##
## simu = simulationClass();                                               - To Create the Simulation Variable
##
## waves = waveClass('<wave type');                                       - To create the Wave Variable and Specify Type
##
## body(<body number>) = bodyClass('<hydrodynamics data file name>.h5');   - To initialize bodyClass:
##
## constraint(<constraint number>) = constraintClass('<Constraint name>'); - To initialize constraintClass
##
## pto(<PTO number>) = ptoClass('<PTO name>');                             - To initialize ptoClass
##
## mooring(<mooring number>) = mooringClass('<Mooring name>');             - To initialize mooringClass (only needed when mooring blocks are used)
##
##
"""
start
"""
import os
import warnings
import numpy as np

import wecSimPythonInputFile

cwd = os.getcwd()
lastDir = os.path.split(cwd)[1]
if lastDir == 'WEC-Sim-Python': # check if we are in desired directory
    cwd = cwd + '/source'

simu = wecSimPythonInputFile.setUpSimu()
waves = wecSimPythonInputFile.setUpWaves()
body = wecSimPythonInputFile.setUpBody()
constraint = wecSimPythonInputFile.setUpConstraint()
pto = wecSimPythonInputFile.setUpPto()


## Read input file
"""
try print('wecSimMCR Case #g\n',imcr); end 
print('\nWEC-Sim Read Input File ...   \n')
evalc('wecSimInputFile');
# Read Inputs for Multiple Conditions Run
if exist('mcr','var') == 1
    for n=1:length(mcr.cases(1,:))
        if iscell(mcr.cases)
            eval([mcr.header{n} '= mcr.cases{imcr,n};']);
        else
            eval([mcr.header{n} '= mcr.cases(imcr,n);']);
        end
    end; clear n combine;
end
"""
# Waves and Simu: check inputs
waves.checkinputs;
simu.checkinputs;
"""
# Constraints: count & set orientation
if exist('constraint','var') == 1
    simu.numConstraints = length(constraint(1,:));
    for ii = 1:simu.numConstraints
        constraint(ii).constraintNum = ii;
        constraint(ii).setOrientation();
    end; clear ii
end
# PTOs: count & set orientation
if exist('pto','var') == 1
    simu.numPtos = length(pto(1,:));
    for ii = 1:simu.numPtos
        pto(ii).ptoNum = ii;
        pto(ii).setOrientation();
    end; clear ii
end
# Mooring Configuration: count
if exist('mooring','var') == 1
    simu.numMoorings = length(mooring(1,:));
    for ii = 1:simu.numMoorings
        mooring(ii).mooringNum = ii;
        mooring(ii).setLoc;
    end; clear ii
end
"""
# Bodies: count, check inputs, read hdf5 file
numHydroBodies = 0;
for ii in range(len(body)):
    body[ii].bodyNumber = ii + 1
    try: 
        body[ii].nhBody
    except AttributeError:
        body[ii].nhBody = 0
    if body[ii].nhBody == 0:
        numHydroBodies = numHydroBodies + 1
    else:
        body[ii].massCalcMethod = 'user';

simu.numWecBodies = numHydroBodies
del numHydroBodies

try:
    mcr,var,imcr
except NameError:
    mcr,var,imcr = 0,0,0
for ii in range(simu.numWecBodies):
    body[ii].checkinputs
    #Determine if hydro data needs to be reloaded from h5 file, or if hydroData
    # was stored in memory from a previous run.
    if mcr != 0 and var != 0 and simu.reloadH5Data == 0 and imcr > 1:
        body[ii].loadHydroData(hydroData[ii])
    else:
        # check for correct h5 file
        #h5Info = cwd +'/'+ body[ii].h5File
        h5Info = os.path.abspath(os.path.realpath(body[ii].h5File))
        if os.path.getsize(h5Info) < 1000:
            warnings.warn('This is not the correct *.h5 file. Please install git-lfs to access the correct *.h5 file, or run \hydroData\bemio.m to generate a new *.h5 file')
        del h5Info
        body[ii].readH5file()
    body[ii].bodyTotal = simu.numWecBodies
    if simu.b2b == 1:
        body[ii].lenJ = np.zeros(6*body[ii].bodyTotal)
    else:
        body[ii].lenJ = np.zeros(6)

# # PTO-Sim: read input, count
# if exist('./ptoSimInputFile.m','file') == 2
#     ptoSimInputFile
#     ptosim.countblocks;
# end

# if simu.yawNonLin==1 && simu.yawThresh==1;
#     warning(['yawNonLin using (default) 1 dg interpolation threshold.' newline 'Ensure this is appropriate for your geometry'])
# end

# toc

# ## Pre-processing start
# tic
print('\nWEC-Sim-Python Pre-processing ...   \n')

## HydroForce Pre-Processing: Wave Setup & HydroForcePre.
# simulation setup
simu.setupSim()

try: 
    waves.wavegauge1loc, waves.wavegauge2loc, waves.wavegauge3loc
except AttributeError:
    waves.wavegauge1loc = [0,0]
    waves.wavegauge2loc = [0,0]
    waves.wavegauge3loc = [0,0]

# wave setup
waves.waveSetup(body[0].hydroData['simulation_parameters']['w'][0], body[0].hydroData['simulation_parameters']['water_depth'], simu.rampTime, simu.dt, simu.maxIt, simu.g, simu.rho, simu.endTime)

# Check that waveDir and freq are within range of hydro data
try: 
    waves.waveDir
except AttributeError:
    waves.waveDir = [0]
if  np.min(waves.waveDir) <  np.min(body[0].hydroData['simulation_parameters']['wave_dir']) or np.max(waves.waveDir) >  np.max(body[0].hydroData['simulation_parameters']['wave_dir']):
    warnings.warn('waves.waveDir outside of range of available hydro data')
if waves.wType !='etaImport' and waves.wType != 'noWave' and waves.wType != 'noWaveCIC':
    if  np.min(waves.w) < np.min(body[0].hydroData['simulation_parameters']['w']) or np.max(waves.w) >  np.max(body[0].hydroData['simulation_parameters']['w']):
        warnings.warn('waves.w outside of range of available hydro data')


# # Non-linear hydro
# if (simu.nlHydro >0) || (simu.paraview == 1)
#     for kk = 1:length(body(1,:))
#         body(kk).bodyGeo(body(kk).geometryFile)
#     end; clear kk
# end

# hydroForcePre
try: 
    waves.numFreq
except AttributeError:
    waves.numFreq = []

for kk in range(simu.numWecBodies):
    try: 
        body[kk].hydroStiffness
    except AttributeError:
        body[kk].hydroStiffness = np.zeros((6, 6))
    try: 
        body[kk].viscDrag 
    except AttributeError:
        body[kk].viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
    try: 
        body[kk].linearDamping 
    except AttributeError:
        body[kk].linearDamping = np.zeros((6, 6))
    try: 
        body[kk].mass
    except AttributeError:
        body[kk].mass = 'equilibrium'
        

for kk in range(simu.numWecBodies):
    body[kk].hydroForcePre(waves.w,waves.waveDir,simu.CIkt,simu.CTTime,waves.numFreq,simu.dt,
        simu.rho,simu.g,waves.wType,waves.waveAmpTime,kk,simu.numWecBodies,simu.ssCalc,simu.nlHydro,simu.b2b)

# Check CITime
if waves.typeNum != 0 and waves.typeNum != 10:
    for iBod in range(simu.numWecBodies):
        if simu.CITime > np.max(body[iBod].hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['t']):
            warnings.warn('simu.CITime is larger than the length of the IRF')


# Check that the hydro data for each body is given for the same frequencies
for ii in range(simu.numWecBodies):
    if np.size(body[1].hydroData['simulation_parameters']['w']) != np.size(body[ii].hydroData['simulation_parameters']['w']):
        warnings.warn('BEM simulations for each body must have the same number of frequencies')
    else:
        for jj in range(np.size(body[1].hydroData['simulation_parameters']['w'])):
            if body[1].hydroData['simulation_parameters']['w'].all() != body[ii].hydroData['simulation_parameters']['w'].all():
                warnings.warn('BEM simulations must be run with the same frequencies.')

# check for etaImport with nlHydro
if waves.wType == 'etaImport' and simu.nlHydro == 1:
    warnings.warn('Cannot run WEC-Sim with non-linear hydro (simu.nlHydro) and "etaImport" wave type')

# check for etaImport with morisonElement
if waves.wType == 'etaImport' and simu.morisonElement == 1:
    warnings.warn('Cannot run WEC-Sim with Morrison Element (simu.morisonElement) and "etaImport" wave type')

# ## Set variant subsystems options
# nlHydro = simu.nlHydro;
# yawNonLin=simu.yawNonLin;
# sv_linearHydro=Simulink.Variant('nlHydro==0');
# sv_nonlinearHydro=Simulink.Variant('nlHydro>0');
# sv_meanFS=Simulink.Variant('nlHydro<2');
# sv_instFS=Simulink.Variant('nlHydro==2');
# # Morrison Element
# morisonElement = simu.morisonElement;
# sv_MEOff=Simulink.Variant('morisonElement==0');
# sv_MEOn=Simulink.Variant('morisonElement==1');
# # Radiation Damping
# if waves.typeNum==0 || waves.typeNum==10 #'noWave' & 'regular'
#     radiation_option = 1;
# elseif simu.ssCalc == 1
#     radiation_option = 3;
# else
#     radiation_option = 2;
# end
# sv_constantCoeff=Simulink.Variant('radiation_option==1');
# sv_convolution=Simulink.Variant('radiation_option==2');
# sv_stateSpace=Simulink.Variant('radiation_option==3');
# # Wave type
# typeNum = waves.typeNum;
# sv_noWave=Simulink.Variant('typeNum<10');
# sv_regularWaves=Simulink.Variant('typeNum>=10 && typeNum<20 && simu.yawNonLin~=1');
# sv_regularWavesNonLinYaw=Simulink.Variant('typeNum>=10 && typeNum<20 && simu.yawNonLin==1');
# sv_irregularWaves=Simulink.Variant('typeNum>=20 && typeNum<30 && simu.yawNonLin~=1');
# sv_irregularWavesNonLinYaw=Simulink.Variant('typeNum>=20 && typeNum<30 && simu.yawNonLin==1');
# sv_udfWaves=Simulink.Variant('typeNum>=30');
# # Body2Body
# B2B = simu.b2b;
# sv_noB2B=Simulink.Variant('B2B==0');
# sv_B2B=Simulink.Variant('B2B==1');
# numBody=simu.numWecBodies;
# # nonHydroBody
# for ii=1:length(body(1,:))
#     eval(['nhbody_' num2str(ii) ' = body(ii).nhBody;'])
#     eval(['sv_b' num2str(ii) '_hydroBody = Simulink.Variant(''nhbody_' num2str(ii) '==0'');'])
#     eval(['sv_b' num2str(ii) '_nonHydroBody = Simulink.Variant(''nhbody_' num2str(ii) '==1'');'])
# #    eval(['sv_b' num2str(ii) '_flexBody = Simulink.Variant(''nhbody_' num2str(ii) '==2'');'])
# #    eval(['sv_b' num2str(ii) '_rigidBody = Simulink.Variant(''nhbody_' num2str(ii) '==0'');'])
# end; clear ii


# ## End Pre-Processing and Output All the Simulation and Model Setting
# toc
# simu.listInfo(waves.typeNum);
# waves.listInfo
# print('\nList of Body: ');
# print('Number of Bodies = #u \n',simu.numWecBodies)
# for i = 1:simu.numWecBodies
#     body(i).listInfo
# end;  clear i
# print('\nList of PTO(s): ');
# if (exist('pto','var') == 0)
#     print('No PTO in the system\n')
# else
#     print('Number of PTOs = #G \n',length(pto(1,:)))
#     for i=1:simu.numPtos
#         pto(i).listInfo
#     end; clear i
# end
# print('\nList of Constraint(s): ');
# if (exist('constraint','var') == 0)
#     print('No Constraint in the system\n')
# else
#     print('Number of Constraints = #G \n',length(constraint(1,:)))
#     for i=1:simu.numConstraints
#         constraint(i).listInfo
#     end; clear i
# end
# print('\n')


# ## Load simMechanics file & Run Simulation
# tic
# print('\nSimulating the WEC device defined in the SimMechanics model #s...   \n',simu.simMechanicsFile)
# # Modify some stuff for simulation
# for iBod = 1:simu.numWecBodies
#     body(iBod).adjustMassMatrix(simu.adjMassWeightFun,simu.b2b);
# end; clear iBod
# warning('off','Simulink:blocks:TDelayTimeTooSmall');
# warning('off','Simulink:blocks:BusSelDupBusCreatorSigNames');
# warning('off','MATLAB:loadlibrary:FunctionNotFound');
# warning('off','MATLAB:loadlibrary:parsewarnings');
# warning('off','Simulink:blocks:DivideByZero');
# warning('off','sm:sli:setup:compile:SteadyStateStartNotSupported')
# set_param(0, 'ErrorIfLoadNewModel', 'off')
# # run simulation
# simu.loadSimMechModel(simu.simMechanicsFile);
# sim(simu.simMechanicsFile, [], simset('SrcWorkspace','parent'));
# # Restore modified stuff
# clear nlHydro sv_linearHydro sv_nonlinearHydro ssCalc radiation_option sv_convolution sv_stateSpace sv_constantCoeff typeNum B2B sv_B2B sv_noB2B;
# clear nhbod* sv_b* sv_noWave sv_regularWaves sv_irregularWaves sv_udfWaves sv_instFS sv_meanFS sv_MEOn sv_MEOff morisonElement flexHydrobody_* sv_irregularWavesNonLinYaw sv_regularWavesNonLinYaw yawNonLin numBody;
# toc

# tic
# ## Post processing and Saving Results
# postProcess
# # User Defined Post-Processing
# if exist('userDefinedFunctions.m','file') == 2
#     userDefinedFunctions;
# end
# # ASCII files
# if simu.outputtxt==1
#     output.writetxt();
# end
# paraViewVisualization

# ## Save files
# clear ans table tout;
# toc
# diary off
# #movefile('simulation.log',simu.logFile)
# if simu.saveMat==1
#     save(simu.caseFile,'-v7.3')
# end
