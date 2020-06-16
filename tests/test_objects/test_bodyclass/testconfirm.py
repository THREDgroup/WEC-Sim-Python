#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 23:01:27 2020

@author: logical
"""

import numpy as np


w = 0.7854
waveDir = [0]
CIkt = 601
CTTime = np.conj(np.transpose(np.loadtxt("./testData/body_1_test/CTTime.txt"))) 
numFreq = []
dt = 0.1
rho = 1000
g = 9.81
waveType = 'regularCIC'
waveAmpTime = np.conj(np.transpose(np.loadtxt("./testData/body_1_test/waveAmpTime.txt"))) 
iBod = 1 #body number
numBod = [] #later change it to 2 to check
ssCalc = 0
nlHydro = 0
B2B = 0

dof_gbm = 0
hydroStiffness = np.conj(np.transpose(np.loadtxt("./testData/body_1_test/hydroStiffness.txt"))) 
hydroForce = {}

# if dof_gbm > 0:
#     #linearDamping = [linearDamping zeros(1,dof-length(linearDamping))]
#     tmp0 = linearDamping
#     tmp1 = np.size(linearDamping)
#     linearDamping = np.zeros(dof)                
#     linearDamping(1:tmp1(1),1:tmp1(2)) = tmp0 

#     tmp0 = viscDrag.Drag
#     tmp1 = size(viscDrag.Drag)
#     viscDrag.Drag = zeros (dof)                
#     viscDrag.Drag(1:tmp1(1),1:tmp1(2)) = tmp0 
    
#     viscDrag.cd   = [viscDrag.cd   zeros(1,dof-length(viscDrag.cd  ))]
#     viscDrag.characteristicArea = [viscDrag.characteristicArea zeros(1,dof-length(viscDrag.characteristicArea))]

#
# if  any(any(viscDrag.Drag)) == 1  #check if viscDrag.Drag is defined
#     hydroForce.visDrag = viscDrag.Drag
# else
#     hydroForce.visDrag = diag(0.5*rho.*viscDrag.cd.*viscDrag.characteristicArea)

hydroForce['linearDamping'] = linearDamping
hydroForce['userDefinedFe'] = np.zeros(np.size(waveAmpTime[1]),dof)   #initializing userDefinedFe for non imported wave cases
if waveType == 'noWave':
    noExcitation()
    constAddedMassAndDamping(w,CIkt,rho,B2B)
elif waveType == 'noWaveCIC':
    noExcitation()
    irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
elif waveType == 'regular':
    regExcitation(w,waveDir,rho,g)
    constAddedMassAndDamping(w,CIkt,rho,B2B)
elif waveType == 'regularCIC':
    regExcitation(w,waveDir,rho,g)
    irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
elif waveType == 'irregular' or waveType == 'spectrumImport':
    irrExcitation(w,numFreq,waveDir,rho,g)
    irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
elif waveType == 'etaImport':
    userDefinedExcitation(waveAmpTime,dt,waveDir,rho,g)
    irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)

# gbmDOF = dof_gbm
# if (gbmDOF>0)
#     hydroForce.gbm.stiffness=hydroData.gbm.stiffness
#     hydroForce.gbm.damping=hydroData.gbm.damping
#     hydroForce.gbm.mass_ff=hydroForce.fAddedMass(7:dof,dof_start+6:dof_end)+hydroData.gbm.mass   # need scaling for hydro part
#     hydroForce.fAddedMass(7:dof,dof_start+6:dof_end) = 0
#     hydroForce.gbm.mass_ff_inv=inv(hydroForce.gbm.mass_ff)
    
#     # state-space formulation for solving the GBM
#     hydroForce.gbm.state_space.A = [zeros(gbmDOF,gbmDOF),...
#         eye(gbmDOF,gbmDOF)...  # move to ... hydroForce sector with scaling .
#         -inv(hydroForce.gbm.mass_ff)*hydroForce.gbm.stiffness,-inv(hydroForce.gbm.mass_ff)*hydroForce.gbm.damping]             # or create a new fun for all flex parameters
#     hydroForce.gbm.state_space.B = eye(2*gbmDOF,2*gbmDOF)
#     hydroForce.gbm.state_space.C = eye(2*gbmDOF,2*gbmDOF)
#     hydroForce.gbm.state_space.D = zeros(2*gbmDOF,2*gbmDOF)
#     flexHydroBody = 1
#     nhBody=0
