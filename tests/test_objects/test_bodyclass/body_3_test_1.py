#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 03:09:20 2020

@author: logical
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#[18508.7581918053,0.223532780361825,1444829.62495620,6.41939080869329,140467.443724649,0.0412557412840770]
r = [18508.7581918053,0.223532780361825,1444829.62495620,6.41939080869329,140467.443724649,0.0412557412840770]

np.testing.assert_allclose(a, r)

Created on Fri Jun 12 17:00:40 2020

@author: logical
"""
import bodyClass_V1
import numpy as np
from scipy import interpolate
import scipy.io as sio

def readData(file):
    matFile = sio.loadmat(file) 
    keys = list(matFile.keys())[-1]
    datas = matFile[keys]
    return datas

body_3 = bodyClass_V1.BodyClass('oswec.h5')

w = np.conj(np.transpose(readData("./testData/body_3_test/w.mat")))
waveDir = [0,30,90] 
CIkt = 301
CTTime =  readData("./testData/body_3_test/CTTime.mat")[0]
dt = 0.1
rho = 1000
g = 9.81
waveType = 'irregular'
waveAmpTime = np.conj(np.transpose(readData("./testData/body_3_test/waveAmpTime.mat"))) 
iBod = 1 #body number
numBod = [] #later change it to 2 to check
ssCalc = 0
nlHydro = 0
B2B = 0
numFreq = 500


body_3.bodyNumber = 1
body_3.readH5file()
body_3.hydroStiffness = np.zeros((6, 6))
body_3.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
body_3.linearDamping = np.zeros((6, 6))
# body_3.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
#body_3.irrExcitation(w,numFreq,waveDir,rho,g)
body_3.irrExcitation(w,numFreq,waveDir,rho,g)

nDOF = int(body_3.dof[0])
if B2B == 1:
    LDOF = int(body_3.bodyTotal[0])*6
else:
    LDOF = int(body_3.dof[0])

# Convolution integral formulation
if B2B == 1:
    body_3.hydroForce['fAddedMass'] = body_3.hydroData['hydro_coeffs']['added_mass']['inf_freq']*rho
else:
    body_3.hydroForce['fAddedMass'] = body_3.hydroData['hydro_coeffs']['added_mass']['inf_freq'][:,int(body_3.dof_start[0])-1:int(body_3.dof_end[0])]*rho

# Radition IRF
body_3.hydroForce['fDamping'] = np.zeros((nDOF,LDOF))
irfk = body_3.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['K']*rho
irft = body_3.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['t'][0]
body_3.hydroForce['irkb'] = np.zeros((len(CTTime),nDOF,LDOF))
if B2B == 1:
    for ii in range(nDOF):
        for jj in range(LDOF):
            s1 = interpolate.CubicSpline(irft, np.squeeze(irfk[ii,jj,:]))
            body_3.hydroForce['irkb'][:,ii,jj] = s1(CTTime)
else:
    for ii in range(nDOF):
        for jj in range(LDOF):
            jjj = int(body_3.dof_start[0])-1+jj
            s1 = interpolate.CubicSpline(irft, np.squeeze(irfk[ii,jjj,:]))
            body_3.hydroForce['irkb'][:,ii,jj] = s1(CTTime)
        
# State Space Formulation
if ssCalc == 1:
    Af = np.zeros((nDOF,LDOF))
    Bf = np.zeros((nDOF,LDOF))
    Cf = np.zeros((nDOF,LDOF))
    if B2B == 1:
        for ii in range(nDOF):
            for jj in range(LDOF):
                arraySize = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'][ii,jj]
                if ii == 0 and jj == 0: # Begin construction of combined state, input, and output matrices
                    Af[:arraySize,:arraySize] = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                    Bf[:arraySize,jj]         = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                    Cf[ii,:arraySize]         = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                else:
                    Af[np.size(Af,0)+1:np.size(Af,0)+arraySize,np.size(Af,1)+1:np.size(Af,1)+arraySize] = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                    Bf[np.size(Bf,0)+1:np.size(Bf,0)+arraySize,jj] = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                    Cf[ii,np.size(Cf,1)+1:np.size(Cf,1)+arraySize] = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
        body_3.hydroForce['ssRadf']['D'] = np.zeros((nDOF,LDOF))
    else:
        for ii in range(nDOF):
            for jj in range(int(body_3.dof_start[0]),int(body_3.dof_end[0])):
                jInd = jj-int(body_3.dof_start[0])+1
                arraySize = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'][ii,jj]
                if ii == 0 and jInd == 0: # Begin construction of combined state, input, and output matrices
                    Af[:arraySize,:arraySize] = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                    Bf[:arraySize,jInd]       = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                    Cf[ii,:arraySize]         = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                else:
                    Af[np.size(Af,0)+1:np.size(Af,0)+arraySize,np.size(Af,1)+1:np.size(Af,1)+arraySize] = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                    Bf[np.size(Bf,0)+1:np.size(Bf,0)+arraySize,jj] = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                    Cf[ii,np.size(Cf,1)+1:np.size(Cf,1)+arraySize]  = body_3.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
        body_3.hydroForce['ssRadf']['D'] = np.zeros((nDOF,nDOF))
    body_3.hydroForce['ssRadf']['A'] = Af
    body_3.hydroForce['ssRadf']['B'] = Bf
    body_3.hydroForce['ssRadf']['C'] = Cf*rho
