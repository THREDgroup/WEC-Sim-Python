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

body_8_1 = bodyClass_V1.BodyClass('rm3.h5')
body_8_2 = bodyClass_V1.BodyClass('rm3.h5')

w = np.conj(np.transpose(readData("./testData/body_8_test/w.mat"))) 
waveDir = [0]
CIkt = 601
CTTime = readData("./testData/body_8_test/CTTime.mat")[0]
dt = 0.1
rho = 1000
g = 9.81
waveType = 'regularCIC'
waveAmpTime = np.conj(np.transpose(readData("./testData/body_8_test/waveAmpTime.mat"))) 
iBod = 1 #body number
numBod = [] #later change it to 2 to check
ssCalc = 1
nlHydro = 0
B2B = 0
numFreq = []
body_8_1.bodyNumber = 1
body_8_1.bodyTotal = [2]
body_8_1.readH5file()
body_8_1.hydroStiffness = np.zeros((6, 6))
body_8_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
body_8_1.linearDamping = np.zeros((6, 6))
#body_8_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)

nDOF = int(body_8_1.dof[0])
if B2B == 1:
    LDOF = int(body_8_1.bodyTotal[0])*6
else:
    LDOF = int(body_8_1.dof[0])

# Convolution integral formulation
if B2B == 1:
    body_8_1.hydroForce['fAddedMass'] = body_8_1.hydroData['hydro_coeffs']['added_mass']['inf_freq']*rho
else:
    body_8_1.hydroForce['fAddedMass'] = body_8_1.hydroData['hydro_coeffs']['added_mass']['inf_freq'][:,int(body_8_1.dof_start[0])-1:int(body_8_1.dof_end[0])]*rho

# Radition IRF
body_8_1.hydroForce['fDamping'] = np.zeros((nDOF,LDOF))
irfk = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['K']*rho
irft = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['t'][0]
body_8_1.hydroForce['irkb'] = np.zeros((len(CTTime),nDOF,LDOF))
if B2B == 1:
    for ii in range(nDOF):
        for jj in range(LDOF):
            s1 = interpolate.CubicSpline(irft, np.squeeze(irfk[ii,jj,:]))
            body_8_1.hydroForce['irkb'][:,ii,jj] = s1(CTTime)
else:
    for ii in range(nDOF):
        for jj in range(LDOF):
            jjj = int(body_8_1.dof_start[0])-1+jj
            s1 = interpolate.CubicSpline(irft, np.squeeze(irfk[ii,jjj,:]))
            body_8_1.hydroForce['irkb'][:,ii,jj] = s1(CTTime)
        
# State Space Formulation
if ssCalc == 1:
    if B2B == 1:
        for ii in range(nDOF):
            for jj in range(LDOF):
                arraySize = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'][ii,jj]
                if ii == 0 and jj == 0: # Begin construction of combined state, input, and output matrices
                    Af = np.zeros((arraySize,arraySize))
                    Bf = np.zeros((arraySize,LDOF))
                    Cf = np.zeros((nDOF,arraySize))
                    Af[:arraySize,:arraySize] = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                    Bf[:arraySize,jj]         = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                    Cf[ii,:arraySize]         = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                else:
                    Ay = np.size(Af,0)
                    Ax = np.size(Af,1)
                    Af = np.pad(Af, ((0,arraySize),(0,arraySize)), mode='constant', constant_values=0)
                    Af[Ay:Ay+arraySize,Ax:Ax+arraySize] = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                    By = np.size(Bf,0)
                    Bf = np.pad(Bf, ((0,arraySize),(0,0)), mode='constant', constant_values=0)
                    Bf[By:By+arraySize,jj] = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                    Cx = np.size(Cf,1)
                    Cf = np.pad(Cf, ((0,0),(0,arraySize)), mode='constant', constant_values=0)
                    Cf[ii,Cx:Cx+arraySize]  = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                    # pass
                    # Af[np.size(Af,0)+1:np.size(Af,0)+arraySize,np.size(Af,1)+1:np.size(Af,1)+arraySize] = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                    # Bf[np.size(Bf,0)+1:np.size(Bf,0)+arraySize,jj] = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                    # Cf[ii,np.size(Cf,1)+1:np.size(Cf,1)+arraySize] = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
        body_8_1.hydroForce['ssRadf']['D'] = np.zeros((nDOF,LDOF))
    else:
        for ii in range(nDOF):
            for jj in range(int(body_8_1.dof_start[0])-1 ,int(body_8_1.dof_end[0])):
                jInd = jj-int(body_8_1.dof_start[0])+1
                arraySize = int(body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'][ii,jj])
                #print(arraySize)
                if ii == 0 and jInd == 0: # Begin construction of combined state, input, and output matrices
                    Af = np.zeros((arraySize,arraySize))
                    Bf = np.zeros((arraySize,LDOF))
                    Cf = np.zeros((nDOF,arraySize))
                    Af[:arraySize,:arraySize] = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                    Bf[:arraySize,jInd]       = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                    Cf[ii,:arraySize]         = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                else:
                    Ay = np.size(Af,0)
                    Ax = np.size(Af,1)
                    Af = np.pad(Af, ((0,arraySize),(0,arraySize)), mode='constant', constant_values=0)
                    Af[Ay:Ay+arraySize,Ax:Ax+arraySize] = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                    By = np.size(Bf,0)
                    Bf = np.pad(Bf, ((0,arraySize),(0,0)), mode='constant', constant_values=0)
                    Bf[By:By+arraySize,jj] = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                    Cx = np.size(Cf,1)
                    Cf = np.pad(Cf, ((0,0),(0,arraySize)), mode='constant', constant_values=0)
                    Cf[ii,Cx:Cx+arraySize]  = body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
        body_8_1.hydroForce['ssRadf']['D'] = np.zeros((nDOF,nDOF))
    body_8_1.hydroForce['ssRadf']['A'] = Af
    body_8_1.hydroForce['ssRadf']['B'] = Bf
    body_8_1.hydroForce['ssRadf']['C'] = Cf*rho

AAA =body_8_1.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][0,0,:2,0]
# result1 = readData("./testData/body_8_test/body1_linearHydroRestCoef.mat")
# assertIsNone(np.testing.assert_allclose(body_8_1.hydroForce['linearHydroRestCoef'], result1))
# result2 = readData("./testData/body_8_test/body1_visDrag.mat")
# assertIsNone(np.testing.assert_allclose(body_8_1.hydroForce['visDrag'], result2))
# result3 = readData("./testData/body_8_test/body1_linearDamping.mat")
# assertIsNone(np.testing.assert_allclose(body_8_1.hydroForce['linearDamping'], result3))
# result4 = readData("./testData/body_8_test/body1_userDefinedFe.mat")
# assertIsNone(np.testing.assert_allclose(body_8_1.hydroForce['userDefinedFe'], result4)) 
# result5 = readData("./testData/body_8_test/body1_re.mat")[0]
# assertIsNone(np.testing.assert_allclose(body_8_1.hydroForce['fExt']['re'], result5))
# result6 = readData("./testData/body_8_test/body1_im.mat")[0]
# assertIsNone(np.testing.assert_allclose(body_8_1.hydroForce['fExt']['im'], result6))
# result7 = readData("./testData/body_8_test/body1_md.mat")[0]
# assertIsNone(np.testing.assert_allclose(body_8_1.hydroForce['fExt']['md'], result7))
# result8 = readData("./testData/body_8_test/body1_fAddedMass.mat")
# assertIsNone(np.testing.assert_allclose(body_8_1.hydroForce['fAddedMass'], result8))
# result9 = readData("./testData/body_8_test/body1_fDamping.mat")
# assertIsNone(np.testing.assert_allclose(body_8_1.hydroForce['fDamping'], result9))
# result10 = readData("./testData/body_8_test/body1_irkb.mat")
# assertIsNone(np.testing.assert_allclose(body_8_1.hydroForce['irkb'], result10))

result11 = readData("./testData/body_8_test/body1_A.mat")
np.testing.assert_allclose(body_8_1.hydroForce['ssRadf']['A'], result11)


result12 = readData("./testData/body_8_test/body1_B.mat")
np.testing.assert_allclose(body_8_1.hydroForce['ssRadf']['B'], result12)






result13 = readData("./testData/body_8_test/body1_C.mat")
np.testing.assert_allclose(body_8_1.hydroForce['ssRadf']['C'], result13)
result14 = readData("./testData/body_8_test/body1_D.mat")
np.testing.assert_allclose(body_8_1.hydroForce['ssRadf']['D'], result14)


# iBod = 2 #body number
# body_8_2.bodyNumber = 2
# body_8_2.bodyTotal = [2]
# body_8_2.readH5file()
# body_8_2.hydroStiffness = np.zeros((6, 6))
# body_8_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
# body_8_2.linearDamping = np.zeros((6, 6))
# body_8_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
# result1 = readData("./testData/body_8_test/body2_linearHydroRestCoef.mat")
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['linearHydroRestCoef'], result1))
# result2 = readData("./testData/body_8_test/body2_visDrag.mat")
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['visDrag'], result2))
# result3 = readData("./testData/body_8_test/body2_linearDamping.mat")
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['linearDamping'], result3))
# result4 = readData("./testData/body_8_test/body2_userDefinedFe.mat")
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['userDefinedFe'], result4)) 
# result5 = readData("./testData/body_8_test/body2_re.mat")[0]
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['fExt']['re'], result5))
# result6 = readData("./testData/body_8_test/body2_im.mat")[0]
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['fExt']['im'], result6))
# result7 = readData("./testData/body_8_test/body2_md.mat")[0]
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['fExt']['md'], result7))
# result8 = readData("./testData/body_8_test/body2_fAddedMass.mat")
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['fAddedMass'], result8))
# result9 = readData("./testData/body_8_test/body2_fDamping.mat")
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['fDamping'], result9))
# result10 = readData("./testData/body_8_test/body2_irkb.mat")
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['irkb'], result10))
# result11 = readData("./testData/body_8_test/body2_A.mat")
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['ssRadf']['A'], result11))
# result12 = readData("./testData/body_8_test/body2_B.mat")
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['ssRadf']['B'], result12))
# result13 = readData("./testData/body_8_test/body2_C.mat")
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['ssRadf']['C'], result13))
# result14 = readData("./testData/body_8_test/body2_D.mat")
# assertIsNone(np.testing.assert_allclose(body_8_2.hydroForce['ssRadf']['D'], result14))

