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

body_1 = bodyClass_V1.BodyClass('rm3.h5')
body_1.bodyNumber = 1
w = 0.785398163397448
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
body_1.bodyNumber = 1
body_1.readH5file()
body_1.hydroStiffness = np.zeros((6, 6))
body_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
body_1.linearDamping = np.zeros((6, 6))
body_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
