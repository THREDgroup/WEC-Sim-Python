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

body_2 = bodyClass_V1.BodyClass('rm3.h5')
body_2.bodyNumber = 1
body_2.readH5file()
w = np.conj(np.transpose(np.loadtxt("./testData/body_2_test/w.txt"))) 
waveDir = [0]
CIkt = 201
CTTime = np.conj(np.transpose(np.loadtxt("./testData/body_2_test/CTTime.txt"))) 
numFreq = 500
dt = 0.1
rho = 1000
g = 9.81
waveType = 'irregular'
waveAmpTime = np.conj(np.transpose(np.loadtxt("./testData/body_2_test/waveAmpTime.txt"))) 
iBod = 1 #body number
numBod = [] #later change it to 2 to check
ssCalc = 0
nlHydro = 0
B2B = 0
body_2.bodyNumber = 1
body_2.readH5file()
body_2.hydroStiffness = np.zeros((6, 6))
body_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
body_2.linearDamping = np.zeros((6, 6))
# body_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
body_2.irrExcitation(w,numFreq,waveDir,rho,g)

