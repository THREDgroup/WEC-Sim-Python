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

body_1 = bodyClass_V1.BodyClass('rm3.h5')
body_1.bodyNumber = 1
w = 0.785398163397448
waveDir = [0]
rho = 1000
g = 9.81
body_1.readH5file()
body_1.regExcitation(w,waveDir,rho,g)

CTTime =  np.array(np.loadtxt("./testData/body_1_test/CTTime.txt"))
ssCalc = 0
CIkt = 601
b2b = 0
# irfk = body_1.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['K']*rho
# irft = body_1.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['t'][0]
# s1 = interpolate.CubicSpline(irft, np.squeeze(irfk[1,1,:]))
# A = s1(CTTime)

body_1.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,b2b)
# a = np.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
# b = a[:,1:2]
