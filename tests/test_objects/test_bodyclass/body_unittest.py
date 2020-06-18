#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 08:43:15 2020

@author: logical
"""

import unittest
import numpy as np

from bodyclass import BodyClass

"""
result =  np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/wk.txt"))) 
self.assertIsNone(np.testing.assert_allclose(self.noWaveCIC_3.wT, [8]))
"""

class TestBody(unittest.TestCase):
    
    def setUpClass():
        print("setupClass")
        
    def tearDownClass():
        print("teardownClass")
    
    def setUp(self):
        print("setUp")
        self.body_1 = BodyClass('rm3.h5')
        
    def tearDown(self):
        print("tearDown\n")
    
    # def test_hydroForcePre(self):
    #     w = 0.7854
    #     waveDir = [0]
    #     CIkt = 601
    #     CTTime = np.conj(np.transpose(np.loadtxt("./testData/body_1_test/CTTime.txt"))) 
    #     numFreq = []
    #     dt = 0.1
    #     rho = 1000
    #     g = 9.81
    #     waveType = 'regularCIC'
    #     waveAmpTime = np.conj(np.transpose(np.loadtxt("./testData/body_1_test/waveAmpTime.txt"))) 
    #     iBod = 1 #body number
    #     numBod = [] #later change it to 2 to check
    #     ssCalc = 0
    #     nlHydro = 0
    #     B2B = 0
        
    #     self.dof_gbm = 0
    #     self.hydroStiffness = np.conj(np.transpose(np.loadtxt("./testData/body_1_test/hydroStiffness.txt"))) 
    #     self.hydroForce = {}
    #     self.body_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
    
    def test_regExcitation(self):
        w = 0.7854
        waveDir = [0]
        rho = 1000
        g = 9.81
        self.body_1.bodyNumber = 1
        self.body_1.readH5file()
        #self.body_1.regExcitation(w,waveDir,rho,g)
        
if __name__ == '__main__':
    unittest.main()