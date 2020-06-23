#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 08:43:15 2020

@author: logical
"""

import unittest
import numpy as np
import scipy.io as sio

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
    
    def test_hydroForcePre(self):
        # RM3 example from WEC-Sim
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
        self.body_1.bodyNumber = 1
        self.body_1.readH5file()
        self.body_1.hydroStiffness = np.zeros((6, 6))
        self.body_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_1.linearDamping = np.zeros((6, 6))
        self.body_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = np.loadtxt("./testData/body_1_test/linearHydroRestCoef.txt")
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['linearHydroRestCoef'], result1))
        result2 = np.zeros((6,6))
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['visDrag'], result2))
        result3 = np.zeros((6,6))
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['linearDamping'], result3))
        result4 = np.loadtxt("./testData/body_1_test/userDefinedFe.txt")
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['userDefinedFe'], result4)) 
        result5 = [18508.7581918053,0.223532780361825,1444829.62495620,6.41939080869329,140467.443724649,0.0412557412840770]
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fExt']['re'], result5))
        result6 = [546623.415011964,-0.0428753721889769,447276.984830190,0.273118282507816,4146814.04776144,0.0138224113420678]
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fExt']['im'], result6))
        result7 = [0,0,0,0,0,0]
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fExt']['md'], result7))
        result8 = np.loadtxt("./testData/body_1_test/fAddedMass.txt")
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fAddedMass'], result8))
        file = "./testData/body_1_test/irkb.mat"
        matFile = sio.loadmat(file) 
        keys = list(matFile.keys())[-1]
        result9 = matFile[keys]
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['irkb'], result9))

        
    def test_regExcitation(self):
        w = 0.785398163397448
        waveDir = [0]
        rho = 1000
        g = 9.81
        self.body_1.bodyNumber = 1
        self.body_1.readH5file()
        self.body_1.regExcitation(w,waveDir,rho,g)
        result1 = [18508.7581918053,0.223532780361825,1444829.62495620,6.41939080869329,140467.443724649,0.0412557412840770]
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fExt']['re'], result1))
        result2 = [546623.415011964,-0.0428753721889769,447276.984830190,0.273118282507816,4146814.04776144,0.0138224113420678]
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fExt']['im'], result2))
        result3 = [0,0,0,0,0,0]
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fExt']['md'], result3))
        
    # def test_constAddedMassAndDamping(self):
    #     """
    #     tested with Regulr CIC data so i need to check this again

    #     """
    #     w = 0.785398163397448
    #     rho = 1000
    #     B2B = 0
    #     CIkt = 601
    #     self.body_1.bodyNumber = 1
    #     self.body_1.readH5file()
    #     self.body_1.constAddedMassAndDamping(w,CIkt,rho,B2B)
    #     result1 =  np.loadtxt("./testData/body_1_test/fAddedMass.txt")
    #     self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fAddedMass'], result1))
    #     result2 =  np.loadtxt("./testData/body_1_test/fDamping.txt")
    #     self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fDamping'], result2))
        
    def test_irfInfAddedMassAndDamping(self):
        CTTime =  np.array(np.loadtxt("./testData/body_1_test/CTTime.txt"))
        ssCalc = 0
        CIkt = 601
        B2B = 0
        rho = 1000
        self.body_1.bodyNumber = 1
        self.body_1.readH5file()
        self.body_1.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        result1 = np.loadtxt("./testData/body_1_test/fAddedMass.txt")
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fAddedMass'], result1))
        file = "./testData/body_1_test/irkb.mat"
        matFile = sio.loadmat(file) 
        keys = list(matFile.keys())[-1]
        result2 = matFile[keys]
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['irkb'], result2))
        
if __name__ == '__main__':
    unittest.main()