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
def readData(file):
    matFile = sio.loadmat(file) 
    keys = list(matFile.keys())[-1]
    datas = matFile[keys]
    return datas

class TestBody(unittest.TestCase):
    
    def setUpClass():
        print("setupClass")
        
    def tearDownClass():
        print("teardownClass")
    
    def setUp(self):
        print("setUp")
        self.body_1 = BodyClass('rm3.h5') #regularCIC
        self.body_2 = BodyClass('rm3.h5') #irregular
        
    def tearDown(self):
        print("tearDown\n")
    
    # def readData(file):
    #     matFile = sio.loadmat(file) 
    #     keys = list(matFile.keys())[-1]
    #     datas = matFile[keys]
    #     return datas
    
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
        result9 = readData("./testData/body_1_test/irkb.mat")
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['irkb'], result9))
                
        w = np.conj(np.transpose(np.loadtxt("./testData/body_2_test/w.txt"))) 
        waveDir = [0]
        CIkt = 201
        CTTime =  np.array(np.loadtxt("./testData/body_2_test/CTTime.txt"))
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
        numFreq = 500
        self.body_2.bodyNumber = 1
        self.body_2.readH5file()
        self.body_2.hydroStiffness = np.zeros((6, 6))
        self.body_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_2.linearDamping = np.zeros((6, 6))
        self.body_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = np.loadtxt("./testData/body_2_test/linearHydroRestCoef.txt")
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['linearHydroRestCoef'], result1))
        result2 = np.zeros((6,6))
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['visDrag'], result2))
        result3 = np.zeros((6,6))
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['linearDamping'], result3))
        result4 = np.loadtxt("./testData/body_2_test/userDefinedFe.txt")
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['userDefinedFe'], result4)) 
        result5 = readData("./testData/body_2_test/re.mat")
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['re'], result5))
        result6 = readData("./testData/body_2_test/im.mat")
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['im'], result6))
        result7 = readData("./testData/body_2_test/md.mat")
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['md'], result7))
        result8 = np.loadtxt("./testData/body_2_test/fAddedMass.txt")
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fAddedMass'], result8))
        result9 = readData("./testData/body_2_test/irkb.mat")
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['irkb'], result9))
        
        
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
        
    def test_irrExcitation(self):
        w = np.conj(np.transpose(np.loadtxt("./testData/body_2_test/w.txt"))) 
        numFreq = 500
        waveDir = [0]
        rho = 1000
        g = 9.81
        self.body_2.bodyNumber = 1
        self.body_2.readH5file()
        self.body_2.irrExcitation(w,numFreq,waveDir,rho,g)
        result1 = readData("./testData/body_2_test/re.mat")
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['re'], result1))
        result2 = readData("./testData/body_2_test/im.mat")
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['im'], result2))
        result3 = readData("./testData/body_2_test/md.mat")
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['md'], result3))
        
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
        result2 = readData("./testData/body_1_test/irkb.mat") 
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['irkb'], result2))
        
        CTTime =  np.array(np.loadtxt("./testData/body_2_test/CTTime.txt"))
        ssCalc = 0
        CIkt = 201
        B2B = 0
        rho = 1000
        self.body_2.bodyNumber = 1
        self.body_2.readH5file()
        self.body_2.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        result1 = np.loadtxt("./testData/body_2_test/fAddedMass.txt")
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fAddedMass'], result1))
        result2 = readData("./testData/body_2_test/irkb.mat") 
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['irkb'], result2))
       
if __name__ == '__main__':
    unittest.main()