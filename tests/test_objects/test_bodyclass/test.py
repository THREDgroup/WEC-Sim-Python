#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Tue May 19 08:43:15 2020

@author: logical
'''

import unittest
import numpy as np
import scipy.io as sio
import os

from bodyclass import BodyClass

global cwd # set current directory as global variabl
cwd = os.getcwd()
lastDir = os.path.split(cwd)[1]
if lastDir == 'WEC-Sim-Python': # check if we are in desired directory
    cwd = cwd + '/tests/test_objects/test_bodyclass'

def readData(file): #read MATLAB file
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
        self.body_1 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regularCIC
        self.body_2 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #irregular
        self.body_3 = BodyClass(cwd + '/testData/hydroData/oswec.h5') #irregular wave dir = [0,30,90]
        self.body_4_1 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regular B2B_Case1:b2b = 0
        self.body_4_2 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regular B2B_Case1:b2b = 0
        self.body_5_1 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regular B2B_Case2:b2b = 1
        self.body_5_2 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regular B2B_Case2:b2b = 1
        self.body_6_1 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regularCIC B2B_Case3:b2b = 0
        self.body_6_2 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regularCIC B2B_Case3:b2b = 0
        self.body_7_1 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regularCIC B2B_Case4:b2b = 1
        self.body_7_2 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regularCIC B2B_Case4:b2b = 1
        self.body_8_1 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regularCIC B2B_Case5:b2b = 0,ssCalc = 1
        self.body_8_2 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regularCIC B2B_Case5:b2b = 0,ssCalc = 1
        self.body_9_1 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regularcIC B2B_Case6:b2b = 1,ssCalc = 1
        self.body_9_2 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regularcIC B2B_Case6:b2b = 1,ssCalc = 1
        self.body_10_1 = BodyClass(cwd + '/testData/hydroData/ellipsoid.h5') #regular nlHydro = 2
        
    def tearDown(self):
        print("tearDown\n")
    
    def test_hydroForcePre(self):
        # RM3 example from WEC-Sim
        # regularCIC
        w = 0.785398163397448
        waveDir = [0]
        CIkt = 601
        CTTime = np.conj(np.transpose(np.loadtxt(cwd + '/testData/body_1_test/CTTime.txt'))) 
        numFreq = []
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'regularCIC'
        waveAmpTime = np.conj(np.transpose(np.loadtxt(cwd + '/testData/body_1_test/waveAmpTime.txt'))) 
        iBod = 1 #body number
        numBod = [] 
        ssCalc = 0
        nlHydro = 0
        B2B = 0
        self.body_1.bodyNumber = 1
        self.body_1.readH5file()
        self.body_1.hydroStiffness = np.zeros((6, 6))
        self.body_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_1.linearDamping = np.zeros((6, 6))
        self.body_1.mass = 'equilibrium'
        self.body_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = np.loadtxt(cwd + '/testData/body_1_test/linearHydroRestCoef.txt')
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['linearHydroRestCoef'], result1))
        result2 = np.zeros((6,6))
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['visDrag'], result2))
        result3 = np.zeros((6,6))
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['linearDamping'], result3))
        result4 = np.loadtxt(cwd + '/testData/body_1_test/userDefinedFe.txt')
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['userDefinedFe'], result4)) 
        result5 = [18508.7581918053,0.223532780361825,1444829.62495620,6.41939080869329,140467.443724649,0.0412557412840770]
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fExt']['re'], result5))
        result6 = [546623.415011964,-0.0428753721889769,447276.984830190,0.273118282507816,4146814.04776144,0.0138224113420678]
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fExt']['im'], result6))
        result7 = [0,0,0,0,0,0]
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fExt']['md'], result7))
        result8 = np.loadtxt(cwd + '/testData/body_1_test/fAddedMass.txt')
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_1_test/irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['irkb'], result9))
        
        #rm3 example with irregular variation    
        # irregular
        w = np.conj(np.transpose(np.loadtxt(cwd + '/testData/body_2_test/w.txt'))) 
        waveDir = [0]
        CIkt = 201
        CTTime =  np.array(np.loadtxt(cwd + '/testData/body_2_test/CTTime.txt'))
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'irregular'
        waveAmpTime = np.conj(np.transpose(np.loadtxt(cwd + '/testData/body_2_test/waveAmpTime.txt'))) 
        iBod = 1 #body number
        numBod = [] 
        ssCalc = 0
        nlHydro = 0
        B2B = 0
        numFreq = 500
        self.body_2.bodyNumber = 1
        self.body_2.readH5file()
        self.body_2.hydroStiffness = np.zeros((6, 6))
        self.body_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_2.linearDamping = np.zeros((6, 6))
        self.body_2.mass = 'equilibrium'
        self.body_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = np.loadtxt(cwd + '/testData/body_2_test/linearHydroRestCoef.txt')
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['linearHydroRestCoef'], result1))
        result2 = np.zeros((6,6))
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['visDrag'], result2))
        result3 = np.zeros((6,6))
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['linearDamping'], result3))
        result4 = np.loadtxt(cwd + '/testData/body_2_test/userDefinedFe.txt')
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_2_test/re.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_2_test/im.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_2_test/md.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['md'], result7))
        result8 = np.loadtxt(cwd + '/testData/body_2_test/fAddedMass.txt')
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_2_test/irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['irkb'], result9))
        self.assertIsNone(np.testing.assert_allclose(self.body_2.mass, 725834))
        
        # OSWEC example from WEC-Sim
        # irregular wave dir = [0,30,90]
        w = np.conj(np.transpose(readData(cwd + '/testData/body_3_test/w.mat'))) 
        waveDir = [0,30,90]
        CIkt = 301
        CTTime =  readData(cwd + '/testData/body_3_test/CTTime.mat')[0]
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'irregular'
        waveAmpTime = np.conj(np.transpose(readData(cwd + '/testData/body_3_test/waveAmpTime.mat'))) 
        iBod = 1 #body number
        numBod = [] 
        ssCalc = 0
        nlHydro = 0
        B2B = 0
        numFreq = 500
        self.body_3.bodyNumber = 1
        self.body_3.readH5file()
        self.body_3.hydroStiffness = np.zeros((6, 6))
        self.body_3.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_3.linearDamping = np.zeros((6, 6))
        self.body_3.mass = 127000
        self.body_3.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_3_test/linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['linearHydroRestCoef'], result1))
        result2 = np.zeros((6,6))
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['visDrag'], result2))
        result3 = np.zeros((6,6))
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_3_test/userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_3_test/re.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_3_test/im.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_3_test/md.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_3_test/fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_3_test/irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['irkb'], result9))
        self.assertIsNone(np.testing.assert_allclose(self.body_3.mass, 127000))
        
        # WEC-Sim-Application B2B 
        # regular B2B_Case1:b2b = 0
        w = np.conj(np.transpose(readData(cwd + '/testData/body_4_test/w.mat'))) 
        waveDir = [0]
        CIkt = 601
        CTTime = readData(cwd + '/testData/body_4_test/CTTime.mat')[0]
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'regular'
        waveAmpTime = np.conj(np.transpose(readData(cwd + '/testData/body_4_test/waveAmpTime.mat'))) 
        iBod = 1 #body number
        numBod = [] 
        ssCalc = 0
        nlHydro = 0
        B2B = 0
        numFreq = []
        self.body_4_1.bodyNumber = 1
        self.body_4_1.bodyTotal = [2]
        self.body_4_1.readH5file()
        self.body_4_1.hydroStiffness = np.zeros((6, 6))
        self.body_4_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_4_1.linearDamping = np.zeros((6, 6))
        self.body_4_1.mass = 'equilibrium'
        self.body_4_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_4_test/body1_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_4_test/body1_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_4_test/body1_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_4_test/body1_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_4_test/body1_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_4_test/body1_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_4_test/body1_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_4_test/body1_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_4_test/body1_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_4_test/body1_totDOF.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['totDOF'], result10))
        iBod = 2 #body number
        self.body_4_2.bodyNumber = 2
        self.body_4_2.bodyTotal = [2]
        self.body_4_2.readH5file()
        self.body_4_2.hydroStiffness = np.zeros((6, 6))
        self.body_4_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_4_2.linearDamping = np.zeros((6, 6))
        self.body_4_2.mass = 'equilibrium'
        self.body_4_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_4_test/body2_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_4_test/body2_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_4_test/body2_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_4_test/body2_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_4_test/body2_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_4_test/body2_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_4_test/body2_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_4_test/body2_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_4_test/body2_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_4_test/body2_totDOF.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['totDOF'], result10))
        
        # WEC-Sim-Application B2B 
        # regular B2B_Case1:b2b = 1
        w = np.conj(np.transpose(readData(cwd + '/testData/body_5_test/w.mat'))) 
        waveDir = [0]
        CIkt = 601
        CTTime = readData(cwd + '/testData/body_5_test/CTTime.mat')[0]
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'regular'
        waveAmpTime = np.conj(np.transpose(readData(cwd + '/testData/body_5_test/waveAmpTime.mat'))) 
        iBod = 1 #body number
        numBod = [] 
        ssCalc = 0
        nlHydro = 0
        B2B = 1
        numFreq = []
        self.body_5_1.bodyNumber = 1
        self.body_5_1.bodyTotal = [2]
        self.body_5_1.readH5file()
        self.body_5_1.hydroStiffness = np.zeros((6, 6))
        self.body_5_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_5_1.linearDamping = np.zeros((6, 6))
        self.body_5_1.mass = 'equilibrium'
        self.body_5_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_5_test/body1_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_5_test/body1_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_5_test/body1_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_5_test/body1_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_5_test/body1_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_5_test/body1_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_5_test/body1_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_5_test/body1_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_5_test/body1_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_5_test/body1_totDOF.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['totDOF'], result10))
        iBod = 2 #body number
        self.body_5_2.bodyNumber = 2
        self.body_5_2.bodyTotal = [2]
        self.body_5_2.readH5file()
        self.body_5_2.hydroStiffness = np.zeros((6, 6))
        self.body_5_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_5_2.linearDamping = np.zeros((6, 6))
        self.body_5_2.mass = 'equilibrium'
        self.body_5_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_5_test/body2_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_5_test/body2_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_5_test/body2_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_5_test/body2_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_5_test/body2_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_5_test/body2_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_5_test/body2_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_5_test/body2_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_5_test/body2_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_5_test/body2_totDOF.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['totDOF'], result10))
        
        # WEC-Sim-Application B2B 
        # regularCIC B2B_Case1:b2b = 0
        w = np.conj(np.transpose(readData(cwd + '/testData/body_6_test/w.mat'))) 
        waveDir = [0]
        CIkt = 601
        CTTime = readData(cwd + '/testData/body_6_test/CTTime.mat')[0]
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'regularCIC'
        waveAmpTime = np.conj(np.transpose(readData(cwd + '/testData/body_6_test/waveAmpTime.mat'))) 
        iBod = 1 #body number
        numBod = [] 
        ssCalc = 0
        nlHydro = 0
        B2B = 0
        numFreq = []
        self.body_6_1.bodyNumber = 1
        self.body_6_1.bodyTotal = [2]
        self.body_6_1.readH5file()
        self.body_6_1.hydroStiffness = np.zeros((6, 6))
        self.body_6_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_6_1.linearDamping = np.zeros((6, 6))
        self.body_6_1.mass = 'equilibrium'
        self.body_6_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_6_test/body1_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_1.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_6_test/body1_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_1.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_6_test/body1_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_1.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_6_test/body1_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_1.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_6_test/body1_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_6_1.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_6_test/body1_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_6_1.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_6_test/body1_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_6_1.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_6_test/body1_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_1.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_6_test/body1_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_1.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_6_test/body1_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_1.hydroForce['irkb'], result10))
        iBod = 2 #body number
        self.body_6_2.bodyNumber = 2
        self.body_6_2.bodyTotal = [2]
        self.body_6_2.readH5file()
        self.body_6_2.hydroStiffness = np.zeros((6, 6))
        self.body_6_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_6_2.linearDamping = np.zeros((6, 6))
        self.body_6_2.mass = 'equilibrium'
        self.body_6_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_6_test/body2_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_2.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_6_test/body2_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_2.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_6_test/body2_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_2.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_6_test/body2_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_2.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_6_test/body2_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_6_2.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_6_test/body2_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_6_2.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_6_test/body2_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_6_2.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_6_test/body2_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_2.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_6_test/body2_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_2.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_6_test/body2_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_6_2.hydroForce['irkb'], result10))
        
        # WEC-Sim-Application B2B 
        # regularCIC B2B_Case1:b2b = 1
        w = np.conj(np.transpose(readData(cwd + '/testData/body_7_test/w.mat'))) 
        waveDir = [0]
        CIkt = 601
        CTTime = readData(cwd + '/testData/body_7_test/CTTime.mat')[0]
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'regularCIC'
        waveAmpTime = np.conj(np.transpose(readData(cwd + '/testData/body_7_test/waveAmpTime.mat'))) 
        iBod = 1 #body number
        numBod = [] 
        ssCalc = 0
        nlHydro = 0
        B2B = 1
        numFreq = []
        self.body_7_1.bodyNumber = 1
        self.body_7_1.bodyTotal = [2]
        self.body_7_1.readH5file()
        self.body_7_1.hydroStiffness = np.zeros((6, 6))
        self.body_7_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_7_1.linearDamping = np.zeros((6, 6))
        self.body_7_1.mass = 'equilibrium'
        self.body_7_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_7_test/body1_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_7_test/body1_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_7_test/body1_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_7_test/body1_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_7_test/body1_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_7_test/body1_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_7_test/body1_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_7_test/body1_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_7_test/body1_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_7_test/body1_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['irkb'], result10))
        iBod = 2 #body number
        self.body_7_2.bodyNumber = 2
        self.body_7_2.bodyTotal = [2]
        self.body_7_2.readH5file()
        self.body_7_2.hydroStiffness = np.zeros((6, 6))
        self.body_7_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_7_2.linearDamping = np.zeros((6, 6))
        self.body_7_2.mass = 'equilibrium'
        self.body_7_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_7_test/body2_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_7_test/body2_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_7_test/body2_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_7_test/body2_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_7_test/body2_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_7_test/body2_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_7_test/body2_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_7_test/body2_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_7_test/body2_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_7_test/body2_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['irkb'], result10))
        
        # WEC-Sim-Application B2B 
        # regularCIC B2B_Case1:b2b = 0,ssCalc =1
        w = np.conj(np.transpose(readData(cwd + '/testData/body_8_test/w.mat'))) 
        waveDir = [0]
        CIkt = 601
        CTTime = readData(cwd + '/testData/body_8_test/CTTime.mat')[0]
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'regularCIC'
        waveAmpTime = np.conj(np.transpose(readData(cwd + '/testData/body_8_test/waveAmpTime.mat'))) 
        iBod = 1 #body number
        numBod = []
        ssCalc = 1
        nlHydro = 0
        B2B = 0
        numFreq = []
        self.body_8_1.bodyNumber = 1
        self.body_8_1.bodyTotal = [2]
        self.body_8_1.readH5file()
        self.body_8_1.hydroStiffness = np.zeros((6, 6))
        self.body_8_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_8_1.linearDamping = np.zeros((6, 6))
        self.body_8_1.mass = 'equilibrium'
        self.body_8_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_8_test/body1_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_8_test/body1_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_8_test/body1_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_8_test/body1_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_8_test/body1_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_8_test/body1_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_8_test/body1_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_8_test/body1_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_8_test/body1_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_8_test/body1_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['irkb'], result10))
        result11 = readData(cwd + '/testData/body_8_test/body1_A.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['ssRadf']['A'], result11))
        result12 = readData(cwd + '/testData/body_8_test/body1_B.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['ssRadf']['B'], result12))
        result13 = readData(cwd + '/testData/body_8_test/body1_C.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['ssRadf']['C'], result13))
        result14 = readData(cwd + '/testData/body_8_test/body1_D.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['ssRadf']['D'], result14))
        iBod = 2 #body number
        self.body_8_2.bodyNumber = 2
        self.body_8_2.bodyTotal = [2]
        self.body_8_2.readH5file()
        self.body_8_2.hydroStiffness = np.zeros((6, 6))
        self.body_8_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_8_2.linearDamping = np.zeros((6, 6))
        self.body_8_2.mass = 'equilibrium'
        self.body_8_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_8_test/body2_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_8_test/body2_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_8_test/body2_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_8_test/body2_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_8_test/body2_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_8_test/body2_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_8_test/body2_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_8_test/body2_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_8_test/body2_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_8_test/body2_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['irkb'], result10))
        result11 = readData(cwd + '/testData/body_8_test/body2_A.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['ssRadf']['A'], result11))
        result12 = readData(cwd + '/testData/body_8_test/body2_B.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['ssRadf']['B'], result12))
        result13 = readData(cwd + '/testData/body_8_test/body2_C.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['ssRadf']['C'], result13))
        result14 = readData(cwd + '/testData/body_8_test/body2_D.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['ssRadf']['D'], result14))
        
        # WEC-Sim-Application B2B 
        # regularCIC B2B_Case1:b2b = 1,ssCalc =1
        w = np.conj(np.transpose(readData(cwd + '/testData/body_9_test/w.mat'))) 
        waveDir = [0]
        CIkt = 601
        CTTime = readData(cwd + '/testData/body_9_test/CTTime.mat')[0]
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'regularCIC'
        waveAmpTime = np.conj(np.transpose(readData(cwd + '/testData/body_9_test/waveAmpTime.mat'))) 
        iBod = 1 #body number
        numBod = []
        ssCalc = 1
        nlHydro = 0
        B2B = 1
        numFreq = []
        self.body_9_1.bodyNumber = 1
        self.body_9_1.bodyTotal = [2]
        self.body_9_1.readH5file()
        self.body_9_1.hydroStiffness = np.zeros((6, 6))
        self.body_9_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_9_1.linearDamping = np.zeros((6, 6))
        self.body_9_1.mass = 'equilibrium'
        self.body_9_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_9_test/body1_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_9_test/body1_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_9_test/body1_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_9_test/body1_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_9_test/body1_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_9_test/body1_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_9_test/body1_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_9_test/body1_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_9_test/body1_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_9_test/body1_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['irkb'], result10))
        result11 = readData(cwd + '/testData/body_9_test/body1_A.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['ssRadf']['A'], result11))
        result12 = readData(cwd + '/testData/body_9_test/body1_B.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['ssRadf']['B'], result12))
        result13 = readData(cwd + '/testData/body_9_test/body1_C.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['ssRadf']['C'], result13))
        result14 = readData(cwd + '/testData/body_9_test/body1_D.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['ssRadf']['D'], result14))
        iBod = 2 #body number
        self.body_9_2.bodyNumber = 2
        self.body_9_2.bodyTotal = [2]
        self.body_9_2.readH5file()
        self.body_9_2.hydroStiffness = np.zeros((6, 6))
        self.body_9_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_9_2.linearDamping = np.zeros((6, 6))
        self.body_9_2.mass = 'equilibrium'
        self.body_9_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_9_test/body2_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_9_test/body2_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_9_test/body2_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_9_test/body2_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_9_test/body2_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_9_test/body2_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_9_test/body2_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_9_test/body2_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_9_test/body2_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_9_test/body2_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['irkb'], result10))
        result11 = readData(cwd + '/testData/body_9_test/body2_A.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['ssRadf']['A'], result11))
        result12 = readData(cwd + '/testData/body_9_test/body2_B.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['ssRadf']['B'], result12))
        result13 = readData(cwd + '/testData/body_9_test/body2_C.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['ssRadf']['C'], result13))
        result14 = readData(cwd + '/testData/body_9_test/body2_D.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['ssRadf']['D'], result14))
        
        # WEC-Sim-Application Non-linear Hydro
        # regular B2B_Case1:nlHydro = 2
        w = np.conj(np.transpose(readData(cwd + '/testData/body_10_test/w.mat'))) 
        waveDir = [0]
        CIkt = 1201
        CTTime = readData(cwd + '/testData/body_10_test/CTTime.mat')[0]
        dt = 0.05
        rho = 1025
        g = 9.81
        waveType = 'regular'
        waveAmpTime = np.conj(np.transpose(readData(cwd + '/testData/body_10_test/waveAmpTime.mat'))) 
        iBod = 1 #body number
        numBod = [] 
        ssCalc = 0
        nlHydro = 2
        B2B = 0
        numFreq = []
        self.body_10_1.bodyNumber = 1
        self.body_10_1.bodyTotal = [1]
        self.body_10_1.readH5file()
        self.body_10_1.hydroStiffness = np.zeros((6, 6))
        self.body_10_1.viscDrag = {'Drag':np.zeros((6, 6)),
                                    'cd':np.array([1,0,1,0,1,0]),
                                    'characteristicArea':np.array([25,0,np.pi*5**2,0,np.pi*5**5,0])}
        self.body_10_1.linearDamping = np.zeros((6, 6))
        self.body_10_1.mass = 'equilibrium'
        self.body_10_1.bodyGeo(cwd + '/testData/body_10_test/elipsoid.stl')
        self.body_10_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        result1 = readData(cwd + '/testData/body_10_test/body10_linearHydroRestCoef.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.hydroForce['linearHydroRestCoef'], result1))
        result2 = readData(cwd + '/testData/body_10_test/body10_visDrag.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.hydroForce['visDrag'], result2))
        result3 = readData(cwd + '/testData/body_10_test/body10_linearDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.hydroForce['linearDamping'], result3))
        result4 = readData(cwd + '/testData/body_10_test/body10_userDefinedFe.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.hydroForce['userDefinedFe'], result4)) 
        result5 = readData(cwd + '/testData/body_10_test/body10_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.hydroForce['fExt']['re'], result5))
        result6 = readData(cwd + '/testData/body_10_test/body10_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.hydroForce['fExt']['im'], result6))
        result7 = readData(cwd + '/testData/body_10_test/body10_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.hydroForce['fExt']['md'], result7))
        result8 = readData(cwd + '/testData/body_10_test/body10_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.hydroForce['fAddedMass'], result8))
        result9 = readData(cwd + '/testData/body_10_test/body10_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.hydroForce['fDamping'], result9))
        result10 = readData(cwd + '/testData/body_10_test/body10_vertex.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.bodyGeometry['vertex'], result10))
        result11 = readData(cwd + '/testData/body_10_test/body10_face.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.bodyGeometry['face'], result11))
        result12 = readData(cwd + '/testData/body_10_test/body10_norm.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.bodyGeometry['norm'], result12))
        result13 = readData(cwd + '/testData/body_10_test/body10_area.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.bodyGeometry['area'], result13))
        result14 = readData(cwd + '/testData/body_10_test/body10_center.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.bodyGeometry['center'], result14))
        self.assertEqual(self.body_10_1.bodyGeometry['numVertex'], 1442)
        self.assertEqual(self.body_10_1.bodyGeometry['numFace'], 2880)
        
    def test_regExcitation(self):
        '''
        test regular excitation

        '''
        # random number test
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
        
        # regular wave case
        w = np.conj(np.transpose(readData(cwd + '/testData/body_4_test/w.mat'))) 
        waveDir = [0]
        rho = 1000
        g = 9.81
        self.body_4_1.bodyNumber = 1
        self.body_4_1.bodyTotal = 2
        self.body_4_1.readH5file()
        self.body_4_1.regExcitation(w,waveDir,rho,g)
        result1 = readData(cwd + '/testData/body_4_test/body1_re.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['fExt']['re'], result1))
        result2 = readData(cwd + '/testData/body_4_test/body1_im.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['fExt']['im'], result2))
        result3 = readData(cwd + '/testData/body_4_test/body1_md.mat')[0]
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['fExt']['md'], result3))
        
    def test_irrExcitation(self):
        '''
        test irregular excitation

        '''
        # single wavDir case
        w = np.conj(np.transpose(np.loadtxt(cwd + '/testData/body_2_test/w.txt'))) 
        numFreq = 500
        waveDir = [0]
        rho = 1000
        g = 9.81
        self.body_2.bodyNumber = 1
        self.body_2.readH5file()
        self.body_2.irrExcitation(w,numFreq,waveDir,rho,g)
        result1 = readData(cwd + '/testData/body_2_test/re.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['re'], result1))
        result2 = readData(cwd + '/testData/body_2_test/im.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['im'], result2))
        result3 = readData(cwd + '/testData/body_2_test/md.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fExt']['md'], result3))
        
        #multiple waveDir case
        w = np.conj(np.transpose(readData(cwd + '/testData/body_3_test/w.mat'))) 
        numFreq = 500
        waveDir = [0,30,90]
        rho = 1000
        g = 9.81
        self.body_3.bodyNumber = 1
        self.body_3.readH5file()
        self.body_3.irrExcitation(w,numFreq,waveDir,rho,g)
        result1 = readData(cwd + '/testData/body_3_test/re.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['fExt']['re'], result1))
        result2 = readData(cwd + '/testData/body_3_test/im.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['fExt']['im'], result2))
        result3 = readData(cwd + '/testData/body_3_test/md.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['fExt']['md'], result3))
        
    def test_constAddedMassAndDamping(self):
        '''
        test constant added mass and damping

        '''
        # B2B = 0 case
        # Regular
        w = np.conj(np.transpose(readData(cwd + '/testData/body_4_test/w.mat'))) 
        rho = 1000
        B2B = 0
        CIkt = 601
        self.body_4_1.bodyNumber = 1
        self.body_4_1.bodyTotal = [2]
        self.body_4_1.readH5file()
        self.body_4_1.constAddedMassAndDamping(w,CIkt,rho,B2B)
        result1 = readData(cwd + '/testData/body_4_test/body1_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['fAddedMass'], result1))
        result2 = readData(cwd + '/testData/body_4_test/body1_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_1.hydroForce['fDamping'], result2))
        self.body_4_2.bodyNumber = 2
        self.body_4_2.bodyTotal = [2]
        self.body_4_2.readH5file()
        self.body_4_2.constAddedMassAndDamping(w,CIkt,rho,B2B)
        result3 = readData(cwd + '/testData/body_4_test/body2_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['fAddedMass'], result3))
        result4 = readData(cwd + '/testData/body_4_test/body2_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_4_2.hydroForce['fDamping'], result4))
        
        # B2B = 1 case
        # Regular
        w = np.conj(np.transpose(readData(cwd + '/testData/body_5_test/w.mat'))) 
        rho = 1000
        B2B = 1
        CIkt = 601
        self.body_5_1.bodyNumber = 1
        self.body_5_1.bodyTotal = [2]
        self.body_5_1.readH5file()
        self.body_5_1.constAddedMassAndDamping(w,CIkt,rho,B2B)
        result1 = readData(cwd + '/testData/body_5_test/body1_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['fAddedMass'], result1))
        result2 = readData(cwd + '/testData/body_5_test/body1_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_1.hydroForce['fDamping'], result2))
        self.body_5_2.bodyNumber = 2
        self.body_5_2.bodyTotal = [2]
        self.body_5_2.readH5file()
        self.body_5_2.constAddedMassAndDamping(w,CIkt,rho,B2B)
        result3 = readData(cwd + '/testData/body_5_test/body2_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['fAddedMass'], result3))
        result4 = readData(cwd + '/testData/body_5_test/body2_fDamping.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_5_2.hydroForce['fDamping'], result4))
        
    def test_irfInfAddedMassAndDamping(self):
        '''
        test irfInfAddedMassAndDamping
        '''
        # ssCalc = 0 case
        # B2B = 0 case
        # Regular CIC
        CTTime =  np.array(np.loadtxt(cwd + '/testData/body_1_test/CTTime.txt'))
        ssCalc = 0
        CIkt = 601
        B2B = 0
        rho = 1000
        self.body_1.bodyNumber = 1
        self.body_1.readH5file()
        self.body_1.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        result1 = np.loadtxt(cwd + '/testData/body_1_test/fAddedMass.txt')
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['fAddedMass'], result1))
        result2 = readData(cwd + '/testData/body_1_test/irkb.mat') 
        self.assertIsNone(np.testing.assert_allclose(self.body_1.hydroForce['irkb'], result2))
        
        # ssCalc = 0 case
        # B2B = 0 case
        # Irregular
        CTTime =  np.array(np.loadtxt(cwd + '/testData/body_2_test/CTTime.txt'))
        ssCalc = 0
        CIkt = 201
        B2B = 0
        rho = 1000
        self.body_2.bodyNumber = 1
        self.body_2.readH5file()
        self.body_2.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        result1 = np.loadtxt(cwd + '/testData/body_2_test/fAddedMass.txt')
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['fAddedMass'], result1))
        result2 = readData(cwd + '/testData/body_2_test/irkb.mat') 
        self.assertIsNone(np.testing.assert_allclose(self.body_2.hydroForce['irkb'], result2))
        
        # ssCalc = 0 case
        # B2B = 0 case
        # Irregular
        CTTime =  readData(cwd + '/testData/body_3_test/CTTime.mat')[0]
        ssCalc = 0
        CIkt = 301
        B2B = 0
        rho = 1000
        self.body_3.bodyNumber = 1
        self.body_3.readH5file()
        self.body_3.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        result1 = readData(cwd + '/testData/body_3_test/fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['fAddedMass'], result1))
        result2 = readData(cwd + '/testData/body_3_test/irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_3.hydroForce['irkb'], result2))
        
        # ssCalc = 0 case
        # B2B = 1 case
        # Regular CIC
        CTTime =  readData(cwd + '/testData/body_7_test/CTTime.mat')[0]
        ssCalc = 0
        CIkt = 601
        B2B = 1
        rho = 1000
        self.body_7_1.bodyNumber = 1
        self.body_7_1.bodyTotal = [2]
        self.body_7_1.readH5file()
        self.body_7_1.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        result1 = readData(cwd + '/testData/body_7_test/body1_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['fAddedMass'], result1))
        result2 = readData(cwd + '/testData/body_7_test/body1_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_1.hydroForce['irkb'], result2))
        self.body_7_2.bodyNumber = 2
        self.body_7_2.bodyTotal = [2]
        self.body_7_2.readH5file()
        self.body_7_2.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        result1 = readData(cwd + '/testData/body_7_test/body2_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['fAddedMass'], result1))
        result2 = readData(cwd + '/testData/body_7_test/body2_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_7_2.hydroForce['irkb'], result2))
        
        # ssCalc = 1 case
        # B2B = 0 case
        # Regular CIC
        CTTime = readData(cwd + '/testData/body_8_test/CTTime.mat')[0]
        ssCalc = 1
        CIkt = 601
        B2B = 0
        rho = 1000
        self.body_8_1.bodyNumber = 1
        self.body_8_1.bodyTotal = [2]
        self.body_8_1.readH5file()
        self.body_8_1.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        result1 = readData(cwd + '/testData/body_8_test/body1_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['fAddedMass'], result1))
        result2 = readData(cwd + '/testData/body_8_test/body1_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['irkb'], result2))
        result3 = readData(cwd + '/testData/body_8_test/body1_A.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['ssRadf']['A'], result3))
        result4 = readData(cwd + '/testData/body_8_test/body1_B.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['ssRadf']['B'], result4))
        result5 = readData(cwd + '/testData/body_8_test/body1_C.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['ssRadf']['C'], result5))
        result6 = readData(cwd + '/testData/body_8_test/body1_D.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_1.hydroForce['ssRadf']['D'], result6))
        self.body_8_2.bodyNumber = 2
        self.body_8_2.bodyTotal = [2]
        self.body_8_2.readH5file()
        self.body_8_2.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        result1 = readData(cwd + '/testData/body_8_test/body2_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['fAddedMass'], result1))
        result2 = readData(cwd + '/testData/body_8_test/body2_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['irkb'], result2))
        result3 = readData(cwd + '/testData/body_8_test/body2_A.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['ssRadf']['A'], result3))
        result4 = readData(cwd + '/testData/body_8_test/body2_B.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['ssRadf']['B'], result4))
        result5 = readData(cwd + '/testData/body_8_test/body2_C.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['ssRadf']['C'], result5))
        result6 = readData(cwd + '/testData/body_8_test/body2_D.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_8_2.hydroForce['ssRadf']['D'], result6))
        
        # ssCalc = 1 case
        # B2B = 1 case
        # Regular CIC
        CTTime = readData(cwd + '/testData/body_9_test/CTTime.mat')[0]
        ssCalc = 1
        CIkt = 601
        B2B = 1
        rho = 1000
        self.body_9_1.bodyNumber = 1
        self.body_9_1.bodyTotal = [2]
        self.body_9_1.readH5file()
        self.body_9_1.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        result1 = readData(cwd + '/testData/body_9_test/body1_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['fAddedMass'], result1))
        result2 = readData(cwd + '/testData/body_9_test/body1_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['irkb'], result2))
        result3 = readData(cwd + '/testData/body_9_test/body1_A.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['ssRadf']['A'], result3))
        result4 = readData(cwd + '/testData/body_9_test/body1_B.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['ssRadf']['B'], result4))
        result5 = readData(cwd + '/testData/body_9_test/body1_C.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['ssRadf']['C'], result5))
        result6 = readData(cwd + '/testData/body_9_test/body1_D.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_1.hydroForce['ssRadf']['D'], result6))
        self.body_9_2.bodyNumber = 2
        self.body_9_2.bodyTotal = [2]
        self.body_9_2.readH5file()
        self.body_9_2.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        result1 = readData(cwd + '/testData/body_9_test/body2_fAddedMass.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['fAddedMass'], result1))
        result2 = readData(cwd + '/testData/body_9_test/body2_irkb.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['irkb'], result2))
        result3 = readData(cwd + '/testData/body_9_test/body2_A.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['ssRadf']['A'], result3))
        result4 = readData(cwd + '/testData/body_9_test/body2_B.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['ssRadf']['B'], result4))
        result5 = readData(cwd + '/testData/body_9_test/body2_C.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['ssRadf']['C'], result5))
        result6 = readData(cwd + '/testData/body_9_test/body2_D.mat')
        self.assertIsNone(np.testing.assert_allclose(self.body_9_2.hydroForce['ssRadf']['D'], result6))
       
    def test_setMassMatrix(self):
        '''
        test setMassMatrix

        '''
        # mass = 'equilibirum' case
        # nlHydro = 0
        rho = 1000
        nlHydro = 0
        self.body_2.bodyNumber = 1
        self.body_2.readH5file()
        self.body_2.mass = 'equilibrium'
        self.body_2.setMassMatrix(rho,nlHydro)
        result1 = 'equilibrium'
        self.assertEqual(self.body_2.massCalcMethod, result1)
        result2 = 725834
        self.assertIsNone(np.testing.assert_allclose(self.body_2.mass, result2))
        
        # mass is given integer case
        # nlHydro = 0
        rho = 1000
        nlHydro = 0
        self.body_2.bodyNumber = 1
        self.body_2.readH5file()
        self.body_2.mass = 127000
        self.body_2.setMassMatrix(rho,nlHydro)
        result1 = 'user'
        self.assertEqual(self.body_2.massCalcMethod, result1)
        result2 = 127000
        self.assertIsNone(np.testing.assert_allclose(self.body_2.mass, result2))
        
        # mass = 'equilibirum' case
        # nlHydro = 2
        rho = 1025
        nlHydro = 2
        self.body_10_1.bodyNumber = 1
        self.body_10_1.readH5file()
        self.body_10_1.mass = 'equilibrium'
        self.body_10_1.bodyGeo(cwd + '/testData/body_10_test/elipsoid.stl')
        self.body_10_1.setMassMatrix(rho,nlHydro)
        result1 = 133376.066729747
        self.assertIsNone(np.testing.assert_allclose(self.body_10_1.mass, result1))
        
    def test_adjustMassMatrix(self):
        '''
        test adjustMassMatrix

        '''
        # B2B = 0 case
        w = np.conj(np.transpose(readData(cwd + '/testData/body_4_test/w.mat'))) 
        waveDir = [0]
        CIkt = 601
        CTTime = readData(cwd + '/testData/body_4_test/CTTime.mat')[0]
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'regular'
        waveAmpTime = np.conj(np.transpose(readData(cwd + '/testData/body_4_test/waveAmpTime.mat'))) 
        iBod = 1 #body number
        numBod = [] #later change it to 2 to check
        ssCalc = 0
        nlHydro = 0
        B2B = 0
        numFreq = []
        self.body_4_1.bodyNumber = 1
        self.body_4_1.bodyTotal = [2]
        self.body_4_1.readH5file()
        self.body_4_1.hydroStiffness = np.zeros((6, 6))
        self.body_4_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_4_1.linearDamping = np.zeros((6, 6))
        self.body_4_1.mass = 'equilibrium'
        self.body_4_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        adjMassWeightFun = 5
        self.body_4_1.momOfInertia = [20907301, 21306090.66, 37085481.11];  
        self.body_4_1.adjustMassMatrix(adjMassWeightFun,B2B)
        result1 = readData(cwd + '/testData/body_4_test/body1_fAddedMass_adjusted.mat')
        np.testing.assert_allclose(self.body_4_1.hydroForce['fAddedMass'], result1)
        
        # B2B = 1 case
        # multiple body case
        w = np.conj(np.transpose(readData(cwd + '/testData/body_5_test/w.mat'))) 
        waveDir = [0]
        CIkt = 601
        CTTime = readData(cwd + '/testData/body_5_test/CTTime.mat')[0]
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'regular'
        waveAmpTime = np.conj(np.transpose(readData(cwd + '/testData/body_5_test/waveAmpTime.mat'))) 
        iBod = 1 #body number
        numBod = [] #later change it to 2 to check
        ssCalc = 0
        nlHydro = 0
        B2B = 1
        numFreq = []
        self.body_5_1.bodyNumber = 1
        self.body_5_1.bodyTotal = [2]
        self.body_5_1.readH5file()
        self.body_5_1.hydroStiffness = np.zeros((6, 6))
        self.body_5_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_5_1.linearDamping = np.zeros((6, 6))
        self.body_5_1.mass = 'equilibrium'
        self.body_5_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        adjMassWeightFun = 5
        self.body_5_1.momOfInertia = [20907301, 21306090.66, 37085481.11]  
        self.body_5_1.adjustMassMatrix(adjMassWeightFun,B2B)
        iBod = 2 #body number
        self.body_5_2.bodyNumber = 2
        self.body_5_2.bodyTotal = [2]
        self.body_5_2.readH5file()
        self.body_5_2.hydroStiffness = np.zeros((6, 6))
        self.body_5_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_5_2.linearDamping = np.zeros((6, 6))
        self.body_5_2.mass = 'equilibrium'
        self.body_5_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        self.body_5_2.momOfInertia = [94419614.57, 94407091.24, 28542224.82]
        self.body_5_2.adjustMassMatrix(adjMassWeightFun,B2B)
        result1 = readData(cwd + '/testData/body_5_test/body1_fAddedMass_adjusted.mat')
        np.testing.assert_allclose(self.body_5_1.hydroForce['fAddedMass'], result1)
        result2 = readData(cwd + '/testData/body_5_test/body2_fAddedMass_adjusted.mat')
        np.testing.assert_allclose(self.body_5_2.hydroForce['fAddedMass'], result2)
        
    def test_restoreMassMatrix(self):
        '''
        test restorMassMatrix

        '''
        # body_4_1 generate maass matrix
        w = np.conj(np.transpose(readData(cwd + '/testData/body_4_test/w.mat'))) 
        waveDir = [0]
        CIkt = 601
        CTTime = readData(cwd + '/testData/body_4_test/CTTime.mat')[0]
        dt = 0.1
        rho = 1000
        g = 9.81
        waveType = 'regular'
        waveAmpTime = np.conj(np.transpose(readData(cwd + '/testData/body_4_test/waveAmpTime.mat'))) 
        iBod = 1 #body number
        numBod = [] 
        ssCalc = 0
        nlHydro = 0
        B2B = 0
        numFreq = []
        self.body_4_1.bodyNumber = 1
        self.body_4_1.bodyTotal = [2]
        self.body_4_1.readH5file()
        self.body_4_1.hydroStiffness = np.zeros((6, 6))
        self.body_4_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
        self.body_4_1.linearDamping = np.zeros((6, 6))
        self.body_4_1.mass = 'equilibrium'
        self.body_4_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
        adjMassWeightFun = 5
        self.body_4_1.momOfInertia = [20907301, 21306090.66, 37085481.11];  
        self.body_4_1.adjustMassMatrix(adjMassWeightFun,B2B) #change mass matrix
        self.body_4_1.restoreMassMatrix() #restore mass matrix
        result1 = readData(cwd + '/testData/body_4_test/body1_fAddedMass.mat')
        np.testing.assert_allclose(self.body_4_1.hydroForce['fAddedMass'], result1)
       
    def test_forceAddedMass(self):
        '''
        test forceAddedMAss

        '''
        # random number test to confirm validity
        # B2B = 0 case
        B2B = 0
        acc =np.array([[80, 33, 90, 62, 34, 89],[72, 14, 56, 93, 22, 19]])
        self.body_4_1.bodyNumber = 1
        self.body_4_1.hydroForce['fAddedMass'] = readData(cwd + '/testData/body_4_test/body1_fAddedMass.mat')
        fam = self.body_4_1.forceAddedMass(acc,B2B)
        result1 = readData(cwd + '/testData/body_4_test/body1_fam.mat')
        np.testing.assert_allclose(fam, result1)
        
        # random number test to confirm validity
        # B2B = 1 case
        B2B = 1
        acc =np.array([[80, 33, 90, 62, 34, 89],[72, 14, 56, 93, 22, 19]])
        self.body_5_2.bodyNumber = 2
        self.body_5_2.hydroForce['fAddedMass'] = readData(cwd + '/testData/body_5_test/body2_fAddedMass.mat')
        fam = self.body_5_2.forceAddedMass(acc,B2B)
        result1 = readData(cwd + '/testData/body_5_test/body2_fam.mat')
        np.testing.assert_allclose(fam, result1)
        
    def test_rotateXYZ(self):
        '''
        test rotateXYZ

        '''
        # random number test to confirm validity
        x = [1, 1, 1]
        ax = [4, 2, -1]
        t = 30
        xn = self.body_4_1.rotateXYZ(x,ax,t)
        result1 = [14.1051275798573, 13.5518950714761, -6.05055454886022]
        np.testing.assert_allclose(xn, result1)
        
        # random number test to confirm validity
        x = [3, -3, -1]
        ax = [2, 7, -4]
        t = 15
        xn = self.body_4_1.rotateXYZ(x,ax,t)
        result2 = [-53.3476667844558, -139.719783953124, 60.6281843944048]
        np.testing.assert_allclose(xn, result2)
        
    def test_setInitDisp(self):
        '''
        test setInitDisp

        '''
        # random number test to confirm validity
        x_rot = np.array([1, 1, 1])
        ax_rot = np.array([4, 2, -1])
        ang_rot = 30
        addLinDisp = 1
        self.body_4_1.cg = np.array([5 , 5, 5])
        self.body_4_1.setInitDisp(x_rot, ax_rot, ang_rot, addLinDisp)
        result1 = [53.4205103194293, 51.2075802859042, -27.2022181954409]
        np.testing.assert_allclose(self.body_4_1.initDisp['initLinDisp'], result1)
        result2 = [4, 2, -1]
        np.testing.assert_allclose(self.body_4_1.initDisp['initAngularDispAxis'], result2)
        result3 = [30]
        np.testing.assert_allclose(self.body_4_1.initDisp['initAngularDispAngle'], result3)
        
        # random number test to confirm validity
        x_rot = np.array([3, -3, -1])
        ax_rot = np.array([2, 7, -4])
        ang_rot = 15
        addLinDisp = 3
        self.body_4_1.cg = np.array([0 , 0, -10])
        self.body_4_1.setInitDisp(x_rot, ax_rot, ang_rot, addLinDisp)
        result1 = [154.602551002163, 645.438156356736, -322.581371323228]
        np.testing.assert_allclose(self.body_4_1.initDisp['initLinDisp'], result1)
        result2 = [2, 7, -4]
        np.testing.assert_allclose(self.body_4_1.initDisp['initAngularDispAxis'], result2)
        result3 = [15]
        np.testing.assert_allclose(self.body_4_1.initDisp['initAngularDispAngle'], result3)
        
if __name__ == '__main__':
    unittest.main()