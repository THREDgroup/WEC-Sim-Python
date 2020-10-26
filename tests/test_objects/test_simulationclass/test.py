#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 08:43:15 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""

import unittest
import numpy as np
import scipy.io as sio
import os

from simulationclass import SimulationClass

global cwd # set current directory as global variabl
cwd = os.getcwd()
lastDir = os.path.split(cwd)[1]
if lastDir == 'WEC-Sim-Python': # check if we are in desired directory
    cwd = cwd + '/tests/test_objects/test_simulationclass'

def readData(file): #read MATLAB file
    matFile = sio.loadmat(file) 
    keys = list(matFile.keys())[-1]
    datas = matFile[keys]
    return datas

class TestBody(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        print("setupClass")
    @classmethod    
    def tearDownClass(cls):
        print("teardownClass")
    
    def setUp(self):
        print("setUp")
        self.simulation_1 = SimulationClass() # RM3 
        self.simulation_2 = SimulationClass() # Randome number test
        
    def tearDown(self):
        print("tearDown\n")
    
    def test_setupSim(self):
        self.simulation_1.startTime = 0                     # Simulation Start Time [s]
        self.simulation_1.rampTime = 100                    # Wave Ramp Time [s]
        self.simulation_1.endTime = 400                       # Simulation End Time [s]
        self.simulation_1.dt = 0.1 	
        self.simulation_1.setupSim()
        result1 = np.conj(np.transpose(readData(cwd + '/testData/simulation_1_test/simulation_time.mat')))[:,0]
        self.assertIsNone(np.testing.assert_allclose(self.simulation_1.time, result1))
        result2 = np.conj(np.transpose(readData(cwd + '/testData/simulation_1_test/simulation_CTTime.mat')))[:,0]
        self.assertIsNone(np.testing.assert_allclose(self.simulation_1.CTTime, result2))
        self.assertEqual(self.simulation_1.maxIt, 4000)
        self.assertEqual(self.simulation_1.CIkt, 601)
        
        self.simulation_2.startTime = 477                     # Simulation Start Time [s]
        self.simulation_2.rampTime = 1750                    # Wave Ramp Time [s]
        self.simulation_2.endTime = 6920                       # Simulation End Time [s]
        self.simulation_2.dt = 0.5 	
        self.simulation_2.setupSim()
        result1 = np.conj(np.transpose(readData(cwd + '/testData/simulation_2_test/simulation_time.mat')))[:,0]
        self.assertIsNone(np.testing.assert_allclose(self.simulation_2.time, result1))
        result2 = np.conj(np.transpose(readData(cwd + '/testData/simulation_2_test/simulation_CTTime.mat')))[:,0]
        self.assertIsNone(np.testing.assert_allclose(self.simulation_2.CTTime, result2))
        self.assertEqual(self.simulation_2.maxIt, 12886)
        self.assertEqual(self.simulation_2.CIkt, 121)
    
    
if __name__ == "__main__":
    unittest.main()