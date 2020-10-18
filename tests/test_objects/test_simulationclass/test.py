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
        self.body_1 = BodyClass(cwd + '/testData/hydroData/rm3.h5') #regularCIC
        
    def tearDown(self):
        print("tearDown\n")
    
    
if __name__ == "__main__":
    unittest.main()