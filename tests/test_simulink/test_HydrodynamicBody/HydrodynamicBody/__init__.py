#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 08:01:43 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
import numpy as np
import os
import warnings

def arange_MATLAB(start, end, step):
    """
    Change np.arange to have same sequence as MATLAB when step is float
    """
    return step*np.arange(np.floor(start/step), np.floor(end/step))

class SimulationClass:
    # This class contains WEC-Sim simulation parameters and settings
    def inputProperties(self):
        
        self.pressureDis         = 0                                            # (`integer`) Option to save pressure distribution: off->0, on->1. Default = ``0``

    def internalProperties(self):
        self.version             = '1.0'                                        # (`string`) WEC-Sim-Python version
        
    def __init__(self):
        """
        Initialize Simulink Class
        Returns
        -------
        None.

        """
        self.inputProperties()
        self.internalProperties()
        print('WEC-Sim-Python: An open-source code for simulating wave energy converters\n')
        self.caseDir = os.getcwd()
        print('\tCase Dir: ',self.caseDir,' \n')

    def setupSim(self):
        """
        Set up simulation property for WEC-Sim-Python
        """
        # Sets simulation properties based on values specified in input file
        self.time = arange_MATLAB(self.startTime,self.endTime+self.dt,self.dt)
        
class SimulationClass:
    # This class contains WEC-Sim simulation parameters and settings
    def inputProperties(self):
        
        self.pressureDis         = 0                                            # (`integer`) Option to save pressure distribution: off->0, on->1. Default = ``0``

    def internalProperties(self):
        self.version             = '1.0'                                        # (`string`) WEC-Sim-Python version
        
    def __init__(self):
        """
        Initialize Simulink Class
        Returns
        -------
        None.

        """
        self.inputProperties()
        self.internalProperties()
        print('WEC-Sim-Python: An open-source code for simulating wave energy converters\n')
        self.caseDir = os.getcwd()
        print('\tCase Dir: ',self.caseDir,' \n')

    def setupSim(self):
        """
        Set up simulation property for WEC-Sim-Python
        """
        # Sets simulation properties based on values specified in input file
        self.time = arange_MATLAB(self.startTime,self.endTime+self.dt,self.dt)

   
        
