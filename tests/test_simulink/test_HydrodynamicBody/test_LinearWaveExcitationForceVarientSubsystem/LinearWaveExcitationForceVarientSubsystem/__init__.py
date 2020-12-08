#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 08:01:43 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
import numpy as np

def arange_MATLAB(start, end, step):
    """
    Change np.arange to have same sequence as MATLAB when step is float
    """
    return step*np.arange(np.floor(start/step), np.floor(end/step))

class LinearWaveExcitationForceVarientSubsystemClass:
    # This class contains WEC-Sim simulation parameters and settings
    def properties(self):
        
        self.body_dof = 0 # (`integer`) Bodyclass degree of freedom. Default = ``0``
        self.waves_A = []
        self.body_hydroForce_fExt_md = []

    def __init__(self):
        """
        Initialize LinearWaveExcitationForceVarientSubsystemClass
        Returns
        -------
        None.

        """
        self.properties()

    def noWaveExcitationForce(self):
        """
        Linear wave excitation for no wave condition.

        Returns
        -------
        output - array of zeros

        """
        self.output = np.zeros(1,self.body_dof)
        return(output)
    
    def regularWaveExcitationForce(self,time):
        """
        Linear wave excitation for regular wave condition.

        Parameters
        ----------
        time : TYPE - array
            DESCRIPTION. time step user deifened in wecSimPythonInputFile

        Returns
        -------
        output - 

        """
        waveAmplitude = self.waves_A * self.waves_A
        realExcitation = waveAmplitude * self.body_hydroForce_fExt_md
        
        


   
        
