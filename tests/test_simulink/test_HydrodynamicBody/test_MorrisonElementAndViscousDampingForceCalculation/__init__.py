#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 08:01:43 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
import numpy as np
import control

def arange_MATLAB(start, end, step):
    """
    Change np.arange to have same sequence as MATLAB when step is float
    """
    return step*np.arange(np.floor(start/step), np.floor(end/step))

class MorrisonElementAndViscousDampingForceCalculationClass:
    # This class contains WEC-Sim simulation parameters and settings
    def constantCoefficients(fDamping,inputA):
        """
       