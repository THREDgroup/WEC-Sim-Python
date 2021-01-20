#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 08:01:43 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
import numpy as np

class SSCIandConstantDampingCoeVariantSubsystemClass:
    # This class contains WEC-Sim simulation parameters and settings
    def noWaveExcitationForce(dof):
        """
        Linear wave excitation for no wave condition.

        Parameters
        ----------
        dof : int
            degree of freedom defined in body class
            
        Returns
        -------
        LE - array
            array of zeros

        """
        LE = np.zeros(1,dof)
        return(LE)
    
    