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

class SSCIandConstantDampingCoeVariantSubsystemClass:
    # This class contains WEC-Sim simulation parameters and settings
    def constantCoefficients(fDamping,inputA):
        """
        Matrix multiplication

        Parameters
        ----------
        fDamping : array
            body_hydroForce_fDamping
        inputA : array
            DESCRIPTION.

        Returns
        -------
        outputA : array
            F_Single Frequency

        """
        outputA = inputA * fDamping
        return(outputA)
        
        
    def convolutionIntegeralCalculation(irkb,dof,CTTime,inputA,dof_gbm):
        """
        

        Parameters
        ----------
        irkb : TYPE
            body_hydroForce_irkb
        dof : TYPE
            body_dof
        CTTime : TYPE
            simu_CTTime
        inputA : TYPE
            v
        dof_gbm : 
            body_dof_gbm

        Returns
        -------
        None.

        """
        v = inputA[:dof_gbm]
        # Function to calculate convolution integral
        
        # define persistent variables
        interp_factor = 1
        
        try: 
            velocity
        except AttributeError:
            CTTime_interp  = arange_MATLAB(0, CTTime, interp_factor)
            velocity       = np.zeros(lenJ,np.size(CTTime_interp))
            IRKB_reordered = np.roll(irkb,1,1)                    # permutes from [601 6 12] to [12 601 6]
            IRKB_interp    = IRKB_reordered[:,arange_MATLAB(0, irkb[0][-1], interp_factor), :]  # what does this do?
        
        # shift velocity
        velocity      = np.roll(velocity,1,1)
        velocity[0,:] = v
        
        # integrate
        time_series = IRKB_interp * velocity
        F_FM = np.squeeze(np.trapz(CTTime_interp,np.sum(time_series, 1)))
        
        def stateSpaceCalculation(A,B,C,D):
        """
        

        Parameters
        ----------
        A : TYPE
            body_hydroForce_ssRadf_A
        B : TYPE
            body_hydroForce_ssRadf_B
        C : TYPE
            body_hydroForce_ssRadf_C
        D : TYPE
            body_hydroForce_ssRadf_D

        Returns
        -------
        OutputA

        """
        outputA = control.ss(A,B,C,D)
        return(outputA)
    