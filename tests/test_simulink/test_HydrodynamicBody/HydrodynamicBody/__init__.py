#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 08:01:43 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
import numpy as np
import NonlinearFKForceVarientSubsystemClass
import LinearWaveExcitationForceVarientSubsystemClass

def arange_MATLAB(start, end, step):
    """
    Change np.arange to have same sequence as MATLAB when step is float
    """
    return step*np.arange(np.floor(start/step), np.floor(end/step))

class HydrodynamicBodyClass:
    def rampFunction(time,rampTime):
        """
        

        Parameters
        ----------
        time : TYPE
            DESCRIPTION.
        rampTime : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        frequency = np.pi / rampTime
        t = time * frequency + 3*np.pi/2
        r = 0.5 * (1+np.sin(t))
        return(r)
        
    def waveDiffractionAndExcitationForceCalculation(rampTime,time,displacement,waveElv):
        """
        

        Parameters
        ----------
        rampTime : TYPE
            DESCRIPTION.
        time : TYPE
            DESCRIPTION.
        displacement : TYPE
            DESCRIPTION.
        waveElv : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        for i in range(np.size(time)):
            if(time[i] < rampTime[i]):
                R = rampFunction(time[i],rampTime[i])
            else:
                R = 1 
        fWave1 = withNonlinearFroudeKrylovForce(displacement,hydroData_properties_cg,dof,elv,center,tnorm,area,rho,g,cg,AH,w,dw,wDepth,deepWaterWave,k,typeNum,t,phaseRand)
        if wType = 'regularWave':
            fWave2 = regularWaveExcitationForce(time,A,w,hydroForce_fExt_md,hydroForce_fExt_re,hydroForce_fExt_im)
        elif wType = 'irregularWave':     
            fWave2 = irregularWaveExcitationForce(A,w,fExtRE,fExtIM,phaseRand,dw,time,WaveDir,WaveSpread,fExtMD)
        fWave = fWave1 + fWave2
        fExcitation = fWave * R