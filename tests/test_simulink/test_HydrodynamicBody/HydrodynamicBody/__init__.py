#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 08:01:43 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
import numpy as np
from WaveDiffractionAndExcitationForceCalculation import NonlinearFKForceVarientSubsystemClass
from WaveDiffractionAndExcitationForceCalculation import LinearWaveExcitationForceVarientSubsystemClass
from WaveRadiationForceCalculation import SSCIandConstantDampingCoeVariantSubsystemClass

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
        if wType == 'noWave' and simu.yawNonLin == 1:
            fWave2 = noWaveExcitationForce(dof)            
        elif wType == 'regularWave' and simu.yawNonLin == 1:
            fWave2 = regularWaveExcitationForce(time,A,w,hydroForce_fExt_md,hydroForce_fExt_re,hydroForce_fExt_im)
        elif wType == 'regularWave' and simu.yawNonLin != 1:
            fWave2 = regularWaveNonLinearYaw(A,w,dofGRD,dirGRD,fEHRE,fEHIM, fEHMD,time,WaveDir,Disp, intThresh, prevYaw, prevCoeffMD, prevCoeffRE, prevCoeffIM)
        elif wType == 'irregularWave' and simu.yawNonLin == 1:     
            fWave2 = irregularWaveExcitationForce(A,w,fExtRE,fExtIM,phaseRand,dw,time,WaveDir,WaveSpread,fExtMD)
        elif wType == 'irregularWave' and simu.yawNonLin != 1:
            fWave2 = irregularWaveNonLinearYaw(A,w,dofGRD,dirGRD,wGRD,fEHRE,fEHIM, fEHMD, phaseRand,dw,time,WaveDir,WaveSpread, Disp, intThresh, prevYaw, prevCoeffMD, prevCoeffRE, prevCoeffIM)
        fWave = fWave1 + fWave2
        fExcitation = fWave * R
    
    def waveRadiationForceCalculation(velocity,time,acceleration):
        """
        

        Parameters
        ----------
        velocity : TYPE
            DESCRIPTION.
        time : TYPE
            DESCRIPTION.
        acceleration : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        F_RadiationDamping = constantCoefficients(velocity,time)
        
        F_AddedMass = fAddedMass * acceleration 
        
        