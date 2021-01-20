#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 08:01:43 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
import numpy as np

class LinearWaveExcitationForceVarientSubsystemClass:
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
    
    def regularWaveExcitationForce(time,A,w,hydroForce_fExt_md,hydroForce_fExt_re,hydroForce_fExt_im):
        """
        Linear wave excitation for regular wave condition.
        Calculates the wave force, F_wave for the case of regular wave.
        F_wave = A*cos(w*t)*Re-A*sin(w*t)*Im

        Parameters
        ----------
        time : array
            time step user deifened in wecSimPythonInputFile
        A : array
            wave amplitude
        w : array
            wave frequency
        hydroForce_fExt_md : array
            body real excitataion
        hydroForce_fExt_re : array
            body real excitataion
        hydroForce_fExt_im : array
            body real excitataion

        Returns
        -------
        LE - array
            linear excitation

        """
        A_A_Re = A * A * hydroForce_fExt_md         
        A_Re_coswt = A * hydroForce_fExt_re * np.cos(w*time)        
        A_Im_sinwt = A * hydroForce_fExt_im * np.sin(w*time)        
        LE = A_A_Re + A_Re_coswt - A_Im_sinwt
        return(LE)
    
    def regularWaveNonLinearYaw(A,w,dofGRD,dirGRD,fEHRE,fEHIM, fEHMD,time,WaveDir,Disp, intThresh, prevYaw, prevCoeffMD, prevCoeffRE, prevCoeffIM):
        """
        

        Parameters
        ----------
        A : array
            wave amplitude
        w : array
            wave frequency
        dofGRD : TYPE
            DESCRIPTION.
        dirGRD : TYPE
            DESCRIPTION.
        fEHRE : TYPE
            DESCRIPTION.
        fEHIM : TYPE
            DESCRIPTION.
        fEHMD : TYPE
            DESCRIPTION.
        time : TYPE
            DESCRIPTION.
        WaveDir : TYPE
            DESCRIPTION.
        Disp : TYPE
            DESCRIPTION.
        intThresh : TYPE
            DESCRIPTION.
        prevYaw : TYPE
            DESCRIPTION.
        prevCoeffMD : TYPE
            DESCRIPTION.
        prevCoeffRE : TYPE
            DESCRIPTION.
        prevCoeffIM : TYPE
            DESCRIPTION.

        Returns
        -------
        Fext array
            excitation force output
        relYawLast
            the angle of relative yaw at which that last interpolation was performed.
            If the current relative yaw angle - last interpolated > threshold, the 
            interpolation is performed again. To interpolate every time step, let threshold=0.
        coeffsLastMD,
        coeffsLastRE,
        coeffsLastIM

        """
        Fext = np.zeros(1,6)
        relYawLast = 0
        coeffsLastMD = np.zeros(np.size(w),6) # should be (1,6) for regular waves
        coeffsLastRE = np.zeros(np.size(w),6)
        coeffsLastIM = np.zeros(np.size(w),6)
        for ii in range(WaveDir): #should be length=1 for regular waves
            relYaw = WaveDir[ii]-(Disp[5]*180/np.pi) # relative yaw angle, convert to deg
            
            # compare relYaw to available WaveDir data
            if relYaw > max(dirGRD[1,:]): # handle interpolation out of BEM range by wrapping other bound
                [~,idx] = min(dirGRD,[],2)
                dirGRD[:,idx] = dirGRD[:,idx]+360
                [dirGRD,I] = np.sort(dirGRD,2)# grid must be in ascending order
            elif relYaw < min(dirGRD[1,:]):
                [~,idx] = max(dirGRD,[],2) 
                dirGRD[:,idx] = dirGRD[:,idx]-360
                [dirGRD,I] = np.sort(dirGRD,2)# grid must be in ascending order
            else:
                I = range(np.size(dirGRD(1,:)))
                    
            if abs(relYaw-prevYaw)>intThresh and min(abs(relYaw-dirGRD(1,:)))>intThresh:# interpolate for nonlinear yaw 
               
                # performs 1D interpolation in wave direction
                fExtMDint = interpn(dofGRD,dirGRD,fEHMD(:,I(1,:)),dofGRD,relYaw*ones(np.size(dirGRD)))
                fExtREint = interpn(dofGRD,dirGRD,fEHRE(:,I(1,:)),dofGRD,relYaw*ones(np.size(dirGRD)))
                fExtIMint = interpn(dofGRD,dirGRD,fEHIM(:,I(1,:)),dofGRD,relYaw*ones(np.size(dirGRD)))
                
                fExtMDint = fExtMDint(:,1).'
                fExtREint = fExtREint(:,1).'
                fExtIMint = fExtIMint(:,1).'
                        
                relYawLast = relYaw
                coeffsLastMD = fExtMDint
                coeffsLastRE = fExtREint
                coeffsLastIM = fExtIMint
          
           elif min(abs(relYaw-dirGRD(1,:)))>intThresh: # significant yaw is present, but close to previous value
                fExtMDint = prevCoeffMD
                fExtREint = prevCoeffRE
                fExtIMint = prevCoeffIM
                
                relYawLast = prevYaw
                coeffsLastMD = prevCoeffMD
                coeffsLastRE = prevCoeffRE
                coeffsLastIM = prevCoeffIM
                      
            else:    # significant yaw may be is present, but nearby BEM calculated data
                [~,idx] = min(abs(relYaw-dirGRD(1,:)))
                fExtMDint = fEHMD(:,I(1,idx)).' # maintains dimensional consistency
                fExtREint = fEHRE(:,I(1,idx)).'
                fExtIMint = fEHIM(:,I(1,idx)).'
             
                relYawLast = dirGRD(1,idx) # indices of dirGRD have already been resorted
                coeffsLastMD = fExtMDint
                coeffsLastRE = fExtREint
                coeffsLastIM = fExtIMint
        
            # regular wave excitation equation (with mean drift)
            Fext = fExtMDint*A**2 + A*fExtREint*cos(w*time)- A*fExtIMint*sin(w*time)
        return(Fext,relYawLast,coeffsLastMD,coeffsLastRE,coeffsLastIM)

    def irregularWaveExcitationForce(A,w,fExtRE,fExtIM,phaseRand,dw,time,WaveDir,WaveSpread,fExtMD):
        """
        Calculates the wave force for irregular wave
        F_wave = sum(F_wave[i])
        where i = each frequency bins

        Parameters
        ----------
        A : array
            wave amplitude
        w : array
            wave frequency
        fExtRE : TYPE
            DESCRIPTION.
        fExtIM : TYPE
            DESCRIPTION.
        phaseRand : TYPE
            DESCRIPTION.
        dw : TYPE
            DESCRIPTION.
        time : TYPE
            DESCRIPTION.
        WaveDir : TYPE
            DESCRIPTION.
        WaveSpread : TYPE
            DESCRIPTION.
        fExtMD : TYPE
            DESCRIPTION.

        Returns
        -------
        Fext array
            excitation force output

        """
        # pversistent A1 B1 B11 C1 D1 D11 E1 E11

        A1=(w*time)+(np.pi/2)
        Fext = np.zeros(1,6)
        for ii in range(WaveDir):
            B1 = np.sin(A1 + phaseRand[ii])
            B11 = np.sin(w*time+phaseRand[ii])
            C0 = A*WaveSpread[ii]*dw
            C1 = np.sqrt(A*WaveSpread[ii]*dw)
            D0 = np.squeeze(fExtMD[ii]*C0)
            D1 = np.squeeze(fExtRE[ii]*C1)
            D11 = np.squeeze(fExtIM[ii]*C1)
            E1 = D0+ B1*D1
            E11 = B11*D11
            Fext = Fext + sum(E1-E11)
        return(Fext)
    
    def irregularWaveNonLinearYaw(A,w,dofGRD,dirGRD,wGRD,fEHRE,fEHIM, fEHMD, phaseRand,dw,time,WaveDir,WaveSpread, Disp, intThresh, prevYaw, prevCoeffMD, prevCoeffRE, prevCoeffIM):
        """
        

        Parameters
        ----------
        A : TYPE
            DESCRIPTION.
        w : TYPE
            DESCRIPTION.
        dofGRD : TYPE
            DESCRIPTION.
        dirGRD : TYPE
            DESCRIPTION.
        wGRD : TYPE
            DESCRIPTION.
        fEHRE : TYPE
            DESCRIPTION.
        fEHIM : TYPE
            DESCRIPTION.
        fEHMD : TYPE
            DESCRIPTION.
        phaseRand : TYPE
            DESCRIPTION.
        dw : TYPE
            DESCRIPTION.
        time : TYPE
            DESCRIPTION.
        WaveDir : TYPE
            DESCRIPTION.
        WaveSpread : TYPE
            DESCRIPTION.
        Disp : TYPE
            DESCRIPTION.
        intThresh : TYPE
            DESCRIPTION.
        prevYaw : TYPE
            DESCRIPTION.
        prevCoeffMD : TYPE
            DESCRIPTION.
        prevCoeffRE : TYPE
            DESCRIPTION.
        prevCoeffIM : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # Fext is the excitation force output, relYawLast is the angle of relative
        # yaw at which that last interpolation was performed. If the current
        # relative yaw angle - last interpolated > threshold, the interpolation is
        # performed again. To interpolate every time step, let threshold=0.
        """
        # caseFlag is for debug, delete when satisfied
        A1=bsxfun(@plus,w*time,pi/2)
        #initialize outputs
        Fext = zeros(1,6)
        relYawLast=0
        coeffsLastMD=zeros(length(w),6)
        coeffsLastRE=zeros(length(w),6)
        coeffsLastIM=zeros(length(w),6)
        
        for ii=1:length(WaveDir)
            relYaw=WaveDir(ii)-(Disp(6)*180/pi) # relative yaw angle
            phaseRandint=phaseRand(:,ii)
            WaveSpreadint=WaveSpread(ii)
            
               # compare relYaw to available WaveDir data
                if relYaw>max(dirGRD(1,:,1)) # handle interpolation out of BEM range by wrapping other bound
                    [~,idx]=min(dirGRD(1,:,1),[],2) # 
                    dirGRD(:,idx,:)=dirGRD(:,idx,:)+360
                    [dirGRD,I]=sort(dirGRD,2) # grid must be in ascending order
                              
                elseif relYaw<min(dirGRD(1,:,1))
                    [~,idx]=max(dirGRD(1,:,1),[],2) # 
                    dirGRD(:,idx,:)=dirGRD(:,idx,:)-360  
                    [dirGRD,I]=sort(dirGRD,2) # grid must be in ascending order
        
                else 
                    I=1:length(dirGRD(1,:,1))
                end
            
            if abs(relYaw-prevYaw)>intThresh && min(abs(relYaw-dirGRD(1,:,1)))>intThresh# interpolate for nonlinear yaw 
                    
                # performs 1D interpolation in wave direction
                fExtMDint=interpn(dofGRD,dirGRD,wGRD,fEHMD(:,I(1,:,1),:),dofGRD,relYaw*ones(size(dirGRD)),wGRD)
                fExtREint=interpn(dofGRD,dirGRD,wGRD,fEHRE(:,I(1,:,1),:),dofGRD,relYaw*ones(size(dirGRD)),wGRD)
                fExtIMint=interpn(dofGRD,dirGRD,wGRD,fEHIM(:,I(1,:,1),:),dofGRD,relYaw*ones(size(dirGRD)),wGRD)
                
                # permute for dimensional consistency with linear yaw
                fExtMDint=squeeze(permute(fExtMDint(:,1,:),[2 3 1]))
                fExtREint=squeeze(permute(fExtREint(:,1,:),[2 3 1]))
                fExtIMint=squeeze(permute(fExtIMint(:,1,:),[2 3 1]))
                
                relYawLast=relYaw
                coeffsLastMD=fExtMDint
                coeffsLastRE=fExtREint
                coeffsLastIM=fExtIMint
                
            elseif min(abs(relYaw-dirGRD(1,:,1)))>intThresh # significant yaw is present, but close to previous value
                fExtMDint=prevCoeffMD
                fExtREint=prevCoeffRE
                fExtIMint=prevCoeffIM
                
                relYawLast=prevYaw
                coeffsLastMD=prevCoeffMD
                coeffsLastRE=prevCoeffRE
                coeffsLastIM=prevCoeffIM
                
            else    # significant yaw may be is present, but nearby BEM calculated data
                [~,idx]=min(abs(relYaw-dirGRD(1,:,1)))
                fExtMDint=squeeze(fEHMD(:,I(1,idx,1),:)).' # maintains dimensional consistency
                fExtREint=squeeze(fEHRE(:,I(1,idx,1),:)).'
                fExtIMint=squeeze(fEHIM(:,I(1,idx,1),:)).'
                
                relYawLast=dirGRD(1,idx,1) # indices of dirGRD have already been resorted
                coeffsLastMD=fExtMDint
                coeffsLastRE=fExtREint
                coeffsLastIM=fExtIMint
            end
            B1= sin(bsxfun(@plus,A1,phaseRandint))
            B11 = sin(bsxfun(@plus,w*time,phaseRandint))
            C0 = bsxfun(@times,A*WaveSpreadint,dw)
            C1 = sqrt(bsxfun(@times,A*WaveSpreadint,dw))
            D0 =bsxfun(@times,fExtMDint,C0)
            D1 =bsxfun(@times,fExtREint,C1)
            D11 = bsxfun(@times,fExtIMint,C1)
            E1 = D0+ bsxfun(@times,B1,D1)
            E11 = bsxfun(@times,B11,D11)
            Fext = Fext + sum(bsxfun(@minus,E1,E11)) 
        """
            
    def userDefinedWaveExcitationForce(waveAmpTime,userDefinedFe):
        """
        Calculates the wave force, F_wave, for the case of user defined Waves
        F_wave = convolution calculation.

        Parameters
        ----------
        waveAmpTime : TYPE
            DESCRIPTION.
        userDefinedFe : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        LE = waveAmpTime + userDefinedFe
