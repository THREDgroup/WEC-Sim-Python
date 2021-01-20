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

class NonlinearFKForceVarientSubsystemClass:
    # This class contains WEC-Sim simulation parameters and settings
    def noNonlinearFroudeKrylovForce(dof):
        """
        non Linear wave excitation.

        Parameters
        ----------
        dof : int
            degree of freedom defined in body class
            
        Returns
        -------
        NE - array
            array of zeros

        """
        NE = np.zeros(1,dof)
        return(NE)
    
    def offsetXYZ(verts,x):
        """
        Function to move the position vertices

        Parameters
        ----------
        verts : TYPE
            DESCRIPTION.
        x : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        verts[1] = verts[1] + x[1]
        verts[2] = verts[2] + x[2]
        verts[3] = verts[3] + x[3]
        return(verts)
        
    def rotateXYZ(x,ax,t):
        """
        % Function to rotate a point about an arbitrary axis
        % x: 3-componenet coordinates
        % ax: axis about which to rotate (must be a vector length 1)
        % t: rotation angle
        % xn: new coordinates after rotation

        Parameters
        ----------
        x : TYPE
            DESCRIPTION.
        ax : TYPE
            DESCRIPTION.
        t : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        rotMat = np.zeros(3)
        rotMat[1,1] = ax[1]*ax[1]*(1-cos(t))    + cos(t)
        rotMat[1,2] = ax[2]*ax[1]*(1-cos(t))    + ax[3]*sin(t)
        rotMat[1,3] = ax[3]*ax[1]*(1-cos(t))    - ax[2]*sin(t)
        rotMat[2,1] = ax[1]*ax[2]*(1-cos(t))    - ax[3]*sin(t)
        rotMat[2,2] = ax[2]*ax[2]*(1-cos(t))    + cos(t)
        rotMat[2,3] = ax[3]*ax[2]*(1-cos(t))    + ax[1]*sin(t)
        rotMat[3,1] = ax[1]*ax[3]*(1-cos(t))    + ax[2]*sin(t)
        rotMat[3,2] = ax[2]*ax[3]*(1-cos(t))    - ax[1]*sin(t)
        rotMat[3,3] = ax[3]*ax[3]*(1-cos(t))    + cos(t)
        xn = x*rotMat
        return(xn)

        def pDis(center,elv,AH,w,dw,wDepth,deepWaterWave,t,k,phaseRand,typeNum,rho,g):
            """
            % Function to calculate pressure distribution            

            Parameters
            ----------
            center : TYPE
                DESCRIPTION.
            elv : TYPE
                DESCRIPTION.
            AH : TYPE
                DESCRIPTION.
            w : TYPE
                DESCRIPTION.
            dw : TYPE
                DESCRIPTION.
            wDepth : TYPE
                DESCRIPTION.
            deepWaterWave : TYPE
                DESCRIPTION.
            t : TYPE
                DESCRIPTION.
            k : TYPE
                DESCRIPTION.
            phaseRand : TYPE
                DESCRIPTION.
            typeNum : TYPE
                DESCRIPTION.
            rho : TYPE
                DESCRIPTION.
            g : TYPE
                DESCRIPTION.

            Returns
            -------
            None.

            """
            f = np.zeros(np.size(center[3]),1)
            z= np.zeros(np.size(center[1]),1)
            if typeNum < 10:
            elif typeNum < 20:
                f = rho*g*AH[1]*np.cos(k[1]*center[1]-w[1]*t)
                if deepWaterWave == 0:
                    z = (center[3]-elv)*wDepth/(wDepth+elv)
                    f = f*(np.cosh(k[1]*(z+wDepth))/np.cosh(k[1]*wDepth))
                else
                    z=(center[3]-elv)
                    f = f*np.exp(k[1]*z)
            elif typeNum < 30:
                for i in range(AH):
                    if deepWaterWave == 0 and wDepth <= 0.5*np.pi/k[i]
                        z=(center[3]-elv)*wDepth/(wDepth+elv)
                        f_tmp = rho*g*np.sqrt(AH[i]*dw[i])*np.cos(k[i]*center[1]-w[i]*t-phaseRand[i])
                        f = f + f_tmp.*(cosh(k(i).*(z+wDepth))./cosh(k(i).*wDepth))
                    else
                        z=(center(:,3)-elv)
                        f_tmp = rho*g*np.sqrt(AH[i]*dw[i]).*np.cos(k[i]*center[1]-w[i]*t-phaseRand[i])
                        f = f + f_tmp*np.exp(k[i].*z)
            if z>0 = 0:
                f = 0
            
        def fk(center,instcg,av,wp):
            """
            

            Parameters
            ----------
            center : TYPE
                DESCRIPTION.
            instcg : TYPE
                DESCRIPTION.
            av : TYPE
                DESCRIPTION.
            wp : TYPE
                DESCRIPTION.

            Returns
            -------
            None.

            """
            # Function to calculate the force and moment about the cog due to Froude-Krylov pressure
            f = np.zeros(6,1)
            
            # Calculate the hydrostatic pressure at each triangle center
            pressureVect = [wp wp wp]*-av
            
            # Compute force about cog
            f(1:3) = np.sum(pressureVect)
            
            # Compute moment about cog
            tmp1 = np.ones(np.size(center(:,1)),1)
            tmp2 = tmp1*np.transpose(instcg)
            center2cgVec = center-tmp2    
            
            f(4:6)= np.sum(np.cross(center2cgVec,pressureVect))

    
    def withNonlinearFroudeKrylovForce(displacement,hydroData_properties_cg,dof,elv,center,tnorm,area,rho,g,cg,AH,w,dw,wDepth,deepWaterWave,k,typeNum,t,phaseRand):
        """
        Function to calculate wave exitation force and moment on a triangulated surface
        NOTE: This function assumes that the STL file is imported with its CG at 0,0,0
        
        
        Parameters
        ----------
        displacement : TYPE
            DESCRIPTION.
        hydroData_properties_cg : TYPE
            DESCRIPTION.
        dof : TYPE
            DESCRIPTION.
        elv : TYPE
            DESCRIPTION.
        center : TYPE
            DESCRIPTION.
        tnorm : TYPE
            DESCRIPTION.
        area : TYPE
            DESCRIPTION.
        rho : TYPE
            DESCRIPTION.
        g : TYPE
            DESCRIPTION.
        cg : TYPE
            DESCRIPTION.
        AH : TYPE
            DESCRIPTION.
        w : TYPE
            DESCRIPTION.
        dw : TYPE
            DESCRIPTION.
        wDepth : TYPE
                DESCRIPTION.
        deepWaterWave : TYPE
            DESCRIPTION.
        k : TYPE
            DESCRIPTION.
        typeNum : TYPE
            DESCRIPTION.
        t : TYPE
            DESCRIPTION.
        phaseRand : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        x = displacement - np.zeros(1,dof-3) * hydroData_properties_cg 
        # Compute new tri coords after cog rotation and translation
        centerMeanFS = offsetXYZ(center,cg)
        avMeanFS     = tnorm * [area, area, area]
        
        # Compute new tri coords after cog rotation and translation
        center = rotateXYZ(center,[1, 0, 0],x[4])
        center = rotateXYZ(center,[0, 1, 0],x[5])
        center = rotateXYZ(center,[0, 0, 1],x[6])
        center = offsetXYZ(center,x)
        center = offsetXYZ(center,cg)
        # Compute new normal vectors coords after cog rotation and translation
        tnorm = rotateXYZ(tnorm,[1, 0, 0],x[4])
        tnorm = rotateXYZ(tnorm,[0, 1, 0],x[5])
        tnorm = rotateXYZ(tnorm,[0, 0, 1],x[6])
        
        # Copute area vectors
        av = tnorm * [area, area, area]
        
        # Calculate the free surface
        wpMeanFS=pDis(centerMeanFS,       0,AH,w,dw,wDepth,deepWaterWave,t,k,phaseRand,typeNum,rho,g)
        wp      =pDis(center      ,     elv,AH,w,dw,wDepth,deepWaterWave,t,k,phaseRand,typeNum,rho,g)
        
        # Calculate forces
        f_linear   =FK(centerMeanFS,       cg,avMeanFS,wpMeanFS)
        f_nonLinear=FK(center      ,x(1:3)+cg,av      ,wp      )
        f = f_nonLinear-f_linear
        return(f, wp, wpMeanFS)
            