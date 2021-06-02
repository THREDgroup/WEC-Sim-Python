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

class LinearAndNonlinearRestoringForceVarientSubsystemClass:
    # This class contains WEC-Sim simulation parameters and settings
    def linearHydrostaticRestoringForce():
        f_gravity = g * mass
        f_buoyancy = rho * g * dispVol
        U = f_gravity - f_buoyancy
        netBuoyancyForce = U #assignment block is used need test
        
        
    def nonLinearHydrostaticRestoringForce(x,elv,center,tnorm,area,rho,g,cg,mass):
        """
        Function to calculate buoyancy force and moment on a triangulated surface
        NOTE: This function assumes that the STL file is imported with its CG at 0,0,0

        Parameters
        ----------
        x : TYPE
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
        mass : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        [f,p] = calc_buoyancy(x,elv,center,tnorm,area,rho,g,cg);
        f = -f + [0 0 g*mass 0 0 0]';
        end
        
        
        function [f,p] = calc_buoyancy(x,elv,center,tnorm,area,rho,g,cg)
        % Function to apply translation and rotation and calculate forces
        
        % Compute new tri coords after cog rotation and translation
        center = rotateXYZ(center,[1 0 0],x(4));
        center = rotateXYZ(center,[0 1 0],x(5));
        center = rotateXYZ(center,[0 0 1],x(6));
        center = offsetXYZ(center,x);
        center = offsetXYZ(center,cg);
        % Compute new normal vectors coords after cog rotation
        tnorm = rotateXYZ(tnorm,[1 0 0],x(4));
        tnorm = rotateXYZ(tnorm,[0 1 0],x(5));
        tnorm = rotateXYZ(tnorm,[0 0 1],x(6));
        
        % Calculate the hydrostatic forces
        av = tnorm .* [area area area];
        [f,p]=fHydrostatic(center,elv,x(1:3)+cg,av,rho,g);
        end
        
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
        
        function [f,p] = fHydrostatic(center,elv,instcg,av,rho,g)
        % Function to calculate the force and moment about the cog due to hydrostatic pressure
        f = zeros(6,1);
        
        % Zeor out regions above the mean free surface
        z=center(:,3); z((z-elv)>0)=0;
        
        % Calculate the hydrostatic pressure at each triangle center
        pressureVect = rho*g.*[-z -z -z].*-av;
        p = rho*g.*-z;
        % Compute force about cog
        f(1:3) = sum(pressureVect);
        
        tmp1 = ones(length(center(:,1)),1);
        tmp2 = tmp1*instcg';
        center2cgVec = center-tmp2;
        
        f(4:6)= sum(cross(center2cgVec,pressureVect));
        end
        