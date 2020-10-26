#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 08:01:43 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

"""
import numpy as np
import warnings

class PtoClass:
    # This class contains WEC-Sim contraint parameters and settings
    def inputProperties(self): #input file 
        self.name                    = 'NOT DEFINED'                                 # Name of the constraint used 
        self.loc                     = [999, 999, 999]                                 # Constraint location. Default = [0 0 0]        
        self.orientation             = {                                    # Structure difining axis orientation parameters
                                         'z': [0, 0, 1],                         # Vector defining the direction of the Z-coordinate for the constraint.
                                         'y': [0, 1, 0],                         # Vector defining the direction of the Y-coordinate for the constraint.
                                         'x': [],                                # Internally calculated vector defining the direction of the X-coordinate for the constraint.
                                         'rotationMatrix':[]}                       # Internally calculated rotation matrix to go form standard coordinate orientation to the constraint's coordinate orientation.
        self.initDisp                = {                                    # Structure defining the initial displacement
                                         'initLinDisp':[0, 0, 0]}           # Initial displacement - used for decay tests (format: [displacment in m], default = [0 0 0])
                                
    
    def internalProperties(self):#internal
        self. constraintNum           = []                                            # PTO number
                                            
    def __init__(self,name):
        """
        Initialize constraint class

        """
        self.inputProperties()
        self.internalProperties()
        self.name = name
        
    def checkLoc(self,action):
        # Checks if location is set and outputs a warning or error.
        # Used in mask Initialization.
        if action == 'W':
            if self.loc == 999: # Because "Allow library block to modify its content" is selected in block's mask initialization, this command runs twice, but warnings cannot be displayed during the first initialization. 
                self.loc = [888, 888, 888]
            elif self.loc == 888:
                self.loc = [0, 0, 0]
                warnings.warn('For ', self.name, ': pto.loc was changed from [9999, 9999, 9999] to [0, 0, 0]')
        elif action == 'E':
            if self.loc == 999:
                warnings.warn('For ', self.name, ': pto.loc needs to be specified in the WEC-Sim input file. pto.loc is the [x y z] location, in meters, for the rotational PTO.')

    def setOrientation(self):
        # Sets orientation based on user input
        self.orientation['z'] = self.orientation['z'] / np.norm(self.orientation['z'])
        self.orientation['y'] = self.orientation['y'] / np.norm(self.orientation['y'])
        z = self.orientation['z']
        y = self.orientation['y']
        if np.abs(np.dot(y,z))>0.001:
            warnings.warn('The Y and Z vectors defining the constraint''s orientation must be orthogonal.')
        x = np.cross(y,z)/np.norm(np.cross(y,z))
        x = np.transpose(x)
        self.orientation['x'] = x;
        self.orientation['rotationMatrix']  = [np.transpose(x),np.transpose(y),np.transpose(z)]


    def setInitDisp(self, x_rot, ax_rot, ang_rot, addLinDisp):
        # Function to set the initial displacement when having initial rotation
        # x_rot: rotation point
        # ax_rot: axis about which to rotate (must be a normal vector)
        # ang_rot: rotation angle in radians
        # addLinDisp: initial linear displacement (in addition to the displacement caused by rotation)
        loc = self.loc
        relCoord = loc - x_rot
        rotatedRelCoord = self.rotateXYZ(relCoord,ax_rot,ang_rot)
        newCoord = rotatedRelCoord + x_rot
        linDisp = newCoord - loc
        self.initDisp['initLinDisp'] = linDisp + addLinDisp
    
    def rotateXYZ(self,x,ax,t):
        """
        Function to rotate a point about an arbitrary axis
        x: 3-componenet coordiantes
        ax: axis about which to rotate (must be a normal vector)
        t: rotation angle
        xn: new coordinates after rotation
        
        """
        rotMat = np.zeros((3,3))
        rotMat[0,0] = ax[0]*ax[0]*(1-np.cos(t))    + np.cos(t)
        rotMat[0,1] = ax[1]*ax[0]*(1-np.cos(t))    + ax[2]*np.sin(t)
        rotMat[0,2] = ax[2]*ax[0]*(1-np.cos(t))    - ax[1]*np.sin(t)
        rotMat[1,0] = ax[0]*ax[1]*(1-np.cos(t))    - ax[2]*np.sin(t)
        rotMat[1,1] = ax[1]*ax[1]*(1-np.cos(t))    + np.cos(t)
        rotMat[1,2] = ax[2]*ax[1]*(1-np.cos(t))    + ax[0]*np.sin(t)
        rotMat[2,0] = ax[0]*ax[2]*(1-np.cos(t))    + ax[1]*np.sin(t)
        rotMat[2,1] = ax[1]*ax[2]*(1-np.cos(t))    - ax[0]*np.sin(t)
        rotMat[2,2] = ax[2]*ax[2]*(1-np.cos(t))    + np.cos(t)
        xn = np.dot(x,rotMat)
        return(xn)
    
    def listInfo(self):
        # This method prints pto information to the MATLAB Command Window.
        print('\n\t***** PTO Name: ', self.name,'*****\n')
        print('\tPTO Stiffness           (N/m;Nm/rad) = ', self.k,'\n')
        print('\tPTO Damping           (Ns/m;Nsm/rad) = ', self.c,'\n')
        