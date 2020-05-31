#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 23:21:05 2020

@author: logical
"""


import numpy as np


wType = 'regular'
bemFreq = [5.19999512307279,0.0199999977946844]
wDepth = 'infinite'
rampTime = 100
dt = 0.1
maxIt = 2000
g = 9.81
rho = 1000
endTime = 200
freqRange = []
freqDisc = 'EqualEnergy'
numFreq = 1001
phaseSeed = 0
wT = 8
wH = 2.5
spectrumType = 'NOT DEFINED'
waveDir = [0]
waveSpread = [1]
wavegauge1loc = [0,0]
wavegauge2loc = [0,0]
wavegauge3loc = [0,0]
wF = []
waveAmpTime1 = []
waveAmpTime2 = []
waveAmpTime3 = []

if wDepth == 'infinite':
    deepWaterWave = 1
    waterDepth = 200
    print('\tInfinite water depth specified in BEM, "waves.waterDepth" set to 200m for vizualisation.\n')
else:
    deepWaterWave = 0
    waterDepth = wDepth
if wType == 'noWave':
    wH = 0
elif wType == 'noWaveCIC':
    wH = 0
    if np.size(wF) == 0 and wT == 'NOT DEFINED':
        wF = np.min(bemFreq)
        wT = 2*pi/wF
    elif np.size(wF) == 0:
        wF = 2*pi/wT
    else:
        wT = 2*pi/wF
elif wType == 'spectrumImport':
    wH = 0
    wT = 0
    freqDisc = 'Imported'
    spectrumType = 'spectrumImport'
    
if np.size(wF) == 0:
    wF = 2*np.pi/wT
wA = wH/2
#WaveNumber
wk= wF**2/g
if deepWaterWave == 0:
    for i in range(1,101):
        wk = wF**2/g/np.tanh(wk*waterDepth)

waveAmpTime = [[],[]]
waveAmpTime1 = [[],[]]
waveAmpTime2 = [[],[]]
waveAmpTime3 = [[],[]]
maxRampIT = int(np.round(rampTime/dt))
if rampTime == 0:
    t = np.arange(maxIt+1)*dt
    waveAmpTime[0]    = t
    waveAmpTime[1]    = wA*np.cos(wF*t)
    waveAmpTime1[0]   = t
    waveAmpTime1[1]   = wA*np.cos(wF*t-wk*(wavegauge1loc[0]*np.cos(waveDir[0]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[0]*np.pi/180)))
    waveAmpTime2[0]   = t
    waveAmpTime2[1]   = wA*np.cos(wF*t-wk*(wavegauge2loc[0]*np.cos(waveDir[0]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[0]*np.pi/180)))
    waveAmpTime3[0]   = t;
    waveAmpTime3[1]   = wA*np.cos(wF*t-wk*(wavegauge3loc[0]*np.cos(waveDir[0]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[0]*np.pi/180)))
else:
    t = np.arange(maxRampIT)*dt
    waveAmpTime[0]    = t
    waveAmpTime[1]    = wA*np.cos(wF*t)*(1+np.cos(np.pi+np.pi*np.arange(maxRampIT)/maxRampIT))/2
    waveAmpTime1[0]   = t
    waveAmpTime1[1]   = wA*np.cos(wF*t-wk*(wavegauge1loc[0]*np.cos(waveDir[0]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[0]*np.pi/180)))*(1+np.cos(np.pi+np.pi*np.arange(maxRampIT)/maxRampIT))/2
    waveAmpTime2[0]   = t
    waveAmpTime2[1]   = wA*np.cos(wF*t-wk*(wavegauge2loc[0]*np.cos(waveDir[0]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[0]*np.pi/180)))*(1+np.cos(np.pi+np.pi*np.arange(maxRampIT)/maxRampIT))/2
    waveAmpTime3[0]   = t
    waveAmpTime3[1]   = wA*np.cos(wF*t-wk*(wavegauge3loc[0]*np.cos(waveDir[0]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[0]*np.pi/180)))*(1+np.cos(np.pi+np.pi*np.arange(maxRampIT)/maxRampIT))/2
    t = np.arange(maxRampIT,maxIt+1)*dt
    waveAmpTime[0] = np.append(waveAmpTime[0], t)
    waveAmpTime[1] = np.append(waveAmpTime[1], wA*np.cos(wF*t))
    waveAmpTime1[0]   = t
    waveAmpTime1[1]   = wA*np.cos(wF*t-wk*(wavegauge1loc[0]*np.cos(waveDir[0]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[0]*np.pi/180)))
    waveAmpTime2[0]   = t
    waveAmpTime2[1]   = wA*np.cos(wF*t-wk*(wavegauge2loc[0]*np.cos(waveDir[0]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[0]*np.pi/180)))
    waveAmpTime3[0]   = t
    waveAmpTime3[1]   = wA*np.cos(wF*t-wk*(wavegauge3loc[0]*np.cos(waveDir[0]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[0]*np.pi/180)))
#wavePowerReg
if deepWaterWave == 1:
    # Deepwater Approximation
    wP = 1/(8*np.pi)*rho*g**(2)*(wA)**(2)*wT
else:
    # Full Wave Power Equation
    wP = rho*g*(wA)**(2)/4*np.sqrt(g/wk*np.tanh(k*waterDepth))*(1+2*wk*waterDepth/np.sinh(wk*waterDepth))
        
            