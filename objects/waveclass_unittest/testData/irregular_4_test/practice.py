#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 23:53:36 2020

@author: logical
"""
import numpy as np
rampTime = 0.2
dt = 0.1
maxIt = 4
dw = np.array([0.0408805215286940,0.0154571054538299,0.00812223235643605,0.0103185502895540])
waveDir = [0]
wA = np.array([0.253488803947639,0.354610720753343,0.413569577232088,0.492917523052480])
waveSpread = [1]
wF = np.array([0.566624163384781,0.582081268838611,0.590203501195047,0.600522051484601])
phase = np.array([2.62022653271779,4.52593227359735,0.000718638171852741,1.89961157824218])
wk = np.array([0.0327281286984203,0.0345380839482943,0.0355086822449431,0.0367611349968679])
wavegauge1loc = [0,0]
wavegauge2loc = [0,0]
wavegauge3loc = [0,0]
waveAmpTime = [np.arange(maxIt+1)*dt,np.zeros((maxIt+1))]
waveAmpTime1 = [np.arange(maxIt+1)*dt,np.zeros((maxIt+1))]
waveAmpTime2 = [np.arange(maxIt+1)*dt,np.zeros((maxIt+1))]
waveAmpTime3 = [np.arange(maxIt+1)*dt,np.zeros((maxIt+1))]
maxRampIT=int(np.round(rampTime/dt))
if rampTime == 0:
    for i in range(maxIt+1):
        for idir in range(np.size(waveDir)):
            t       = (i)*dt
            tmp     = np.sqrt(wA*dw*waveSpread[idir])
            tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(wF*t + phase[idir])))
            tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge1loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[idir])))
            tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge2loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[idir])))
            tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge3loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[idir])))
            waveAmpTime[1][i]    = waveAmpTime[1][i] + sum(tmp1)
            waveAmpTime1[1][i]   = waveAmpTime1[1][i] + sum(tmp11)
            waveAmpTime2[1][i]   = waveAmpTime2[1][i] + sum(tmp12)
            waveAmpTime3[1][i]   = waveAmpTime3[1][i] + sum(tmp13)
else:
    for i in range(maxRampIT):
        for idir in range(np.size(waveDir)):
            t       = (i)*dt
            tmp     = np.sqrt(wA*dw*waveSpread[idir])
            tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(wF*t + phase[idir])))
            tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge1loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[idir])))
            tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge2loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[idir])))
            tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge3loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[idir])))
            waveAmpTime[1][i]    = waveAmpTime[1][i] + sum(tmp1)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
            waveAmpTime1[1][i]   = waveAmpTime1[1][i] + sum(tmp11)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
            waveAmpTime2[1][i]   = waveAmpTime2[1][i] + sum(tmp12)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
            waveAmpTime3[1][i]   = waveAmpTime3[1][i] + sum(tmp13)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
    for i in range(maxRampIT, maxIt+1):
        for idir in range (np.size(waveDir)):
            t       = (i)*dt
            tmp     = np.sqrt(wA*dw*waveSpread[idir])
            tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(wF*t + phase[idir])))
            tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge1loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[idir])))
            tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge2loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[idir])))
            tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge3loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[idir])))
            waveAmpTime[1][i]    = waveAmpTime[1][i] + sum(tmp1)
            waveAmpTime1[1][i]   = waveAmpTime1[1][i] + sum(tmp11)
            waveAmpTime2[1][i]   = waveAmpTime2[1][i] + sum(tmp12)
            waveAmpTime3[1][i]   = waveAmpTime3[1][i] + sum(tmp13)
print(waveAmpTime[1])
