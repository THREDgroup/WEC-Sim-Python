#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 23:21:05 2020

@author: logical
"""

from scipy import integrate
import numpy as np


from sys import maxsize
from numpy import set_printoptions

wType = 'spectrumImport'
bemFreq = [5.19999512307279,0.0199999977946844]
wDepth = 'infinite'
rampTime = 100
dt = 0.1
maxIt = 2000
g = 9.81
rho = 1000
endTime = 400
freqRange = []
numFreq = 0
phaseSeed = 1
waveDir = [0]
waveSpread = [1]
wavegauge1loc = [0,0]
wavegauge2loc = [0,0]
wavegauge3loc = [0,0]
wF = []
dw = []
waveAmpTime = []
waveAmpTime1 = []
waveAmpTime2 = []
waveAmpTime3 = []
spectrumDataFile = "spectrumData.txt"

if wDepth == 'infinite':
    deepWaterWave = 1
    waterDepth = 200
    print('\tInfinite water depth specified in BEM, "waves.waterDepth" set to 200m for vizualisation.\n')
else:
    deepWaterWave = 0
    waterDepth = wDepth
wH = 0
wT = 0
freqDisc = 'Imported'
spectrumType = 'spectrumImport'
#waveSetup
WFQSt=np.min(bemFreq)
WFQEd=np.max(bemFreq)

data = np.conj(np.transpose(np.loadtxt(spectrumDataFile)))
freq_data = data[0]
wF = np.array([i*2*np.pi for i in freq_data if i>=min(bemFreq)/2/np.pi and i<=max(bemFreq)/2/np.pi])
numFreq = len(wF)
dw = np.zeros(numFreq)
dw[0]= np.array(wF[1]-wF[0])
dw[1:numFreq-1]= np.array((wF[2:]-wF[:-2])/2)
dw[numFreq-1]= np.array(wF[-1]-wF[-2])


#setWavePhase
if phaseSeed != 0:
    np.random.seed(phaseSeed) #Phase seed = 1,2,3,...,etc
else:
    np.random.seed(np.random.shuffle(phaseSeed))
if (freqDisc == 'EqualEnergy') or (freqDisc == 'Traditional'):
    phase = 2*np.pi*np.random.rand(1,numFreq)[0]
elif (freqDisc == 'Imported'):
    data = np.conj(np.transpose(np.loadtxt(spectrumDataFile)))
    if len(data) == 3:
        freq_data = data[0]
        phase = np.array([x for x,i in zip(data[2],freq_data) if i>=min(bemFreq)/2/np.pi and i<=max(bemFreq)/2/np.pi])
    else:
        phase = 2*np.pi*np.random.rand(1,numFreq)[0]


#irregwavespectrum
freq = wF/(2*np.pi)
Tp = wT
Hs = wH
data = np.conj(np.transpose(np.loadtxt(spectrumDataFile)))
freq_data = data[0]
S_data = data[1]
S_f = np.array([x for x,i in zip(S_data,freq_data) 
                if i>=min(bemFreq)/2/np.pi and i<=max(bemFreq)/2/np.pi])       # Wave Spectrum [m^2-s] for 'EqualEnergy'
wS = S_f/(2*np.pi)                                      # Wave Spectrum [m^2-s/rad] for 'Traditional'
print('\t"spectrumImport" uses the number of imported wave frequencies (not "Traditional" or "EqualEnergy")\n')
# Power per Unit Wave Crest
wk= wF**2/g
if deepWaterWave == 0:
    wk = np.zeros(100)
    for i in range(100):
        wk[i] = wF**2/g/np.tanh(wk*waterDepth)                                        # Calculate Wave Number for Larger Number of Frequencies Before Down Sampling in Equal Energy Method
if deepWaterWave == 1:
    # Deepwater Approximation
    Pw = np.sum(1/2*rho*g**(2)*S_f*dw/wF)
else:
    # Full Wave Power Equation
    Pw = np.sum((1/2)*rho*g*S_f*dw*np.sqrt(9.81/wk*np.tanh(wk*waterDepth))*(1 + 2*wk*waterDepth/np.sinh(2*wk*waterDepth)))
wA =2*wS
wk= wF**2/g
if deepWaterWave == 0:
    for i in range(1,101):
        wk = wF**2/g/np.tanh(wk*waterDepth)

waveAmpTime = np.zeros((maxIt+1,2))
waveAmpTime1 = np.zeros((maxIt+1,2))
waveAmpTime2 = np.zeros((maxIt+1,2))
waveAmpTime3 = np.zeros((maxIt+1,2))
maxRampIT=int(np.round(rampTime/dt))
if rampTime == 0:
    for i in range(maxIt+1):
        for idir in range(np.size(waveDir)):
            t       = (i)*dt
            tmp     = np.sqrt(wA*dw*wavepSread[idir])
            tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(wF*t + phase[:idir])))
            tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge1loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[:idir])))
            tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge2loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[:idir])))
            tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge3loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[:idir])))
            waveAmpTime[i][0]    = t
            waveAmpTime[i][1]    = waveAmpTime[i][1] + sum(tmp1)
            waveAmpTime1[i][0]   = t
            waveAmpTime1[i][1]   = waveAmpTime1[i][1] + sum(tmp11)
            waveAmpTime2[i][0]   = t
            waveAmpTime2[i][1]   = waveAmpTime2[i][1] + sum(tmp12)
            waveAmpTime3[i][0]   = t
            waveAmpTime3[i][1]   = waveAmpTime3[i][1] + sum(tmp13)
else:
    for i in range(maxRampIT):
        for idir in range(np.size(waveDir)):
            t       = (i)*dt
            tmp     = np.sqrt(wA*np.array(dw)*waveSpread[idir])
            tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(wF*t + phase)))
            tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge1loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge2loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge3loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            waveAmpTime[i][0]    = t
            waveAmpTime[i][1]    = waveAmpTime[i][1] + sum(tmp1)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
            waveAmpTime1[i][0]   = t
            waveAmpTime1[i][1]   = waveAmpTime1[i][1] + sum(tmp11)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
            waveAmpTime2[i][0]   = t
            waveAmpTime2[i][1]   = waveAmpTime2[i][1] + sum(tmp12)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
            waveAmpTime3[i][0]   = t
            waveAmpTime3[i][1]   = waveAmpTime3[i][1] + sum(tmp13)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
    for i in range(maxRampIT, maxIt+1):
        for idir in range (np.size(waveDir)):
            t       = (i)*dt
            tmp     = np.sqrt(wA*dw*waveSpread[idir])
            tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(wF*t + phase)))
            tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge1loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge2loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge3loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            waveAmpTime[i][0]    = t
            waveAmpTime[i][1]    = waveAmpTime[i][1] + sum(tmp1)
            waveAmpTime1[i][0]   = t
            waveAmpTime1[i][1]   = waveAmpTime1[i][1] + sum(tmp11)
            waveAmpTime2[i][0]   = t
            waveAmpTime2[i][1]   = waveAmpTime2[i][1] + sum(tmp12)
            waveAmpTime3[i][0]   = t
            waveAmpTime3[i][1]   = waveAmpTime3[i][1] + sum(tmp13)

                
                
                