#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 23:21:05 2020

@author: logical
"""

from scipy import integrate
import numpy as np

import time
t11 = time.time()
wType = 'irregular'
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
numFreq = 500
phaseSeed = 1
wT = 8
wH = 2.5
spectrumType = 'BS'
waveDir = [30,60]
waveSpread = [30,60]
wavegauge1loc = [0,0]
wavegauge2loc = [0,0]
wavegauge3loc = [0,0]

if wDepth == 'infinite':
    deepWaterWave = 1
    waterDepth = 200
    print('\tInfinite water depth specified in BEM, "waves.waterDepth" set to 200m for vizualisation.\n')
else:
    deepWaterWave = 0
    waterDepth = wDepth

WFQSt=np.min(bemFreq)
WFQEd=np.max(bemFreq)
numFreq_interp = 500000
wF = np.arange(WFQSt,WFQEd+((WFQEd-WFQSt)/numFreq_interp),(WFQEd-WFQSt)/numFreq_interp)
wF1 = wF



dw = np.mean(np.diff(wF))
if np.size(numFreq) == 0:
    numFreq = 500
elapsed = time.time() - t11
print(elapsed)


if phaseSeed != 0:
    np.random.seed(phaseSeed) #Phase seed = 1,2,3,...,etc
else:
    np.random.seed(np.random.shuffle(phaseSeed))
if (freqDisc == 'EqualEnergy') or (freqDisc == 'Traditional'):
    phase = 2*np.pi*np.conj(np.transpose(np.random.rand(numFreq,np.size(waveDir))))
#np.conj(np.transpose(
elapsed2 = time.time() - t11
print(elapsed2)


freq = wF/(2*np.pi)
Tp = wT
Hs = wH
if spectrumType == 'BS':
    # Bretschneider Sprectrum from Tucker and Pitt (2001)
    B_BS = (1.057/Tp)**4
    A_BS = B_BS*(Hs/2)**2
    S_f = (A_BS*freq**(-5)*np.exp(-B_BS*freq**(-4)))               # Wave Spectrum [m^2-s]
    wS = S_f/(2*np.pi)                                      # Wave Spectrum [m^2-s/rad]
wk= wF**2/g
if deepWaterWave == 0:
    for i in range(1,101):
        wk = wF**2/g/np.tanh(wk*waterDepth)                                         # Calculate Wave Number for Larger Number of Frequencies Before Down Sampling in Equal Energy Method
if deepWaterWave == 1:
    # Deepwater Approximation
    Pw = np.sum(1/2*rho*g**(2)*S_f*dw/wF)
else:
    # Full Wave Power Equation
    Pw = np.sum((1/2)*rho*g*S_f*dw*np.sqrt(9.81/wk*np.tanh(wk*waterDepth))*(1 + 2*wk*waterDepth/np.sinh(2*wk*waterDepth)))
elapsed3 = time.time() - t11
print(elapsed3)
if freqDisc == 'EqualEnergy':
    m0 = np.trapz(np.abs(S_f),freq)
    numBins = numFreq+1
    a_targ = m0/numBins
    SF = integrate.cumtrapz(S_f,freq,initial = 0)
    wn = [1]
    tmpa = {}
    a_targ_real = []
    wna = []
    a_bins = []
    for kk in range(numBins):
        jj = 1
        tmpa.update({kk:[0]})
        while (tmpa[kk][jj-1]-(kk+1)*a_targ) < 0:
            tmpa[kk].append(SF[wn[kk]+jj-1])
            jj += 1
            if wn[kk]+jj >= len(S_f):
                break
        a_targ_value = np.abs([x-(kk+1)*a_targ for x in tmpa[kk]])
        a_targ_real.append(np.min(a_targ_value))
        wna.append(np.argmin(a_targ_value))
        #a_targ_real.append(np.min(np.abs([x-(kk+1)*a_targ for x in tmpa[kk]])))
        #wna.append(np.argmin(np.abs([x-(kk+1)*a_targ for x in tmpa[kk]])))
        wn.append(wna[kk]+wn[kk])
        a_bins.append(np.trapz(np.abs(S_f[np.arange(wn[kk],wn[kk]+1)-1]),freq[np.arange(wn[kk],wn[kk]+1)-1]))   
    wF = 2*np.pi*freq[wn[1:len(wn)-1]]
    dw = np.append(wF[0]-2*np.pi*freq[wn[0]-1],np.diff(wF)) 
    wS = wS[wn[1:len(wn)-1]]                         # Wave Spectrum [m^2-s/rad] 
wA = 2 * wS 
elapsed4 = time.time() - t11
print(elapsed4)

wk= wF**2/g
if deepWaterWave == 0:
    for i in range(1,101):
        wk = wF**2/g/np.tanh(wk*waterDepth)
elapsed5 = time.time() - t11
print(elapsed5)


timestep = np.arange(maxIt+1)*dt
initialize = np.zeros((maxIt+1))
waveAmpTime = [timestep,initialize]
waveAmpTime1 = [timestep,initialize]
waveAmpTime2 = [timestep,initialize]
waveAmpTime3 = [timestep,initialize]
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


elapsed6 = time.time() - t11
print(elapsed6)
#result8 = np.conj(np.transpose(np.loadtxt('./freq.txt')))
#np.testing.assert_allclose(freq, result8)
#result7 = np.conj(np.transpose(np.loadtxt('./S_f.txt')))
#np.testing.assert_allclose(S_f, result7)
#result6 = np.conj(np.transpose(np.loadtxt('./wF1.txt')))
#np.testing.assert_allclose(wF1, result6)
#result5 = np.conj(np.transpose(np.loadtxt('./SF.txt')))
#np.testing.assert_allclose(SF, result5)
#result4 = np.conj(np.transpose(np.loadtxt('./wF.txt')))
#np.testing.assert_allclose(wF, result4)
#result3 = np.conj(np.transpose(np.loadtxt('./dw.txt')))
#np.testing.assert_allclose(dw, result3)
#result2 = np.conj(np.transpose(np.loadtxt('./wA.txt')))
#np.testing.assert_allclose(wA, result2)
#result1 = np.conj(np.transpose(np.loadtxt('./phase.txt')))
#np.testing.assert_allclose(phase, result1)
result =  np.conj(np.transpose(np.loadtxt('./waveAmpTime.txt')))
np.testing.assert_allclose(waveAmpTime, result)