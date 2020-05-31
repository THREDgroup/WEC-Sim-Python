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
waveDir = [0]
waveSpread = [1]
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



WFQSt=np.min(bemFreq)
WFQEd=np.max(bemFreq)
numFreq_interp = 500000
wF = np.arange(WFQSt,WFQEd,(WFQEd-WFQSt)/numFreq_interp) 
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
    phase = 2*np.pi*np.random.rand(1,numFreq)[0]
elif (freqDisc == 'Imported'):
    data = pd.read_csv(spectrumDataFile, sep='\t', header=None)
    if np.size(data,2) == 3:
        freq_data = data[0]
        freq_loc = [1 if i>=min(bemFreq)/2/np.pi and i<=max(bemFreq)/2/np.pi *0 else 0 for i in freq_data]
        phase_data = data[freq_loc][3]
        phase = np.conj(np.transpose(phase_data))
    else:
        phase = 2*np.pi*np.random.rand(1,numFreq)[0]
elapsed2 = time.time() - t11
print(elapsed2)



freq = wF/(2*np.pi)
Tp = wT
Hs = wH
if spectrumType == 'PM':
    # Pierson-Moskowitz Spectrum from Tucker and Pitt (2001)
    B_PM = (5/4)*(1/Tp)**(4)
    A_PM = 0.0081*g**2*(2*np.pi)^(-4)
    S_f  = (A_PM*freq**(-5)*exp(-B_PM*freq**(-4)))              # Wave Spectrum [m^2-s] for 'EqualEnergy'
    wS = S_f/(2*np.pi)                                      # Wave Spectrum [m^2-s/rad] for 'Traditional'
    S_f = wS*2*np.pi
elif spectrumType == 'BS':
    # Bretschneider Sprectrum from Tucker and Pitt (2001)
    B_BS = (1.057/Tp)**4
    A_BS = B_BS*(Hs/2)**2
    S_f = (A_BS*freq**(-5)*np.exp(-B_BS*freq**(-4)))               # Wave Spectrum [m^2-s]
    wS = S_f/(2*np.pi)                                      # Wave Spectrum [m^2-s/rad]
elif spectrumType == 'JS':
    # JONSWAP Spectrum from Hasselmann et. al (1973)
    fp = 1/Tp
    siga = 0.07
    sigb = 0.09                                                 # cutoff frequencies for gamma function
    lind = np.argwhere(np.array(freq)<=fp)
    hind = np.argwhere(np.array(freq)>fp)
    Gf = np.zeros((np.size(freq)))
    Gf[lind] = gamma**np.exp(-(freq(lind)-fp)**2/(2*siga**2*fp**2))
    Gf[hind] = gamma**np.exp(-(freq(hind)-fp)**2/(2*sigb**2*fp**2))
    S_temp = g**2*(2*np.pi)**(-4)*freq**(-5)*np.exp(-(5/4)*(freq/fp)**(-4))
    alpha_JS = Hs**(2)/16/np.trapz(freq,S_temp*Gf)
    S_f = alpha_JS*S_temp*Gf                                 # Wave Spectrum [m^2-s]
    wS = S_f/(2*pi)                                       # Wave Spectrum [m^2-s/rad]
elif spectrumType == 'spectrumImport':
    # Imported Wave Spectrum
    data = pd.read_csv(spectrumDataFile, sep='\t', header=None)
    freq_data = data[0]
    S_data = data[1]
    freq_loc = [1 if i>=min(bemFreq)/2/np.pi and i<=max(bemFreq)/2/np.pi *0 else 0 for i in freq_data]
    S_f = S_data[freq_loc-1]                                    # Wave Spectrum [m^2-s] for 'EqualEnergy'
    wS = S_f/(2*np.pi)                                      # Wave Spectrum [m^2-s/rad] for 'Traditional'
    print('\t"spectrumImport" uses the number of imported wave frequencies (not "Traditional" or "EqualEnergy")\n')
# Power per Unit Wave Crest
wk= wF**2/g
if deepWaterWave == 0:
    for i in range(1,101):
        wk = wF**2/g/np.tanh(wk*waterDepth)                                         # Calculate Wave Number for Larger Number of Frequencies Before Down Sampling in Equal Energy Method
if deepWaterWave == 1:
    # Deepwater Approximation
    Pw = np.sum(1/2*rho*g**(2)*S_f*dw/wF)
else:
    # Full Wave Power Equation
    Pw = np.sum((1/2)*rho*g*S_f*dw*np.sqrt(9.81/wk*np.tanh(wk*waterDepth))*(1 + 2*k*waterDepth/np.sinh(2*k*waterDepth)))
if freqDisc == 'EqualEnergy':
    m0 = np.trapz(np.abs(S_f),freq)
    numBins = numFreq+1
    a_targ = m0/numBins
    SF = np.insert(integrate.cumtrapz(S_f,freq),0,0)
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
            jj = jj + 1
            if wn[kk]+jj >= len(S_f):
                break
        a_targ_real.append(np.min(np.abs([x-(kk+1)*a_targ for x in tmpa[kk]])))
        wna.append(np.argmin(np.abs([x-(kk+1)*a_targ for x in tmpa[kk]])))
        wn.append(wna[kk]+wn[kk])
        a_bins.append(np.trapz(np.abs(S_f[np.arange(wn[kk],wn[kk]+1)-1]),freq[np.arange(wn[kk],wn[kk]+1)-1]))   
    wF = 2*np.pi*freq[wn[1:len(wn)-1]]
    dw = np.append(wF[0]-2*np.pi*freq[wn[0]-1],np.diff(wF)) 
    wS = wS[wn[1:len(wn)-1]]                         # Wave Spectrum [m^2-s/rad] 
wA = 2 * wS                                             # Wave Amplitude [m]
elapsed3 = time.time() - t11
print(elapsed3)

wk= wF**2/g
if deepWaterWave == 0:
    for i in range(1,101):
        wk = wF**2/g/np.tanh(wk*waterDepth)
elapsed4 = time.time() - t11
print(elapsed4)


waveAmpTime = [[],[]]
#waveAmpTime1 = np.zeros((maxIt+1,2))
#waveAmpTime2 = np.zeros((maxIt+1,2))
#waveAmpTime3 = np.zeros((maxIt+1,2))
maxRampIT=int(np.round(rampTime/dt))
if rampTime == 0:
    for i in range(maxIt+1):
        for idir in range(np.size(waveDir)):
            t       = np.arange(maxIt+1)*dt
            tmp     = np.sqrt(wA*dw*wavepSread[idir])
            tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(wF*t + phase)))
            #tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge1loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[:idir])))
            #tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge2loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[:idir])))
            #tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge3loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase[:idir])))
            waveAmpTime[i][0]    = t
            waveAmpTime[i][1]    = waveAmpTime[i][1] + sum(tmp1)
            #waveAmpTime1[i][0]   = t
            #waveAmpTime1[i][1]   = waveAmpTime1[i][1] + sum(tmp11)
            #waveAmpTime2[i][0]   = t
            #waveAmpTime2[i][1]   = waveAmpTime2[i][1] + sum(tmp12)
            #waveAmpTime3[i][0]   = t
            #waveAmpTime3[i][1]   = waveAmpTime3[i][1] + sum(tmp13)
else:
    for i in range(maxRampIT):
        for idir in range(np.size(waveDir)):
            t       = (i)*dt
            tmp     = np.sqrt(wA*np.array(dw)*waveSpread[idir])
            tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(wF*t + phase)))
            #tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge1loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            #tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge2loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            #tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge3loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            waveAmpTime[i][0]    = t
            waveAmpTime[i][1]    = waveAmpTime[i][1] + sum(tmp1)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
            #waveAmpTime1[i][0]   = t
            #waveAmpTime1[i][1]   = waveAmpTime1[i][1] + sum(tmp11)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
            #waveAmpTime2[i][0]   = t
            #waveAmpTime2[i][1]   = waveAmpTime2[i][1] + sum(tmp12)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
            #waveAmpTime3[i][0]   = t
            #waveAmpTime3[i][1]   = waveAmpTime3[i][1] + sum(tmp13)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
    for i in range(maxRampIT, maxIt+1):
        for idir in range (np.size(waveDir)):
            t       = (i)*dt
            tmp     = np.sqrt(wA*dw*waveSpread[idir])
            tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(wF*t + phase)))
            #tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge1loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge1loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            #tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge2loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge2loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            #tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(wF*t - wk*(wavegauge3loc[0]*np.cos(waveDir[idir]*np.pi/180) + wavegauge3loc[1]*np.sin(waveDir[idir]*np.pi/180)) + phase)))
            waveAmpTime[i][0]    = t
            waveAmpTime[i][1]    = waveAmpTime[i][1] + sum(tmp1)
            #waveAmpTime1[i][0]   = t
            #waveAmpTime1[i][1]   = waveAmpTime1[i][1] + sum(tmp11)
            #waveAmpTime2[i][0]   = t
            #waveAmpTime2[i][1]   = waveAmpTime2[i][1] + sum(tmp12)
            #waveAmpTime3[i][0]   = t
            #waveAmpTime3[i][1]   = waveAmpTime3[i][1] + sum(tmp13)
elapsed5 = time.time() - t11
print(elapsed5)

   