# -*- coding: utf-8 -*-
"""
WEC toy problem
Created on Wed Feb 20 10:11:42 2020

@author: SungJun Won
"""

import numpy as np 
import matplotlib.pyplot as plt 
  
def initialCondition():
    """
    Return integer vlaues of the
    """
    wHs = 1.1#significant wave hight - unit: m
    wTp = 5.7#Wave Watch III hindcast peak period - unit : s
    return wHs, wTp

def bretschneiderParameters():
    """
    Return integer value of the Bretschneider parameters
    Type I Sea State (developing wind sea)
    """
    alpha = 1#coefficient for energy period
    gamma = 1#spectral peakedness parameter
    n = 5#spectral width parameter
    return alpha,gamma,n

def pierson_moskowitzParameters(u10,wTp,wf,kb):
    """
    incomplete I need to work more on it.
    Return list of integer value of the Pierson-Moskowitz parameters
    Type I&II Sea State 
    (developing wind sea)&(swells, decaying and fully developed wind sea)
    u10: WaveWatch III hindcast wind speed at 10m above sea level
    wf: WaveWatch III hindcast wind fraction of the partition
    kb: Shape coefficient
    """
    alpha = 0.86#coefficient for energy period
    g = 9.8067#gravity - unit:m/s^2
    wTpf = []#peak wave period of fully developed sea
    gamma = []
    n = []
    for i in range(0,len(u10)):
       wTpf.append((2*np.pi*u10[i]*1.95**(1/7))/(0.87*g))
    for i in range(0,len(wTpf)):
        if wTp < wTpf[i]:
            gamma.append(3.3)#spectral peakedness parameter
            n.append(5)#spectral width parameter
        else:
            gamma.append(1)#spectral peakedness parameter
            n.append(5*wf+kb+wTp*(1-wf))
    return alpha,gamma,n

def jonswapParameters():
    """
    Return integer value of the JOHNSWAP parameters
    Type I Sea State (developing wind sea)
    """
    alpha = 0.9#coefficient for energy period
    gamma = 3.3#spectral peakedness parameter
    n = 5#spectral width parameter
    return alpha,gamma,n

def frequency():
    """
    Return integer list of frequency Hz
    """
    f = np.linspace(0.001,0.45,450)
    return f

def gammaSpectrum(wHs,wTp,f,n):
    """
    Return list of float spectral values
    """
    a = n*wHs**2/wTp**4
    b = n/wTp**4
    wFp = 1/wTp#peak frequency
    s = []
    sigma = []
    for i in range(0,len(f)):
        if f[i] <= wFp:
            sigma.append(0.07)
        else:
            sigma.append(0.09)
    for i in range(0,len(f)):
       s.append((a/f[i]**n)*np.exp(-b/f[i]**(n - 1))*
                gamma**np.exp(-(f[i]-wFp)**2/(sigma[i]**2*wFp**2)))
    return s

def waveEnergyPeriod(alpha,wTp):
    """
    Return integer value of wave energy period
    """
    wTpe = alpha*wTp
    return wTpe

def wavePowerDensity(wTpe,wHs):
    """
    Return integer value of wave energy density in deep water
    """
    p0 = 490*wTpe*wHs**2
    return p0

def plot(f,s):
    """
    Return plot of the frequency vs Spectrum for wave and response
    """
    plt.plot(f, s, color = 'red', marker = "o")
    plt.xlabel("Frequency (Hz)") 
    plt.ylabel("Spectrum (m^2/Hz)") 
    plt.show() 

wHs, wTp = initialCondition()
alpha,gamma,n = jonswapParameters()
f = frequency()
s = gammaSpectrum(wHs,wTp,f,n)
wTpe = waveEnergyPeriod(alpha,wTp)
p0 = wavePowerDensity(wTpe,wHs)
plot = plot(f,s)
print(np.mean(np.nan_to_num(s)))
print("Power", p0, "kW/m")
