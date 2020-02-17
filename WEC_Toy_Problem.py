# -*- coding: utf-8 -*-
"""
WEC toy problem
Created on Wed Feb 12 12:42:51 2020

@author: SungJun Won
"""

import numpy as np 
import matplotlib.pyplot as plt 
  
def initialCondition():
    """
    Return integer vlaues of the equivalent mass, equivalent stiffness, and 
    damping ratio.
    """
    meq = 100#units (kg)
    keq = 20#units (N/m)
    zeta = 0.1#Damping Ratio
    r = 1#Frequency Ratio
    return meq,keq,zeta,r
    
def drivingFrequency():
    """
    Return integer value of the sinusoidal force amplitude.
    """
    fo =3#Force Amplitude
    return fo

def naturalFrequency(meq,keq):
    """
    Return integer value of the natural frequency.
    """
    wn = np.sqrt(keq/meq)
    return wn

def excitationFrequency(r,wn):
    """
    Return integer value of the excitation frequency.
    """
    w = r*wn
    return w

def period(wn):
    """
    Return integer value of the period.
    """
    T = 1/(wn/np.pi)
    return T

def criticalDamping(meq,wn):
    """
    Return integer value of the citical damping.
    """
    cc = 2*meq*wn
    return cc

def equivalentDamping(zeta,cc):
    """
    Return integer value of the equivalent damping.
    """
    ceq = zeta*cc
    return ceq

def amplitude(fo,keq,zeta,r):
    """
    Return integer value of the amplitude.
    """
    a = (fo/keq)/np.sqrt((1-r**2)**2 +(2*zeta*r)**2)
    return a

def phase(zeta,r):
    """
    Return integer value of the phase in radian.
    """
    if r == 1:
        theta = np.pi/2
    else:
        theta = np.arctan((2*zeta*r)/(1-r**2))
    return theta

def time(w):
    """
    Return list of float value of time and float value of end time of n cycle
    """
    cycle = 3
    end = 2*np.pi/w*cycle
    t = [i for i in np.linspace(0,end,60)]
    return t, end


def particularSolution(w,a,t,theta):
    """
    Return list of float value of the particular solution.
    """
    xp = []
    for i in range(0,len(t)):
        xp.append(a*np.sin(w*t[i]-theta))
    return xp

def response(fo,keq,xp):
    """
    Return list of float value of the response force.
    """
    d = fo/keq
    r = [i/d for i in xp]
    return r

def force(fo,w,t):
    """
    Return list of float value of the force.
    """
    f = []
    for i in range(0,len(t)):
        f.append(fo*np.sin(w*t[i]))
    return f
    
def power(fo,w,t,end,theta):
    """
    Return integer value of the average power per hour.
    """
    p = []
    for i in range(0,len(t)):
        p.append((fo**2)*w*np.sin(w*t[i])*np.cos(w*t[i]-theta))
    avgP = np.mean(p)*(3600/end)
    return avgP

def plot(t,r,f):
    """
    Return plot of the displacement vs time for wave and response
    """
    plt.plot(t, r, color = 'red', marker = "o")
    plt.plot(t, f, color = 'blue', marker = "o")
    plt.xlabel("Time (s)") 
    plt.ylabel("Displacement (m)") 
    plt.legend(['Response', 'Wave'], loc='upper left')
    plt.show() 
    


meq,keq,zeta,r = initialCondition()
fo = drivingFrequency()
wn = naturalFrequency(meq,keq)
w = excitationFrequency(r,wn)
T = period(wn)
cc = criticalDamping(meq,wn)
ceq = equivalentDamping(zeta,cc)
a = amplitude(fo,keq,zeta,r)
theta = phase(zeta,r)
t,end = time(w)
xp = particularSolution(w,a,t,theta)
r = response(fo,keq,xp)
f = force(fo,w,t)
avgP = power(fo,w,t,end,theta)
x = plot(t,r,f)
print("Average Power per hour:",avgP,"W")