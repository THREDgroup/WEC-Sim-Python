# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 07:52:07 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

Note: readData method is included to import data type of both .txt and .mat

Note: irregularWaveSpectrum method has been modified for the faster computation.
    The original method from WEC-sim is commented.

Note: MATLAB column vectors are changed to array for faster computation

Note: waveElevationGrid method and write_paraview_vtp_wave method is moved to
    paraviewClass.py 

Note: Values are equal to tolerance rtol=1e-07, atol=0
    Values will not be exact to WEC-sim due to difference in significant figures 
    or the way some methods are used.
    (e.g) 
    integrate.cumtrapz(S_f,freq) and MATLAB cumtrapz(freq,S_f) does not produce
    same results but silimar values to rtol=1e-06

Note: "RuntimeWarning: overflow encountered in sinh"
    When using irregular wave or spectrumImport, Equal Energy with wDepth value 
    over 100 can generate "RuntimeWarning: overflow encountered in sinh" as the 
    value of sinh in waveSetup method reaches infinity.

"""
from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib 
import warnings
import scipy.io as sio

class WaveClass:
    def inputProperties(self):
        """
        Set up necessary properties of the wave condition based on input file
    
        Returns 
        -------
        wType, T, H, spectrumType, gamma, phaseSeed, spectrumDataFile,
        etaDataFile, numFreq, waveDir, waveSpread, viz, statisticsDataLoad,
        freqDisc, wavegauge1loc, wavegauge2loc, wavegauge3loc, currentSpeed, 
        currentDirection, currentOption, currentDepth 
    
        """
        
        #wType - String defining the type of waves to be generated.
        #Can be one of: 'noWave', 'noWaveCIC', 'regular', 'regularCIC','irregular',
        #'spectrumImport', and 'etaImport' 
        #(Default = 'NOT DEFINED').
        self.wType = 'NOT DEFINED'
        
        #T - [s] Wave period, peak wave period or BEM period.
        #Wave period (regular waves), peak period (irregular waves), or period of 
        #BEM data used for hydrodynamic coefficients ('noWave') 
        #(Default = 'NOT DEFINED').     
        self.T = 'NOT DEFINED'
        
        #H - % [m] Wave height or significant wave height for irregular
        #Wave height (regular waves) or significant wave height (irregular waves) 
        #(Default = 'NOT DEFINED').   
        self.H = 'NOT DEFINED'
        
        #spectrumType -  String containing the wave spectrum type
        #Can be one of : 'PM', 'BS', and 'JS' 
        #(Default = 'NOT DEFINED'). 
        self.spectrumType = 'NOT DEFINED'
        
        #gamma - Only used for 'JS' spectrum type to define gamma 
        #(Default = 3.3)
        self.gamma = 3.3
        
        #phaseSeed - Only used for irregular waves 
        #if equal to 1,2,3,...,etc, the waves phase is seeded.
        #(Default = 0)
        self.phaseSeed = 0
        
        #spectrumDataFile - Data file that contains the spectrum data file
        #(Default = 'NOT DEFINED')
        self.spectrumDataFile = 'NOT DEFINED'
        
        #etaDataFile - Data file that contains the times-series data file
        #(Default = 'NOT DEFINED')
        self.etaDataFile = 'NOT DEFINED'
        
        #freqRange - Min and max frequency for irregular waves. 
        #array with two values, rad/s, (default = frequency range in BEM data)
        #(Default = [])
        #eg. [0.3,0.5]
        self.freqRange = []
        
        #numFreq - of interpolated wave frequencies 
        #Number of frequencies used, varies depending on method:
        #Traditional = 1000, EqualEnergy = 500 or 'Imported'
        #(Default = 0)
        self.numFreq = 0
        
        #waveDir - [deg] Incident wave direction(s)
        #Should be defined as a column vector for more than one wave direction 
        #For multiple wave direction it will be an array of angle in degree
        #[20,45,60,80]
        #(Default = 0)
        self.waveDir = [0]
        
        #waveSpread - Wave Spread probability associated with wave direction(s)
        #Should be defined as a column vector for more than one wave direction
        #For multiple wave direction it will be an array that corresponds to 
        #wave direction [20,45,60,80]
        #(Default = [1])
        self.waveSpread = [1]
        
        #viz - Dictionary defining visualization options
        #Should be a dictionary containing the fields 'numPointsX' and 'numPointsY'. 
        #numPointsX is the number of visualization points in x direction, and 
        #numPointsY the number of visualization points in y direction
        self.viz = {'numPointsX': 50,'numPointsY': 50 }
        
        #statisticsDataLoad - File name from which to load wave statistics data
        #(Default = [])
        self.statisticsDataLoad = []
        
        #freqDisc - Method of frequency discretization for irregular waves. 
        #Options for this variable are 'EqualEnergy' or 'Traditional'.
        #(Default = 'EqualEnergy').
        self.freqDisc = 'EqualEnergy'
     
        #wavegauge1loc - [m] Wave gauge 1 [x,y] location
        #(Default = [0,0]).
        self.wavegauge1loc = [0,0]
        
        #wavegauge2loc - [m] Wave gauge 2 [x,y] location
        #(Default = [0,0]).
        self.wavegauge2loc = [0,0]
    
        #wavegauge3loc - [m] Wave gauge 3 [x,y] location
        #(Default = [0,0]).
        self.wavegauge3loc = [0,0]
        
        #currentSpeed - [m/s] Surface current speed that is uniform along the 
        #water column.
        #(Default = 0).
        self.currentSpeed = 0
      
        #currentDirection - [deg] Surface current direction.
        #(Default = 0).
        self.currentDirection = 0
        
        #currentOption - [-] Define the sub-surface current model to be used in 
        #WEC-Sim.
        #(Default = 0)
        #0 : Depth-independent model
        #1 : 1/7 power law variation with depth
        #2 : linear variation with depth
        #3 : no current
        self.currentOption = 3
        
        #currentDepth - [m] Define the depth over which the sub-surface current is 
        #modeled.
        #For options (1) and (2) the currentDepth must be defined. The
        #current is not calculated for any depths greater than the
        #specified currentDepth. (Default = 0).
        self.currentDepth = 0
    
    def _internalProperties(self):
        """
        The following properties are for internal use
    
        Returns
        -------
        typeNum, bemFreq, waterDepth, deepWaterWave, waveAmpTime,waveAmpTime1, 
        waveAmpTime2, waveAmpTime3, A, w, phase, dw, k, S, Pw
    
        """
        
        #typeNum - Number to represent different type of waves
        self.typeNum = []
        
        #bemFreq - Number of wave frequencies from BEM
        self.bemFreq = []
        
        #waterDepth - [m] Water depth (from BEM)
        self.waterDepth = []
        
        #deepWaterWave - Deep water or not, depending on input from WAMIT, 
        #NEMOH and AQWA
        self.deepWaterWave = []
        
        #waveAmpTime - [m] Wave elevation time history
        self.waveAmpTime = []
        
        #waveAmpTime1 - [m] Wave elevation time history at a wave gauge 1 location 
        #specified by user
        self.waveAmpTime1 = []
        
        #waveAmpTime2 - [m] Wave elevation time history at a wave gauge 2 location 
        #specified by user
        self.waveAmpTime2 = []
        
        #waveAmpTime3 - [m] Wave elevation time history at a wave gauge 3 location 
        #specified by user
        self.waveAmpTime3 = []
        
        #A - [m] Wave amplitude for regular waves or 2*(wave spectrum vector) for 
        #irregular waves
        self.A = [] 
        
        #w - [rad/s] Wave frequency (regular waves) or 
        #wave frequency vector (irregular waves)
        self.w = []
        
        #phase - [rad] Wave phase (only used for irregular waves)
        self.phase = 0
        
        #dw - [rad] Frequency spacing for irregular waves.
        self.dw = 0 
        
        #k - Wave Number
        self.k = []
        
        #S - Wave Spectrum [m^2-s/rad] for 'Traditional'
        self.S = []
        
        #Pw - Wave Power Per Unit Wave Crest [W/m]
        self.Pw = []
            
    def __init__(self,wType):
        """
        waveClass constructor function
        Parameters
        ----------
        Returns
        -------
        self - waveClass object
    
        """
        
        self.wType = wType
        self.gamma = 3.3
        self.w = []
        
        if self.wType == 'noWave': #No Waves with Constant Hydrodynamic Coefficients
            self.typeNum = 0
        elif self.wType == 'noWaveCIC': #No Waves w/Convolution Integral Calculation
            self.typeNum = 1
        elif self.wType == 'regular': #Regular Waves with Constant Hydrodynamic Coefficients
            self.typeNum = 10
        elif self.wType == 'regularCIC': #Regular Waves w/Convolution Integral Calculation
            self.typeNum = 11
        elif self.wType == 'irregular': #Irregular Waves with 'PM', BS' or 'JS' wave spectrum
            self.typeNum = 20
        elif self.wType == 'spectrumImport': #Irregular waves with imported wave spectrum
            self.typeNum = 21
        elif self.wType == 'etaImport': #Waves with imported wave elevation time-history
            self.typeNum = 30
        
    def plotEta(self, *rampTime):
        """
        Plot wave elevation time-history
        This method is for our visual understanding. 
        Does not get used in any method
        rampTime input argument is optional

        """
        plt.figure(figsize=(10,8))
        plt.plot(self.waveAmpTime[0],self.waveAmpTime[1])
        plt.title('Wave Surfave Elevation')
        if np.size(rampTime) == 1: # If ramptime is given we can specify it but will generate similar graph
            plt.plot(np.array(rampTime),np.array(1.5*min(self.waveAmpTime[1]),1.5*max(self.waveAmpTime[1])))
            plt.title(['Wave Surfave Elevation, Ramp Time ' + str(rampTime) +' (s)'])
        plt.xlabel('Time (s)')
        plt.ylabel('Eta (m)')
    
    def plotSpectrum(self):
        """
        Plot wave spetrum
        This method is for our visual understanding. 
        Does not get used in any method

        """
        m0 = np.trapz(self.S, x = self.w)
        HsTest = 4*np.sqrt(m0)
        I = np.argmax(np.abs(self.S))
        wp = self.w[I]
        TpTest = 2*np.pi/wp
        plt.figure(figsize=(10,8))
        plt.plot(self.w,self.S,'s-')
        plt.plot(wp,np.array(max(self.S)))
        plt.xlim([0, max(self.w)])
        plt.title([self.spectrumType, ' Spectrum, T_p= ' + str(TpTest) + ' [s], H_m_0= ' + str(HsTest) + ', [m]'])
        if self.spectrumType == 'JS':
            plt.title([self.spectrumType, ' Spectrum, T_p= ' + str(TpTest) + ' [s], H_m_0= ' + str(HsTest) + ', [m], gamma = ' + str(self.gamma)])
        
        plt.xlabel('Frequency (rad/s)')
        plt.ylabel('Spectrum (m^2-s/rad)')
        
    def waveSetup(self,bemFreq,wDepth,rampTime,dt,maxIt,g,rho,endTime):
        """
        Set up wave for all wave type. Input parameters are from simulationClass.py
        This method is the most important method to check out in this class

        """
        self.bemFreq = bemFreq  # array of the beam frequency from BEMIO. Only the max and min values gets used.
        self.setWaveProps(wDepth) # method called to set up wave properties.
        if (self.wType == 'noWave') or (self.wType == 'noWaveCIC'):
            if np.size(self.w) == 0:
                self.w = 2*np.pi/self.T
            self.waveNumber(g) # used to set self.k
            self.A = self.H/2
            self.waveElevNowave(maxIt,dt) # method called to set wave elevation
        elif (self.wType == 'regular') or (self.wType == 'regularCIC'):
            if np.size(self.w) == 0:
                self.w = 2*np.pi/self.T
            self.A = self.H/2
            self.waveNumber(g) # used to set self.k
            self.waveElevReg(rampTime, dt, maxIt) # method called to set wave elevation
            self.wavePowerReg(g,rho) # method called to find wave power
        elif (self.wType == 'irregular') or (self.wType == 'spectrumImport'):
            WFQSt = np.min(bemFreq)
            WFQEd = np.max(bemFreq)
            if np.size(self.freqRange) != 0: # frequency range that can be set bu user. eg.[0.2, 0.5]
                if self.freqRange[0] > WFQSt and self.freqRange[0] > 0: # if minimum frequency range value set by user is greater than minimum beam frequency and 0
                    WFQSt = self.freqRange[0]   # set new minimum beam frequency provided by user
                else:
                    warnings.warn("Min frequency range outside BEM data, min frequency set to min BEM frequency",DeprecationWarning)
                if self.freqRange[1] < WFQEd and self.freqRange[1] > WFQSt: # if maximum frequency range value set by user is lower than maximum beam frequency but greater than minimum beam frequency
                    WFQEd = self.freqRange[1]   # set new maximum beam frequncy provided by user
                else:
                    warnings.warn("Max frequency range outside BEM data, max frequency set to max BEM frequency",DeprecationWarning)
            if self.freqDisc == 'Traditional':    # Traditional method of computing. Refer to theory of waveclass provided by WEC-Sim to understand the theory.
                if np.size(self.numFreq) == 0:  # numfreq for Traditional is 1000 for default
                    self.numFreq = 1000
                self.w = np.arange(WFQSt,WFQEd+((WFQEd-WFQSt)/(self.numFreq-1)),(WFQEd-WFQSt)/(self.numFreq-1))
                self.dw = np.ones(shape=(self.numFreq,1))*(WFQEd-WFQSt)/(self.numFreq-1)
            elif self.freqDisc == 'EqualEnergy':    # Default way of computing irregular wave. Refer to theory of waveclass provided by WEC-Sim to understand the theory.
                numFreq_interp = 500000     # number of interpolation that will set array size for SF and S_f used in irregWaveSpectrum method. Lowering this value might decrease the run time but accuracy will decrease
                self.w = np.arange(WFQSt,WFQEd+((WFQEd-WFQSt)/numFreq_interp),(WFQEd-WFQSt)/numFreq_interp)
                self.dw = np.mean(np.diff(self.w))
                if np.size(self.numFreq) == 0:  # numfreq for EqualEnergy is 500 for default
                    self.numFreq = 500
            elif self.freqDisc == 'Imported': # set from setWaveProps method
                data = self.readData(self.spectrumDataFile) # call on readData method to get files in both mat file and txt file
                freq_data = data[0] # the first row out of the three rows in spectrum data file is frequency
                self.w = np.array([i*2*np.pi for i in freq_data 
                                    if i>=min(bemFreq)/2/np.pi and i<=max(bemFreq)/2/np.pi])
                self.numFreq = len(self.w)
                self.dw = np.zeros(self.numFreq)
                self.dw[0]= np.array(self.w[1]-self.w[0])
                self.dw[1:self.numFreq-1]= np.array((self.w[2:]-self.w[:-2])/2)
                self.dw[self.numFreq-1]= np.array(self.w[-1]-self.w[-2])
            self.setWavePhase() # method called to generate wave phase. there can be multiple phase if there are more than one wave direction
            self.irregWaveSpectrum(g,rho) # method called to calculate for different kinds of irregular wave calculation methods
            self.waveNumber(g) # used to set self.k
            self.waveElevIrreg(rampTime, dt, maxIt, self.dw) # method called to set wave elevation
        elif self.wType == 'etaImport':
            #Import 'etaImport' time-series here and interpolate
            data = self.readData(self.etaDataFile) #Import time-series
            t = np.arange(0,endTime+dt,dt)      #simulation time
            self.waveElevUser(rampTime, dt, maxIt, data, t) # method called to set wave elevation
            t2 = np.arange(maxIt+1)*dt
            initialZeros = np.zeros((maxIt+1))
            self.waveAmpTime1 = [t2,initialZeros] # set wave elevation as zeros since we just wants to look at imported data
            self.waveAmpTime2 = self.waveAmpTime1 # set wave elevation as zeros since we just wants to look at imported data
            self.waveAmpTime3 = self.waveAmpTime1 # set wave elevation as zeros since we just wants to look at imported data

    def setWaveProps(self,wDepth):
        """
        Sets wave properties
        check for wave depth
        check for wave type of noWave, noWaveCIC, and spectrumImport

        """
        if wDepth == 'infinite': # Can be 'infinite' or number. From BEMIO
            self.deepWaterWave = 1  # means deep water
            self.waterDepth = 200   # minimu water depth for deep water set
            print('\tInfinite water depth specified in BEM, "waves.waterDepth" set to 200m for vizualisation.\n')
        else:
            self.deepWaterWave = 0  # means shallow water
            self.waterDepth = wDepth # set it to specific water depth. Can cause warning when some calculations reach infinity
        if self.wType == 'noWave':
            self.H = 0     # set wave height as 0
        elif self.wType == 'noWaveCIC':
            self.H = 0     # set wave height as 0
            if np.size(self.w) == 0 and self.T == 'NOT DEFINED':
                self.w = np.min(self.bemFreq)
                self.T = 2*np.pi/self.w
            elif np.size(self.w) == 0:
                self.w = 2*np.pi/self.T
            else:
                self.T = 2*np.pi/self.w
        elif self.wType == 'spectrumImport':
            self.H = 0
            self.T = 0
            self.freqDisc = 'Imported' # one of the type of freqDisc that used later in waveSetup method 
            self.spectrumType = 'spectrumImport'
             
    def waveNumber(self, g):
        """
        Calculate wave number

        """
        self.k= self.w**2/g # general method of calculating k
        if self.deepWaterWave == 0: # calculate k directly from specific water depth
            for i in range(100):
                self.k = self.w**2/g/np.tanh(self.k*self.waterDepth)

    def readData(self,file):
        """
        Return Data from the file
        Supported file type: .txt, .mat
        
        """
        if file.endswith('.txt'):
            data = np.conj(np.transpose(np.loadtxt(file))) # transforms data in to array no matter it was in vector form or array form
        elif file.endswith('.mat'): # specific for MATLAB data file. Allows collaboration between MATLAB user and Python user.
            matFile = sio.loadmat(file) 
            keys = list(matFile.keys())[-1]
            data = np.conj(np.transpose(matFile[keys])) #  this transforms data in to array no matter it was in vector form or array form
        return data          

    def setWavePhase(self):
        """
        Sets the irregular wave's random phase
        MATLAB and Python use same random number generator
        multiple arrays of phase is not supported for regular wave as it was not supported in the WEC-Sim

        """
        if self.phaseSeed != 0:
            np.random.seed(self.phaseSeed) #Phase seed = 1,2,3,...,etc
        else:
            np.random.seed(np.random.shuffle(self.phaseSeed)) # shuffle phase seed
        if (self.freqDisc == 'EqualEnergy') or (self.freqDisc == 'Traditional'): 
            self.phase = 2*np.pi*np.conj(np.transpose(np.random.rand(self.numFreq,np.size(self.waveDir)))) # for multiple wave direction, multiple arrays of phase will be made
        elif (self.freqDisc == 'Imported'):
            data = self.readData(self.spectrumDataFile)
            if len(data) == 3: # if imported spectrum data file is correct it should have 3 rows of data
                freq_data = data[0]
                self.phase = np.array([[x for x,i in zip(data[2],freq_data) 
                                       if i>=min(self.bemFreq)/2/np.pi and i<=max(self.bemFreq)/2/np.pi]])
            else:
                self.phase = 2*np.pi*np.random.rand(1,self.numFreq) # if imported spectrum data is faulty, phase will be calculated randomly
            
    def waveElevNowave(self,maxIt,dt):
        """
        Set noWave elevation time-history

        """
        t = np.arange(maxIt+1)*dt
        initialZeros = np.zeros((maxIt+1))
        self.waveAmpTime = [t,initialZeros]     # since this is for no wave type wave height will be all zeros
        self.waveAmpTime1 = self.waveAmpTime    # since this is for no wave type wave height will be all zeros
        self.waveAmpTime2 = self.waveAmpTime    # since this is for no wave type wave height will be all zeros
        self.waveAmpTime3 = self.waveAmpTime    # since this is for no wave type wave height will be all zeros
        
    def waveElevReg(self, rampTime, dt, maxIt):
        """
        Calculate regular wave elevation time history
        modified so that it would compute waveAmpTime1, waveAmpTime2, 
        waveAmpTime3 only when wave guage location is specified
        
        """
        t = np.arange(maxIt+1)*dt # array of time with dt time steps
        self.waveAmpTime = [t,[]]
        if rampTime == 0:
            c1 = self.w*t
            self.waveAmpTime[1]    = self.A*np.cos(c1)
        else:
            maxRampIT = int(np.round(rampTime/dt))
            t = np.arange(maxRampIT)*dt # array of time with dt time steps until maxRampIT
            t2 = np.arange(maxRampIT,maxIt+1)*dt # array of time with dt time steps from maxRampIT to the end
            c1 = self.w*t
            c2 = self.w*t2
            ad = (1+np.cos(np.pi+np.pi*np.arange(maxRampIT)/maxRampIT))/2
            self.waveAmpTime[1]    = np.append(self.A*np.cos(c1)*ad, self.A*np.cos(c2))
        self.waveAmpTime1 = self.waveAmpTime # if wave guage location is not set, wave elevation is same as waveAmpTime
        self.waveAmpTime2 = self.waveAmpTime # if wave guage location is not set, wave elevation is same as waveAmpTime
        self.waveAmpTime3 = self.waveAmpTime # if wave guage location is not set, wave elevation is same as waveAmpTime
        if self.wavegauge1loc[0] != 0 or self.wavegauge1loc[1] != 0 or self.wavegauge2loc[0] != 0 or self.wavegauge2loc[1] != 0 or self.wavegauge3loc[0] != 0 or self.wavegauge3loc[1] != 0:
            t = np.arange(maxIt+1)*dt # array of time with dt time steps
            self.waveAmpTime1 = [t,[]] # set to empty array of wave elevation. If it is not set, error occurs
            self.waveAmpTime2 = [t,[]] # set to empty array of wave elevation. If it is not set, error occurs
            self.waveAmpTime3 = [t,[]] # set to empty array of wave elevation. If it is not set, error occurs
            if rampTime == 0:
                c1 = self.w*t # multiple of array of frequency and time with dt time steps
                c_cos = np.cos(self.waveDir[0]*np.pi/180)
                c_sin = np.sin(self.waveDir[0]*np.pi/180)
                self.waveAmpTime1[1]   = self.A*np.cos(c1-self.k*(self.wavegauge1loc[0]*c_cos + self.wavegauge1loc[1]*c_sin))
                self.waveAmpTime2[1]   = self.A*np.cos(c1-self.k*(self.wavegauge2loc[0]*c_cos + self.wavegauge2loc[1]*c_sin))
                self.waveAmpTime3[1]   = self.A*np.cos(c1-self.k*(self.wavegauge3loc[0]*c_cos + self.wavegauge3loc[1]*c_sin))
            else:
                c_cos = np.cos(self.waveDir[0]*np.pi/180)
                c_sin = np.sin(self.waveDir[0]*np.pi/180)
                self.waveAmpTime1[1]   = np.append(self.A*np.cos(c1-self.k*(self.wavegauge1loc[0]*c_cos + self.wavegauge1loc[1]*c_sin))*ad, 
                                                   self.A*np.cos(c2-self.k*(self.wavegauge1loc[0]*c_cos + self.wavegauge1loc[1]*c_sin)))
                self.waveAmpTime2[1]   = np.append(self.A*np.cos(c1-self.k*(self.wavegauge2loc[0]*c_cos + self.wavegauge2loc[1]*c_sin))*ad, 
                                                   self.A*np.cos(c2-self.k*(self.wavegauge2loc[0]*c_cos + self.wavegauge2loc[1]*c_sin)))
                self.waveAmpTime3[1]   = np.append(self.A*np.cos(c1-self.k*(self.wavegauge3loc[0]*c_cos + self.wavegauge3loc[1]*c_sin))*ad, 
                                                   self.A*np.cos(c2-self.k*(self.wavegauge3loc[0]*c_cos + self.wavegauge3loc[1]*c_sin)))

    def wavePowerReg(self,g,rho):
        """
        Calculate wave power per unit wave crest for regular waves
        
        """
        if self.deepWaterWave == 1:
            # Deepwater Approximation
            self.Pw = 1/(8*np.pi)*rho*g**(2)*(self.A)**(2)*self.T
        else:
            # Full Wave Power Equation
            self.Pw = rho*g*(self.A)**(2)/4*np.sqrt(g/self.k*np.tanh(self.k*self.waterDepth))*(1+2*self.k*self.waterDepth/np.sinh(self.k*self.waterDepth))
        
    def irregWaveSpectrum(self,g,rho):
        """
        Calculate irregular wave spectrum vector
        Check theory section in WEC-Sim documentation for a reference
        
        """
        freq = self.w/(2*np.pi)
        Tp = self.T
        Hs = self.H
        if self.spectrumType == 'PM':
            # Pierson-Moskowitz Spectrum from Tucker and Pitt (2001)
            B_PM = (5/4)*(1/Tp)**(4)
            A_PM = 0.0081*g**2*(2*np.pi)**(-4)
            S_f  = (A_PM*freq**(-5)*np.exp(-B_PM*freq**(-4)))              # Wave Spectrum [m^2-s] for 'EqualEnergy'
            self.S = S_f/(2*np.pi)                                      # Wave Spectrum [m^2-s/rad] for 'Traditional'
            S_f = self.S*2*np.pi
        elif self.spectrumType == 'BS':
            # Bretschneider Sprectrum from Tucker and Pitt (2001)
            B_BS = (1.057/Tp)**4
            A_BS = B_BS*(Hs/2)**2
            S_f = (A_BS*freq**(-5)*np.exp(-B_BS*freq**(-4)))               # Wave Spectrum [m^2-s]
            self.S = S_f/(2*np.pi)                                      # Wave Spectrum [m^2-s/rad]
        elif self.spectrumType == 'JS':
            # JONSWAP Spectrum from Hasselmann et. al (1973)
            fp = 1/Tp
            siga = 0.07
            sigb = 0.09                                                 # cutoff frequencies for gamma function
            Gf = np.zeros((np.size(freq)))
            for lind in np.argwhere(np.array(freq)<=fp):
                Gf[lind] = self.gamma**np.exp(-(freq[lind]-fp)**2/(2*siga**2*fp**2))
            for hind in np.argwhere(np.array(freq)>fp):
                Gf[hind] = self.gamma**np.exp(-(freq[hind]-fp)**2/(2*sigb**2*fp**2))
            S_temp = g**2*(2*np.pi)**(-4)*freq**(-5)*np.exp(-(5/4)*(freq/fp)**(-4))
            alpha_JS = Hs**(2)/16/np.trapz(S_temp*Gf,freq)
            S_f = alpha_JS*S_temp*Gf                                 # Wave Spectrum [m^2-s]
            self.S = S_f/(2*np.pi)                                       # Wave Spectrum [m^2-s/rad]
        elif self.spectrumType == 'spectrumImport':
            # Imported Wave Spectrum
            data = self.readData(self.spectrumDataFile)
            freq_data = data[0]
            S_data = data[1]
            S_f = np.array([x for x,i in zip(S_data,freq_data) 
                            if i>=min(self.bemFreq)/2/np.pi and i<=max(self.bemFreq)/2/np.pi])       # Wave Spectrum [m^2-s] for 'EqualEnergy'
            self.S = S_f/(2*np.pi)                                      # Wave Spectrum [m^2-s/rad] for 'Traditional'
            print('\t"spectrumImport" uses the number of imported wave frequencies (not "Traditional" or "EqualEnergy")\n')
            # Power per Unit Wave Crest
        self.waveNumber(g)                                          # Calculate Wave Number for Larger Number of Frequencies Before Down Sampling in Equal Energy Method
        if self.deepWaterWave == 1:
            # Deepwater Approximation
            self.Pw = np.sum(1/2*rho*g**(2)*S_f*self.dw/self.w)
        else:
            # Full Wave Power Equation
            self.Pw = np.sum((1/2)*rho*g*S_f*self.dw*np.sqrt(9.81/self.k*np.tanh(self.k*self.waterDepth))*(1 + 2*self.k*self.waterDepth/np.sinh(2*self.k*self.waterDepth)))
        if self.freqDisc == 'EqualEnergy':
            m0 = np.trapz(np.abs(S_f),freq)
            numBins = self.numFreq+1
            a_targ = m0/numBins
            SF = np.insert(integrate.cumtrapz(S_f,freq),0,0) # integrate cumtrapz produce similar but not exact number that MATLAB produces. From my test E-6 decimal were same
            # Modified method for solving wn
            # start of the midified code:
            
            wn = [np.argmin(np.abs(SF-(x+1)*a_targ))+1 for x in range(numBins)]
            self.w = 2*np.pi*freq[wn[:len(wn)-1]]
            self.dw = np.append(self.w[0]-2*np.pi*freq[0],np.diff(self.w)) 
            self.S = self.S[wn[:len(wn)-1]]                         # Wave Spectrum [m^2-s/rad] 
            
            # end of the modified code:
            
            # The original method provided from WEC-Sim
            # tested but commented due to slow computation speed. wn is only required variable for future computation
            # if you want to check each values, comment the modified method uncomment the original method
            
            # start of the original code:
            #
            # wn = [1]
            # tmpa = {}                 
            # a_targ_real = []                                       
            # wna = []
            # a_bins = []
            # for kk in range(numBins):
            #     jj = 0
            #     tmpa.update({kk:[0]})
            #     while (tmpa[kk][jj]-(kk+1)*a_targ) < 0:
            #         tmpa[kk].append(SF[wn[kk]+jj])
            #         jj += 1
            #         if wn[kk]+jj+1 >= len(S_f):
            #             break
            #     wn.append(np.argmin(np.abs(tmpa[kk]-(kk+1)*a_targ))+wn[kk])
            #     a_targ_value = np.abs([x-(kk+1)*a_targ for x in tmpa[kk]])     
            #     a_targ_real.append(np.min(a_targ_value))
            #     wna.append(np.argmin(a_targ_value))
            #     wn.append(wna[kk]+wn[kk])
            #     a_bins.append(np.trapz(np.abs(S_f[np.arange(wn[kk],wn[kk]+1)-1]),freq[np.arange(wn[kk],wn[kk]+1)-1]))   
            # self.w = 2*np.pi*freq[wn[1:len(wn)-1]]
            # self.dw = np.append(self.w[0]-2*np.pi*freq[wn[0]-1],np.diff(self.w)) 
            # self.S = self.S[wn[1:len(wn)-1]]                         # Wave Spectrum [m^2-s/rad] 
            #
            # end of the original code:
                
        self.A = 2 * self.S   # Wave Amplitude [m]
           
    def waveElevIrreg(self,rampTime,dt,maxIt,df):
        """
        Calculate irregular wave elevetaion time history
        modified so that it would compute waveAmpTime1, waveAmpTime2, 
        waveAmpTime3 only when wave guage location is specified
        Unlike waveElevReg, waveElevIrreg can compute wave elevation for
        multiple wave direction. multple wave direction produces multiple 
        phase arrays and uses wave spread value that matches wave direction in 
        each calculation.
        
        """
        t = np.arange(maxIt+1)*dt # array of time with dt time steps
        initialZeros = np.zeros((maxIt+1))
        self.waveAmpTime = [t,initialZeros]
        maxRampIT=int(np.round(rampTime/dt))
        iiter = np.size(self.waveDir)
        tmp = np.sqrt(np.matlib.repmat(self.A,iiter,1)*np.matlib.repmat(df,iiter,1)*np.transpose([self.waveSpread,]))
        c1 = np.matlib.repmat(self.w,iiter,1) # matlib.repmat method that repeats arrays which results in matrix form of arrays
        if rampTime == 0:
            for i in range(maxIt+1): # keeping for loop was faster than changing it to all matrix computation even when there were multiple wave direction
                t       = (i)*dt
                tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(c1*t + self.phase)))  
                self.waveAmpTime[1][i]    = np.sum(tmp1)
        else:
            for i in range(maxRampIT):
                t       = (i)*dt
                tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(c1*t + self.phase)))  
                ad = (1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
                self.waveAmpTime[1][i]    = np.sum(np.sum(tmp1, axis=1)*ad)
            for i in range(maxRampIT, maxIt+1):
                t       = (i)*dt
                tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(c1*t + self.phase)))  
                self.waveAmpTime[1][i]    = np.sum(tmp1)
        self.waveAmpTime1 = self.waveAmpTime # if wave guage location is not set, wave elevation is same as waveAmpTime
        self.waveAmpTime2 = self.waveAmpTime # if wave guage location is not set, wave elevation is same as waveAmpTime
        self.waveAmpTime3 = self.waveAmpTime # if wave guage location is not set, wave elevation is same as waveAmpTime   
        if self.wavegauge1loc[0] != 0 or self.wavegauge1loc[1] != 0 or self.wavegauge2loc[0] != 0 or self.wavegauge2loc[1] != 0 or self.wavegauge3loc[0] != 0 or self.wavegauge3loc[1] != 0:
            c2 = np.matlib.repmat(self.k,iiter,1)
            c_cos = np.cos(np.transpose([self.waveDir,])*np.pi/180)
            c_sin = np.sin(np.transpose([self.waveDir,])*np.pi/180)
            t = np.arange(maxIt+1)*dt # array of time with dt time steps
            self.waveAmpTime1 = [t,np.zeros((maxIt+1))] # set to arrays of zero get an error if it is not set individually
            self.waveAmpTime2 = [t,np.zeros((maxIt+1))] # set to arrays of zero get an error if it is not set individually
            self.waveAmpTime3 = [t,np.zeros((maxIt+1))] # set to arrays of zero get an error if it is not set individually
            if rampTime == 0:
                for i in range(maxIt+1):
                    t       = (i)*dt
                    tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(c1*t - c2*(self.wavegauge1loc[0]*c_cos + self.wavegauge1loc[1]*c_sin) + self.phase)))
                    tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(c1*t - c2*(self.wavegauge2loc[0]*c_cos + self.wavegauge2loc[1]*c_sin) + self.phase)))
                    tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(c1*t - c2*(self.wavegauge3loc[0]*c_cos + self.wavegauge3loc[1]*c_sin) + self.phase)))
                    self.waveAmpTime1[1][i]   = np.sum(tmp11)
                    self.waveAmpTime2[1][i]   = np.sum(tmp12)
                    self.waveAmpTime3[1][i]   = np.sum(tmp13)
            else:
                for i in range(maxRampIT):
                    t       = (i)*dt
                    tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(c1*t - c2*(self.wavegauge1loc[0]*c_cos + self.wavegauge1loc[1]*c_sin) + self.phase)))
                    tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(c1*t - c2*(self.wavegauge2loc[0]*c_cos + self.wavegauge2loc[1]*c_sin) + self.phase)))
                    tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(c1*t - c2*(self.wavegauge3loc[0]*c_cos + self.wavegauge3loc[1]*c_sin) + self.phase)))
                    ad = (1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
                    self.waveAmpTime1[1][i]   = np.sum(np.sum(tmp11, axis=1)*ad)
                    self.waveAmpTime2[1][i]   = np.sum(np.sum(tmp12, axis=1)*ad)
                    self.waveAmpTime3[1][i]   = np.sum(np.sum(tmp13, axis=1)*ad)
                for i in range(maxRampIT, maxIt+1):
                    t       = (i)*dt
                    tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(c1*t - c2*(self.wavegauge1loc[0]*c_cos + self.wavegauge1loc[1]*c_sin) + self.phase)))
                    tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(c1*t - c2*(self.wavegauge2loc[0]*c_cos + self.wavegauge2loc[1]*c_sin) + self.phase)))
                    tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(c1*t - c2*(self.wavegauge3loc[0]*c_cos + self.wavegauge3loc[1]*c_sin) + self.phase)))
                    self.waveAmpTime1[1][i]   = np.sum(tmp11)
                    self.waveAmpTime2[1][i]   = np.sum(tmp12)
                    self.waveAmpTime3[1][i]   = np.sum(tmp13)
    
    def waveElevUser(self,rampTime,dt,maxIt,data,t):
        """
        Calculate imported wave elevation time history
        used in ETAimport wave type
        
        """
        self.waveAmpTime = [[],[]]
        maxRampIT = int(np.round(rampTime/dt))
        data_t = data[0]                    # Data Time [s]
        data_x = data[1]                    # Wave Surface Elevation [m] 
        #t = np.arange(0,endTime+dt,dt)  
        self.waveAmpTime[0] = t
        self.waveAmpTime[1] = np.interp(t,data_t,data_x)
        if rampTime != 0:
            for i in range(maxRampIT+1):
                self.waveAmpTime[0][i] = t[i]
                self.waveAmpTime[1][i] *= (1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
  
    def listInfo(self):
        """
        List wave info

        """
        print('\n Wave Environment: \n')
        if self.wType == 'noWave':
            print('\tWave Type                            = No Wave (Constant Hydrodynamic Coefficients)\n')
            print('\tHydro Data Wave Period, T (sec)    	= {:f}\n'.format(self.T))
        if self.wType == 'regular':
            print('\tWave Type                            = Regular Waves (Constant Hydrodynamic Coefficients)\n')
            print('\tWave Height, H (m)                   = {:f}\n'.format(self.H))
            print('\tWave Period, T (sec)                 = {:f}\n'.format(self.T))
        if self.wType == 'noWaveCIC':
            print('\tWave Type                            = No Wave (Convolution Integral Calculation)\n')
        if self.wType == 'regularCIC':
            print('\tWave Type                            = Regular Waves (Convolution Integral Calculation)\n')
            print('\tWave Height, H (m)                   = {:f}\n'.format(self.H))
            print('\tWave Period, T (sec)                 = {:f}\n'.format(self.T))
        if self.wType == 'irregular':
            if self.phaseSeed == 0:
                print('\tWave Type                            = Irregular Waves (Arbitrary Random Phase)\n')
            else:
                print('\tWave Type                            = Irregular Waves (Predefined Random Phase)\n')
            self.printWaveSpectrumType
            print('\tSignificant Wave Height, Hs      (m) = {:f}\n'.format(self.H))
            if self.spectrumType == 'PM':
                print('\tNOTE: Pierson-Moskowitz does not use Hs to define spectrum\n')
            print('\tPeak Wave Period, Tp           (sec) = {:f}\n'.format(self.T))
        if self.wType == 'spectrumImport':
            if np.size(np.loadtxt(self.spectrumDataFile),1) == 3:
                print('\tWave Type                            = Irregular waves with imported wave spectrum (Imported Phase)\n')
            elif self.phaseSeed == 0:
                print('\tWave Type                            = Irregular waves with imported wave spectrum (Random Phase)\n')
            else:
                print('\tWave Type                            = Irregular waves with imported wave spectrum (Seeded Phase)\n')
            self.printWaveSpectrumType
        if self.wType == 'etaImport':
            print( '\tWave Type                           = Waves with imported wave elevation time-history\n')
            print(['\tWave Elevation Time-Series File    	= ' + self.etaDataFile + '  \n'])
                     
    def checkinputs(self):
        """
        Check user inputs

        """
        # 'noWave' period undefined for hydro data
        if (self.wType == 'noWave') and (self.T == 'NOT DEFINED'):
            warnings.warn('"waves.T" must be defined for the hydrodynamic data period when using the "noWave" wave type',DeprecationWarning)
        # spectrumDataFile defined for 'spectrumImport' case
        if (self.wType == 'spectrumImport') and (self.spectrumDataFile == 'NOT DEFINED'):
            warnings.warn('The spectrumDataFile variable must be defined when using the "spectrumImport" wave type',DeprecationWarning)
        # check waves types
        types = ['noWave', 'noWaveCIC', 'regular', 'regularCIC', 'irregular', 'spectrumImport', 'etaImport']
        if [1 for x in types if x == self.wType] != [1]:
            warnings.warn('Unexpected wave environment type setting, choose from: ' \
                    '"noWave", "noWaveCIC", "regular", "regularCIC", "irregular", "spectrumImport", and "etaImport".',DeprecationWarning)
            
    def printWaveSpectrumType(self):
        """
        Lists the wave spectrum type

        """
        if self.spectrumType == 'BS':
            print('\tSpectrum Type                        = Bretschneider \n')
        elif self.spectrumType == 'JS':
            print('\tSpectrum Type                        = JONSWAP \n')
        elif self.spectrumType == 'PM':
            print('\tSpectrum Type                        = Pierson-Moskowitz  \n')
        elif self.spectrumType == 'spectrumImport':
            print('\tSpectrum Type                        = Imported Spectrum \n')
            