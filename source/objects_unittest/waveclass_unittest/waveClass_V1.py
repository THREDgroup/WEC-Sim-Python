# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 07:52:07 2020

@author: SungJun Won

This code is written based on WEC-sim. 
wonsungjun0000@gmail.com
"""
from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np
import warnings
import os

class WaveClass:
    def inputProperties():
        """
        Set up necessary properties of the wave condition based on input file
    
        Returns 
        -------
        wType, wT, wH, spectrumType, gamma, phaseSeed, spectrumDataFile,
        etaDataFile, numFreq, waveDir, waveSpread, viz, statisticsDataLoad,
        freqDisc, wavegauge1loc, wavegauge2loc, wavegauge3loc, currentSpeed, 
        currentDirection, currentOption, currentDepth 
    
        """
        
        #wType - String defining the type of waves to be generated.
        #Can be one of: 'noWave', 'noWaveCIC', 'regular', 'regularCIC','irregular',
        #'spectrumImport', and 'etaImport' 
        #(Default = 'NOT DEFINED').
        wType = 'NOT DEFINED'
        
        #wT - [s] Wave period, peak wave period or BEM period.
        #Wave period (regular waves), peak period (irregular waves), or period of 
        #BEM data used for hydrodynamic coefficients ('noWave') 
        #(Default = 'NOT DEFINED').     
        wT = 'NOT DEFINED'
        
        #wH - % [m] Wave height or significant wave height for irregular
        #Wave height (regular waves) or significant wave height (irregular waves) 
        #(Default = 'NOT DEFINED').   
        wH = 'NOT DEFINED'
        
        #spectrumType -  String containing the wave spectrum type
        #Can be one of : 'PM', 'BS', and 'JS' 
        #(Default = 'NOT DEFINED'). 
        spectrumType = 'NOT DEFINED'
        
        #gamma - Only used for 'JS' spectrum type to define gamma 
        #(Default = 3.3)
        gamma = 3.3
        
        #phaseSeed - Only used for irregular waves 
        #if equal to 1,2,3,...,etc, the waves phase is seeded.
        #(Default = 0)
        phaseSeed = 0
        
        #spectrumDataFile - Data file that contains the spectrum data file
        #(Default = 'NOT DEFINED')
        """
        #For MATLAB file use scipi.io to read .mat file instead of .txt file
        #example:
        import scipy.io as sio
        mat_contents = sio.loadmat('spectrumData.mat')
        
        """
        spectrumDataFile = 'NOT DEFINED'
        
        #etaDataFile - Data file that contains the times-series data file
        #(Default = 'NOT DEFINED')
        etaDataFile = 'NOT DEFINED'
        
        #freqRange - Min and max frequency for irregular waves. 
        #2x1 vector, rad/s, (default = frequency range in BEM data)
        #(Default = [])
        freqRange = []
        
        #numFreq - of interpolated wave frequencies 
        #Number of frequencies used, varies depending on method:
        #Traditional = 1000, EqualEnergy = 500 or 'Imported'
        #(Default = 0)
        numFreq = 0
        
        #waveDir - [deg] Incident wave direction(s)
        #Should be defined as a column vector for more than one wave direction        
        #(Default = 0)
        waveDir = 0
        
        #waveSpread - Wave Spread probability associated with wave direction(s)
        #Should be defined as a column vector for more than one wave direction
        #(Default = 1)
        waveSpread = 1
        
        #viz - Dictionary defining visualization options
        #Should be a structure containing the fields 'numPointsX' and 'numPointsY'. 
        #numPointsX is the number of visualization points in x direction, and 
        #numPointsY the number of visualization points in y direction
        viz = {'numPointsX': 50,'numPointsY': 50 }
        
        #statisticsDataLoad - File name from which to load wave statistics data
        #(Default = [])
        statisticsDataLoad = []
        
        #freqDisc - Method of frequency discretization for irregular waves. 
        #Options for this variable are 'EqualEnergy' or 'Traditional'.
        #(Default = 'EqualEnergy').
        freqDisc = 'EqualEnergy'
     
        #wavegauge1loc - [m] Wave gauge 1 [x,y] location
        #(Default = [0,0]).
        wavegauge1loc = [0,0]
        
        #wavegauge2loc - [m] Wave gauge 2 [x,y] location
        #(Default = [0,0]).
        wavegauge2loc = [0,0]
    
        #wavegauge3loc - [m] Wave gauge 3 [x,y] location
        #(Default = [0,0]).
        wavegauge3loc = [0,0]
        
        #currentSpeed - [m/s] Surface current speed that is uniform along the 
        #water column.
        #(Default = 0).
        currentSpeed = 0
      
        #currentDirection - [deg] Surface current direction.
        #(Default = 0).
        currentDirection = 0
        
        #currentOption - [-] Define the sub-surface current model to be used in 
        #WEC-Sim.
        #(Default = 0)
        #0 : Depth-independent model
        #1 : 1/7 power law variation with depth
        #2 : linear variation with depth
        #3 : no current
        currentOption = 3
        
        #currentDepth - [m] Define the depth over which the sub-surface current is 
        #modeled.
        #For options (1) and (2) the currentDepth must be defined. The
        #current is not calculated for any depths greater than the
        #specified currentDepth. (Default = 0).
        currentDepth = 0
        
        return(wType, wT, wH, spectrumType, gamma, phaseSeed, spectrumDataFile,\
               etaDataFile, freqRange, numFreq, waveDir, waveSpread, viz, \
               statisticsDataLoad, freqDisc, wavegauge1loc, wavegauge2loc, \
               wavegauge3loc, currentSpeed, currentDirection, currentOption, \
               currentDepth)
    
    def _internalProperties():
        """
        The following properties are for internal use
    
        Returns
        -------
        typeNum, bemFreq, waterDepth, deepWaterWave, waveAmpTime,waveAmpTime1, 
        waveAmpTime2, waveAmpTime3, wA, wF, phase, dw, wk, wS, wP
    
        """
        
        #typeNum - Number to represent different type of waves
        typeNum = []
        
        #bemFreq - Number of wave frequencies from BEM
        bemFreq = []
        
        #waterDepth - [m] Water depth (from BEM)
        waterDepth = []
        
        #deepWaterWave - Deep water or not, depending on input from WAMIT, 
        #NEMOH and AQWA
        deepWaterWave = []
        
        #waveAmpTime - [m] Wave elevation time history
        waveAmpTime = []
        
        #waveAmpTime1 - [m] Wave elevation time history at a wave gauge 1 location 
        #specified by user
        waveAmpTime1 = []
        
        #waveAmpTime2 - [m] Wave elevation time history at a wave gauge 2 location 
        #specified by user
        waveAmpTime2 = []
        
        #waveAmpTime3 - [m] Wave elevation time history at a wave gauge 3 location 
        #specified by user
        waveAmpTime3 = []
        
        #A - [m] Wave amplitude for regular waves or 2*(wave spectrum vector) for 
        #irregular waves
        wA = [] 
        
        #w - [rad/s] Wave frequency (regular waves) or 
        #wave frequency vector (irregular waves)
        wF = []
        
        #phase - [rad] Wave phase (only used for irregular waves)
        phase = 0
        
        #dw - [rad] Frequency spacing for irregular waves.
        dw = 0 
        
        #k - Wave Number
        wk = []
        
        #S - Wave Spectrum [m^2-s/rad] for 'Traditional'
        wS = []
        
        #Pw - Wave Power Per Unit Wave Crest [W/m]
        wP = []
        
        return(typeNum, bemFreq, waterDepth, deepWaterWave, waveAmpTime,\
               waveAmpTime1, waveAmpTime2, waveAmpTime3, wA, wF, phase, dw, wk,\
               wS, wP)
            
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
        self.wF = []
        
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
    
    def plotEta(self, rampTime):
        """
        Plot wave elevation time-history

        Parameters
        ----------
        rampTime : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        plt.figure(figsize=(10,8))
        plt.plot(self.waveAmpTime[:,0],self.waveAmpTime[:,1])
        plt.title('Wave Surfave Elevation')
        if np.size(rampTime) == 2:
            plt.plot(rampTime,[1.5*min(self.waveAmpTime[:,1]),1.5*max(self.waveAmpTime[:,1])],'Color','k')
            plt.title(['Wave Surfave Elevation, Ramp Time ' + str(rampTime) +' (s)'])
        plt.xlabel('Time (s)')
        plt.ylabel('Eta (m)')
    
    def plotSpectrum(self):
        """
        Plot wave spetrum

        Returns
        -------
        None.

        """
        m0 = np.trapz(self.wS, x = self.wF)
        HsTest = 4*np.sqrt(m0)
        I = np.argmax(np.abs(self.wS))
        wp = self.wF[I]
        TpTest = 2*np.pi/wp
        
        plt.figure(figsize=(10,8))
        plt.plot(self.wF,self.wS,'s-')
        plt.plot(wp,[0,max(self.wS)],'Color','k')
        plt.xlim([0, max(self.wF)])
        plt.title([self.spectrumType, ' Spectrum, T_p= ' + str(TpTest) + ' [s], '  'H_m_0= ' + str(HsTest) + ', [m]'])
        if self.spectrumType == 'JS':
            plt.title([self.spectrumType, ' Spectrum, T_p= ' + str(TpTest) + ' [s], H_m_0= ' + str(HsTest) + ', [m], gamma = ' + str(self.gamma)])
        
        plt.xlabel('Frequency (rad/s)')
        plt.ylabel('Spectrum (m^2-s/rad)')
        
    def waveSetup(self,bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime):
        """
        

        Parameters
        ----------
        bemFreq : TYPE
            DESCRIPTION.
        wDepth : TYPE
            DESCRIPTION.
        rampTime : TYPE
            DESCRIPTION.
        dt : TYPE
            DESCRIPTION.
        maxIt : TYPE
            DESCRIPTION.
        g : TYPE
            DESCRIPTION.
        rho : TYPE
            DESCRIPTION.
        endTime : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.bemFreq = bemFreq
        self.setWaveProps(wDepth)
        if (self.wType == 'noWave') or (self.wType == 'noWaveCIC'):
            if np.size(self.wF) == 0:
                self.wF = 2*np.pi/self.wT
            self.waveNumber(g)
            self.wA = self.wH/2
            self.waveElevNowave(maxIt,dt)
        elif (self.wType == 'regular') or (self.wType == 'regularCIC'):
            if np.size(self.wF) == 0:
                self.wF = 2*np.pi/self.wT
            self.wA = self.wH/2
            self.waveNumber(g)
            self.waveElevReg(rampTime, dt, maxIt)
            self.wavePowerReg(g,rho)
        elif (self.wType == 'irregular') or (self.wType == 'spectrumImport'):
            WFQSt=np.min(bemFreq)
            WFQEd=np.max(bemFreq)
            if np.size(self.freqRange) != 0:
                if self.freqRange[0] > WFQSt and self.freqRange[0] > 0:
                    WFQSt = self.freqRange[0]
                else:
                    warnings.warn("Min frequency range outside BEM data, min frequency set to min BEM frequency",DeprecationWarning)
                if self.freqRange[1] < WFQEd and self.freqRange[1] > WFQSt:
                    WFQEd = self.freqRange[1]
                else:
                    warnings.warn("Max frequency range outside BEM data, max frequency set to max BEM frequency",DeprecationWarning)
            if self.freqDisc == 'Traditional':    
                if np.size(self.numFreq) == 0:
                    self.numFreq = 1000
                self.wF = np.arange(WFQSt,WFQEd+((WFQEd-WFQSt)/(self.numFreq-1)),(WFQEd-WFQSt)/(self.numFreq-1))
                self.dw = np.ones(shape=(self.numFreq,1))*(WFQEd-WFQSt)/(self.numFreq-1)
            elif self.freqDisc == 'EqualEnergy':
                numFreq_interp = 500000
                self.wF = np.arange(WFQSt,WFQEd+((WFQEd-WFQSt)/numFreq_interp),(WFQEd-WFQSt)/numFreq_interp)
                self.dw = np.mean(np.diff(self.wF))
                if np.size(self.numFreq) == 0:
                    self.numFreq = 500
            elif self.freqDisc == 'Imported':
                data = np.conj(np.transpose(np.loadtxt(self.spectrumDataFile)))
                freq_data = data[0]
                self.wF = np.array([i*2*np.pi for i in freq_data 
                                    if i>=min(bemFreq)/2/np.pi and i<=max(bemFreq)/2/np.pi])
                self.numFreq = len(self.wF)
                self.dw = np.zeros(self.numFreq)
                self.dw[0]= np.array(self.wF[1]-self.wF[0])
                self.dw[1:self.numFreq-1]= np.array((self.wF[2:]-self.wF[:-2])/2)
                self.dw[self.numFreq-1]= np.array(self.wF[-1]-self.wF[-2])
            self.setWavePhase()
            self.irregWaveSpectrum(g,rho)
            self.waveNumber(g)
            self.waveElevIrreg(rampTime, dt, maxIt, self.dw)
        elif self.wType == 'etaImport':
            #Import 'etaImport' time-series here and interpolate
            data = np.conj(np.transpose(np.loadtxt(self.etaDataFile))) #Import time-series
            #data = np.loadtxt(self.etaDataFile) #Import time-series 
            t = np.arange(0,endTime+dt,dt)      #simulation time
            self.waveElevUser(rampTime, dt, maxIt, data, t)
            self.waveAmpTime1 = np.zeros((maxIt+1,2))
            self.waveAmpTime2 = np.zeros((maxIt+1,2))
            self.waveAmpTime3 = np.zeros((maxIt+1,2))
            tmp = np.array(np.arange(0,maxIt+1)*dt)
            for i in range(maxIt+1):
                self.waveAmpTime1[i][0] = tmp[i]
                self.waveAmpTime2[i][0] = tmp[i]
                self.waveAmpTime3[i][0] = tmp[i]
            
    def listInfo(self):
        """
        List wave info

        Returns
        -------
        None.

        """
        print('\n Wave Environment: \n')
        if self.wType == 'noWave':
            print('\tWave Type                            = No Wave (Constant Hydrodynamic Coefficients)\n')
            print('\tHydro Data Wave Period, T (sec)    	= {:d}\n'.format(self.wT))
        if self.wType == 'regular':
            print('\tWave Type                            = Regular Waves (Constant Hydrodynamic Coefficients)\n')
            print('\tWave Height, H (m)                   = {:d}\n'.format(self.wH))
            print('\tWave Period, T (sec)                 = {:d}\n'.format(self.wT))
        if self.wType == 'noWaveCIC':
            print('\tWave Type                            = No Wave (Convolution Integral Calculation)\n')
        if self.wType == 'regularCIC':
            print('\tWave Type                            = Regular Waves (Convolution Integral Calculation)\n')
            print('\tWave Height, H (m)                   = {:d}\n'.format(self.wH))
            print('\tWave Period, T (sec)                 = {:d}\n'.format(self.wT))
        if self.wType == 'irregular':
            if self.phaseSeed == 0:
                print('\tWave Type                            = Irregular Waves (Arbitrary Random Phase)\n')
            else:
                print('\tWave Type                            = Irregular Waves (Predefined Random Phase)\n')
            self.printWaveSpectrumType
            print('\tSignificant Wave Height, Hs      (m) = {:d}\n'.format(self.wH))
            if self.spectrumType == 'PM':
                print('\tNOTE: Pierson-Moskowitz does not use Hs to define spectrum\n')
            print('\tPeak Wave Period, Tp           (sec) = {:d}\n'.format(self.wT))
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
    
    def waveNumber(self, g):
        """
        Calculate wave number

        Parameters
        ----------
        g : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.wk= self.wF**2/g
        if self.deepWaterWave == 0:
            for i in range(100):
                self.wk = self.wF**2/g/np.tanh(self.wk*self.waterDepth)
    
    def checkinputs(self):
        """
         Check user inputs
         

        Returns
        -------
        None.

        """
        # 'noWave' period undefined for hydro data
        if (self.wType == 'noWave') and (self.wT == 'NOT DEFINED'):
            warnings.warn('"waves.T" must be defined for the hydrodynamic data period when using the "noWave" wave type',DeprecationWarning)
        # spectrumDataFile defined for 'spectrumImport' case
        if (self.wType == 'spectrumImport') and (self.spectrumDataFile == 'NOT DEFINED'):
            warnings.warn('The spectrumDataFile variable must be defined when using the "spectrumImport" wave type',DeprecationWarning)
        # check waves types
        types = ['noWave', 'noWaveCIC', 'regular', 'regularCIC', 'irregular', 'spectrumImport', 'etaImport']
        if [1 for x in types if x == self.wType] != [1]:
            warnings.warn('Unexpected wave environment type setting, choose from: ' \
                    '"noWave", "noWaveCIC", "regular", "regularCIC", "irregular", "spectrumImport", and "etaImport".',DeprecationWarning)
    
    def write_paraview_vtp(self, t, numPointsX, numPointsY, domainSize, model, simdate, mooring):
        """
        Write vtp files for visualization using Paraview        

        Parameters
        ----------
        t : TYPE
            DESCRIPTION.
        numPointsX : TYPE
            DESCRIPTION.
        numPointsY : TYPE
            DESCRIPTION.
        domainSize : TYPE
            DESCRIPTION.
        model : TYPE
            DESCRIPTION.
        simdate : TYPE
            DESCRIPTION.
        mooring : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        #ground plane
        hDir = os.getcwd()
        os.mkdir('vtk')
        vDir = os.path.join(hDir, 'vtk')
        os.chdir(vDir)
        filename = 'ground.txt'
        fid = open(filename,"w+")
        fid.write(str(domainSize) + "\n")
        fid.write(str(self.waterDepth) + "\n")
        fid.write(str(mooring) + "\n")
        fid.close()
        #wave
        x = np.linspace(-domainSize, domainSize, numPointsX)
        y = np.linspace(-domainSize, domainSize, numPointsY)
        [X,Y] = np.meshgrid(x,y)
        lx = np.size(x,1)
        ly = np.size(y,1)
        numVertex = lx * ly
        numFace = (lx-1) * (ly-1)
        for it in range(0,len(t)):
            #open file
            os.mkdir('waves')
            wDir = os.path.join(hDir, 'waves')
            os.chdir(wDir)
            filenames = 'waves_'+str(it)+'.vtp'
            fids = open(filenames,"w+")
            # calculate wave elevation
            Z = self.waveElevationGrid(t[it], X, Y)
            # write header
            fids.write('<?xml version="1.0"?>\n')
            fids.write('<!-- WEC-Sim Visualization using ParaView -->\n')
            fids.write('<!--   model: ' , model , ' - ran on ' , simdate , ' -->\n')
            fids.write('<!--   wave:  ' , self.wType ,' -->\n')
            fids.write('<!--   time:  ' , str(t(it)) , ' -->\n')
            fids.write('<VTKFile type="PolyData" version="0.1">\n')
            fids.write('  <PolyData>\n')
            # write wave info
            fids.write('    <Piece NumberOfPoints="' + str(numVertex) + '" NumberOfPolys="' + str(numFace) + '">\n')
            # write points
            fids.write('      <Points>\n')
            fids.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
            for jj in range(1,len(y)+1):
                for ii in range(1,len(x)+1):
                    pt = [X[jj][ii], Y[jj][ii], Z[jj][ii]]
                    """
                    if not work try np.array on X,Y,Z
                    """
                    fids.write('          {:5.5f} {:5.5f} {:5.5f}\n'.format(pt[0],pt[1],pt[2]))
            fids.write('        </DataArray>\n')
            fids.write('      </Points>\n')
            # write squares connectivity
            fids.write('      <Polys>\n')
            fids.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
            for jj in range(1,ly):
                for ii in range(1,lx):
                    p1 = (jj-1)*lx + (ii-1)
                    p2 = p1+1
                    p3 = p2 + lx
                    p4 = p1 + lx
                    fids.write('          {} {} {} {}\n'.format(p1,p2,p3,p4))
            fids.write('        </DataArray>\n')
            fids.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
            fids.write('         ')
            for ii in range(1,numFace+1):
                n = ii * 4
                fids.write(' {}'.format(n))
            fids.write('\n')
            fids.write('        </DataArray>\n')
            fids.write('      </Polys>\n')
            # end file
            fids.write('    </Piece>\n')
            fids.write('  </PolyData>\n')
            fids.write('</VTKFile>')
            #close file
            fids.close()
            
    def waveElevationGrid(self, t, X, Y):
        """
        

        Parameters
        ----------
        t : int
            the current time
        X : list
            (m x n) matrix of X coordinates at which to calculate
            the wave elevation
        Y : list
            (m x n) matrix of Y coordinates at which to calculate
            the wave elevation.

        Returns
        -------
        Z - (m x n) matrix of Z coordinates of the wave elevation

        """
        if (self.wType == 'noWave') and (self.wType == 'noWaveCIC') and (self.wType == 'etaImport'):
            Z = np.zeros((np.size(X),np.size(X)))
        elif (self.wType == 'regular') and (self.wType == 'regularCIC'):
            Xt = X*np.cos(self.waveDir*np.pi/180) + Y*np.sin(self.waveDir*np.pi/180)
            Z = self.A*np.cos(-1*self.k*Xt + self.wF*t)
        elif (self.wType == 'irregular') and (self.wType == 'spectrumImport'):
            Z = np.zeros((np.size(X),np.size(X)))
            Xt = X*np.cos(self.waveDir*np.pi/180) + Y*np.sin(self.waveDir*np.pi/180)
            for iw in range(1,len(self.wF)+1):
                Z = Z + np.sqrt(self.wA(iw)*self.dw(iw))*np.cos(-1*self.k(iw)*Xt + self.wF(iw)*t + self.phase(iw))
        return(Z)
    
    def setWavePhase(self):
        """
        Sets the irregular wave's random phase
        Used by waveSetup
        Returns
        -------
        None.

        """
        if self.phaseSeed != 0:
            np.random.seed(self.phaseSeed) #Phase seed = 1,2,3,...,etc
        else:
            np.random.seed(np.random.shuffle(self.phaseSeed))
        if (self.freqDisc == 'EqualEnergy') or (self.freqDisc == 'Traditional'):
            self.phase = 2*np.pi*np.conj(np.transpose(np.random.rand(self.numFreq,np.size(self.waveDir))))
        elif (self.freqDisc == 'Imported'):
            data = np.conj(np.transpose(np.loadtxt(self.spectrumDataFile)))
            if len(data) == 3:
                freq_data = data[0]
                self.phase = np.array([x for x,i in zip(data[2],freq_data) 
                                       if i>=min(self.bemFreq)/2/np.pi and i<=max(self.bemFreq)/2/np.pi])
            else:
                self.phase = 2*np.pi*np.random.rand(1,self.numFreq)[0]
        
    def setWaveProps(self,wDepth):
        """
        Sets global and type-specific properties
        Used by waveSetup
        
        Parameters
        ----------
        wDepth : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if wDepth == 'infinite':
            self.deepWaterWave = 1
            self.waterDepth = 200
            print('\tInfinite water depth specified in BEM, "waves.waterDepth" set to 200m for vizualisation.\n')
        else:
            self.deepWaterWave = 0
            self.waterDepth = wDepth
        if self.wType == 'noWave':
            self.wH = 0
        elif self.wType == 'noWaveCIC':
            self.wH = 0
            if np.size(self.wF) == 0 and self.wT == 'NOT DEFINED':
                self.wF = np.min(self.bemFreq)
                self.wT = 2*np.pi/self.wF
            elif np.size(self.wF) == 0:
                self.wF = 2*np.pi/self.wT
            else:
                self.wT = 2*np.pi/self.wF
        elif self.wType == 'spectrumImport':
            self.wH = 0
            self.wT = 0
            self.freqDisc = 'Imported'
            self.spectrumType = 'spectrumImport'
    
    def waveElevNowave(self,maxIt,dt):
        """
        Set noWave levation time-history

        Parameters
        ----------
        maxIt : TYPE
            DESCRIPTION.
        dt : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.waveAmpTime = np.zeros((maxIt+1,2))
        self.waveAmpTime1 = np.zeros((maxIt+1,2))
        self.waveAmpTime2 = np.zeros((maxIt+1,2))
        self.waveAmpTime3 = np.zeros((maxIt+1,2))
        for i in range(maxIt+1):
            self.waveAmpTime[i][0] = i*dt
            self.waveAmpTime1[i][0] = i*dt
            self.waveAmpTime2[i][0] = i*dt
            self.waveAmpTime3[i][0] = i*dt
        
    def waveElevReg(self, rampTime, dt, maxIt):
        """
        Calculate regular wave elevation time history
        Used by waveSetup

        Parameters
        ----------
        rampTime : TYPE
            DESCRIPTION.
        dt : TYPE
            DESCRIPTION.
        maxIt : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.waveAmpTime = np.zeros((maxIt+1,2))
        self.waveAmpTime1 = np.zeros((maxIt+1,2))
        self.waveAmpTime2 = np.zeros((maxIt+1,2))
        self.waveAmpTime3 = np.zeros((maxIt+1,2))
        maxRampIT = int(np.round(rampTime/dt))
        if rampTime == 0:
            for i in range(maxIt+1):
                t = i*dt
                self.waveAmpTime[i][0]    = t
                self.waveAmpTime[i][1]    = self.wA*np.cos(self.wF*t)
                self.waveAmpTime1[i][0]   = t
                self.waveAmpTime1[i][1]   = self.wA*np.cos(self.wF*t-self.wk*(self.wavegauge1loc[0]*np.cos(self.waveDir[0]*np.pi/180) + self.wavegauge1loc[1]*np.sin(self.waveDir[0]*np.pi/180)))
                self.waveAmpTime2[i][0]   = t
                self.waveAmpTime2[i][1]   = self.wA*np.cos(self.wF*t-self.wk*(self.wavegauge2loc[0]*np.cos(self.waveDir[0]*np.pi/180) + self.wavegauge2loc[1]*np.sin(self.waveDir[0]*np.pi/180)))
                self.waveAmpTime3[i][0]   = t;
                self.waveAmpTime3[i][1]   = self.wA*np.cos(self.wF*t-self.wk*(self.wavegauge3loc[0]*np.cos(self.waveDir[0]*np.pi/180) + self.wavegauge3loc[1]*np.sin(self.waveDir[0]*np.pi/180)))
        else:
            for i in range(maxRampIT):
                t = i*dt
                self.waveAmpTime[i][0]    = t
                self.waveAmpTime[i][1]    = self.wA*np.cos(self.wF*t)*(1+np.cos(np.pi+np.pi*i/maxRampIT))/2
                self.waveAmpTime1[i][0]   = t
                self.waveAmpTime1[i][1]   = self.wA*np.cos(self.wF*t-self.wk*(self.wavegauge1loc[0]*np.cos(self.waveDir[0]*np.pi/180) + self.wavegauge1loc[1]*np.sin(self.waveDir[0]*np.pi/180)))*(1+np.cos(np.pi+np.pi*i/maxRampIT))/2
                self.waveAmpTime2[i][0]   = t
                self.waveAmpTime2[i][1]   = self.wA*np.cos(self.wF*t-self.wk*(self.wavegauge2loc[0]*np.cos(self.waveDir[0]*np.pi/180) + self.wavegauge2loc[1]*np.sin(self.waveDir[0]*np.pi/180)))*(1+np.cos(np.pi+np.pi*i/maxRampIT))/2
                self.waveAmpTime3[i][0]   = t
                self.waveAmpTime3[i][1]   = self.wA*np.cos(self.wF*t-self.wk*(self.wavegauge3loc[0]*np.cos(self.waveDir[0]*np.pi/180) + self.wavegauge3loc[1]*np.sin(self.waveDir[0]*np.pi/180)))*(1+np.cos(np.pi+np.pi*i/maxRampIT))/2
            for i in range(maxRampIT,maxIt+1):
                t = i*dt
                self.waveAmpTime[i][0]    = t
                self.waveAmpTime[i][1]    = self.wA*np.cos(self.wF*t)
                self.waveAmpTime1[i][0]   = t
                self.waveAmpTime1[i][1]   = self.wA*np.cos(self.wF*t-self.wk*(self.wavegauge1loc[0]*np.cos(self.waveDir[0]*np.pi/180) + self.wavegauge1loc[1]*np.sin(self.waveDir[0]*np.pi/180)))
                self.waveAmpTime2[i][0]   = t
                self.waveAmpTime2[i][1]   = self.wA*np.cos(self.wF*t-self.wk*(self.wavegauge2loc[0]*np.cos(self.waveDir[0]*np.pi/180) + self.wavegauge2loc[1]*np.sin(self.waveDir[0]*np.pi/180)))
                self.waveAmpTime3[i][0]   = t
                self.waveAmpTime3[i][1]   = self.wA*np.cos(self.wF*t-self.wk*(self.wavegauge3loc[0]*np.cos(self.waveDir[0]*np.pi/180) + self.wavegauge3loc[1]*np.sin(self.waveDir[0]*np.pi/180)))
    
    def wavePowerReg(self,g,rho):
        """
        Calculate wave power per unit wave crest for regular waves

        Parameters
        ----------
        g : TYPE
            DESCRIPTION.
        rho : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if self.deepWaterWave == 1:
            # Deepwater Approximation
            self.wP = 1/(8*np.pi)*rho*g**(2)*(self.wA)**(2)*self.wT
        else:
            # Full Wave Power Equation
            self.wP = rho*g*(self.wA)**(2)/4*np.sqrt(g/self.wk*np.tanh(self.wk*self.waterDepth))*(1+2*self.wk*self.waterDepth/np.sinh(self.wk*self.waterDepth))
        
    def irregWaveSpectrum(self,g,rho):
        """
        Calculate wave spectrum vector (self.A)
        Used by wavesIrreg (wavesIrreg used by waveSetup)

        Parameters
        ----------
        g : TYPE
            DESCRIPTION.
        rho : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        freq = self.wF/(2*np.pi)
        Tp = self.wT
        Hs = self.wH
        if self.spectrumType == 'PM':
            # Pierson-Moskowitz Spectrum from Tucker and Pitt (2001)
            B_PM = (5/4)*(1/Tp)**(4)
            A_PM = 0.0081*g**2*(2*np.pi)**(-4)
            S_f  = (A_PM*freq**(-5)*np.exp(-B_PM*freq**(-4)))              # Wave Spectrum [m^2-s] for 'EqualEnergy'
            self.wS = S_f/(2*np.pi)                                      # Wave Spectrum [m^2-s/rad] for 'Traditional'
            S_f = self.wS*2*np.pi
        elif self.spectrumType == 'BS':
            # Bretschneider Sprectrum from Tucker and Pitt (2001)
            B_BS = (1.057/Tp)**4
            A_BS = B_BS*(Hs/2)**2
            S_f = (A_BS*freq**(-5)*np.exp(-B_BS*freq**(-4)))               # Wave Spectrum [m^2-s]
            self.wS = S_f/(2*np.pi)                                      # Wave Spectrum [m^2-s/rad]
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
            self.wS = S_f/(2*np.pi)                                       # Wave Spectrum [m^2-s/rad]
        elif self.spectrumType == 'spectrumImport':
            # Imported Wave Spectrum
            data = np.conj(np.transpose(np.loadtxt(self.spectrumDataFile)))
            freq_data = data[0]
            S_data = data[1]
            S_f = np.array([x for x,i in zip(S_data,freq_data) 
                            if i>=min(self.bemFreq)/2/np.pi and i<=max(self.bemFreq)/2/np.pi])       # Wave Spectrum [m^2-s] for 'EqualEnergy'
            self.wS = S_f/(2*np.pi)                                      # Wave Spectrum [m^2-s/rad] for 'Traditional'
            print('\t"spectrumImport" uses the number of imported wave frequencies (not "Traditional" or "EqualEnergy")\n')
            # Power per Unit Wave Crest
        self.waveNumber(g)                                          # Calculate Wave Number for Larger Number of Frequencies Before Down Sampling in Equal Energy Method
        if self.deepWaterWave == 1:
            # Deepwater Approximation
            self.Pw = np.sum(1/2*rho*g**(2)*S_f*self.dw/self.wF)
        else:
            # Full Wave Power Equation
            self.Pw = np.sum((1/2)*rho*g*S_f*self.dw*np.sqrt(9.81/self.wk*np.tanh(self.wk*self.waterDepth))*(1 + 2*self.wk*self.waterDepth/np.sinh(2*self.wk*self.waterDepth)))
        if self.freqDisc == 'EqualEnergy':
            m0 = np.trapz(np.abs(S_f),freq)
            numBins = self.numFreq+1
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
            self.wF = 2*np.pi*freq[wn[1:len(wn)-1]]
            self.dw = np.append(self.wF[0]-2*np.pi*freq[wn[0]-1],np.diff(self.wF)) 
            self.wS = self.wS[wn[1:len(wn)-1]]                         # Wave Spectrum [m^2-s/rad] 
        self.wA = 2 * self.wS                                             # Wave Amplitude [m]
           
    def waveElevIrreg(self,rampTime,dt,maxIt,df):
        """
        Calculate irregular wave elevetaion time history
        Used by waveSetup
        
        Parameters
        ----------
        rampTime : TYPE
            DESCRIPTION.
        dt : TYPE
            DESCRIPTION.
        maxIt : TYPE
            DESCRIPTION.
        df : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.waveAmpTime = np.zeros((maxIt+1,2))
        self.waveAmpTime1 = np.zeros((maxIt+1,2))
        self.waveAmpTime2 = np.zeros((maxIt+1,2))
        self.waveAmpTime3 = np.zeros((maxIt+1,2))
        maxRampIT=int(np.round(rampTime/dt))
        if rampTime == 0:
            for i in range(maxIt+1):
                for idir in range(np.size(self.waveDir)):
                    t       = (i)*dt
                    tmp     = np.sqrt(self.wA*df*self.wavepSread[idir])
                    tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t + self.phase[:idir])))
                    tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t - self.wk*(self.wavegauge1loc[0]*np.cos(self.waveDir[idir]*np.pi/180) + self.wavegauge1loc[1]*np.sin(self.waveDir[idir]*np.pi/180)) + self.phase[:idir])))
                    tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t - self.wk*(self.wavegauge2loc[0]*np.cos(self.waveDir[idir]*np.pi/180) + self.wavegauge2loc[1]*np.sin(self.waveDir[idir]*np.pi/180)) + self.phase[:idir])))
                    tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t - self.wk*(self.wavegauge3loc[0]*np.cos(self.waveDir[idir]*np.pi/180) + self.wavegauge3loc[1]*np.sin(self.waveDir[idir]*np.pi/180)) + self.phase[:idir])))
                    self.waveAmpTime[i][0]    = t
                    self.waveAmpTime[i][1]    = self.waveAmpTime[i][1] + sum(tmp1)
                    self.waveAmpTime1[i][0]   = t
                    self.waveAmpTime1[i][1]   = self.waveAmpTime1[i][1] + sum(tmp11)
                    self.waveAmpTime2[i][0]   = t
                    self.waveAmpTime2[i][1]   = self.waveAmpTime2[i][1] + sum(tmp12)
                    self.waveAmpTime3[i][0]   = t
                    self.waveAmpTime3[i][1]   = self.waveAmpTime3[i][1] + sum(tmp13)
        else:
            for i in range(maxRampIT):
                for idir in range(np.size(self.waveDir)):
                    t       = (i)*dt
                    tmp     = np.sqrt(self.wA*np.array(df)*self.waveSpread[idir])
                    tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t + self.phase[idir])))
                    tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t - self.wk*(self.wavegauge1loc[0]*np.cos(self.waveDir[idir]*np.pi/180) + self.wavegauge1loc[1]*np.sin(self.waveDir[idir]*np.pi/180)) + self.phase[idir])))
                    tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t - self.wk*(self.wavegauge2loc[0]*np.cos(self.waveDir[idir]*np.pi/180) + self.wavegauge2loc[1]*np.sin(self.waveDir[idir]*np.pi/180)) + self.phase[idir])))
                    tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t - self.wk*(self.wavegauge3loc[0]*np.cos(self.waveDir[idir]*np.pi/180) + self.wavegauge3loc[1]*np.sin(self.waveDir[idir]*np.pi/180)) + self.phase[idir])))
                    self.waveAmpTime[i][0]    = t
                    self.waveAmpTime[i][1]    = self.waveAmpTime[i][1] + sum(tmp1)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
                    self.waveAmpTime1[i][0]   = t
                    self.waveAmpTime1[i][1]   = self.waveAmpTime1[i][1] + sum(tmp11)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
                    self.waveAmpTime2[i][0]   = t
                    self.waveAmpTime2[i][1]   = self.waveAmpTime2[i][1] + sum(tmp12)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
                    self.waveAmpTime3[i][0]   = t
                    self.waveAmpTime3[i][1]   = self.waveAmpTime3[i][1] + sum(tmp13)*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
            for i in range(maxRampIT, maxIt+1):
                for idir in range (np.size(self.waveDir)):
                    t       = (i)*dt
                    tmp     = np.sqrt(self.wA*df*self.waveSpread[idir])
                    tmp1    = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t + self.phase[idir])))
                    tmp11   = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t - self.wk*(self.wavegauge1loc[0]*np.cos(self.waveDir[idir]*np.pi/180) + self.wavegauge1loc[1]*np.sin(self.waveDir[idir]*np.pi/180)) + self.phase[idir])))
                    tmp12   = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t - self.wk*(self.wavegauge2loc[0]*np.cos(self.waveDir[idir]*np.pi/180) + self.wavegauge2loc[1]*np.sin(self.waveDir[idir]*np.pi/180)) + self.phase[idir])))
                    tmp13   = tmp*np.real(np.exp(((-1)**0.5)*(self.wF*t - self.wk*(self.wavegauge3loc[0]*np.cos(self.waveDir[idir]*np.pi/180) + self.wavegauge3loc[1]*np.sin(self.waveDir[idir]*np.pi/180)) + self.phase[idir])))
                    self.waveAmpTime[i][0]    = t
                    self.waveAmpTime[i][1]    = self.waveAmpTime[i][1] + sum(tmp1)
                    self.waveAmpTime1[i][0]   = t
                    self.waveAmpTime1[i][1]   = self.waveAmpTime1[i][1] + sum(tmp11)
                    self.waveAmpTime2[i][0]   = t
                    self.waveAmpTime2[i][1]   = self.waveAmpTime2[i][1] + sum(tmp12)
                    self.waveAmpTime3[i][0]   = t
                    self.waveAmpTime3[i][1]   = self.waveAmpTime3[i][1] + sum(tmp13)
    
    def waveElevUser(self,rampTime,dt,maxIt,data,t):
        """
        Calculate imported wave elevation time history
        Used by waveSetup

        Parameters
        ----------
        rampTime : TYPE
            DESCRIPTION.
        dt : TYPE
            DESCRIPTION.
        maxIt : TYPE
            DESCRIPTION.
        data : TYPE
            DESCRIPTION.
        t : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.waveAmpTime = np.zeros((maxIt+1,2))
        maxRampIT = int(np.round(rampTime/dt))
        data_t = data[0]                    # Data Time [s]
        data_x = data[1]                    # Wave Surface Elevation [m] 
        for i in range(len(t)):
            self.waveAmpTime[i][0] = t[i]
            self.waveAmpTime[i][1] = np.interp(t,data_t,data_x)[i]
        if rampTime != 0:
            for i in range(maxRampIT+1):
                self.waveAmpTime[i][0] = t[i]
                self.waveAmpTime[i][1] = self.waveAmpTime[i][1]*(1+np.cos(np.pi+np.pi*(i)/maxRampIT))/2
                    
    
    def printWaveSpectrumType(self):
        """
        Lists the wave spectrum type
        Used by listInfo

        Returns
        -------
        None.

        """
        if self.spectrumType == 'BS':
            print('\tSpectrum Type                        = Bretschneider \n')
        elif self.spectrumType == 'JS':
            print('\tSpectrum Type                        = JONSWAP \n')
        elif self.spectrumType == 'PM':
            print('\tSpectrum Type                        = Pierson-Moskowitz  \n')
        elif self.spectrumType == 'spectrumImport':
            print('\tSpectrum Type                        = Imported Spectrum \n')
            