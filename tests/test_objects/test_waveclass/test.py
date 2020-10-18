#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 08:43:15 2020

@author: logical
"""

import unittest
import numpy as np

from waveclass import WaveClass

class TestWave(unittest.TestCase):
    
    def setUpClass():
        print("setupClass")
        
    def tearDownClass():
        print("teardownClass")
    
    def setUp(self):
        """
        set up test cases for each wave type.
        each test cases have dependent test data files in the testData directory
        """
        print("setUp")
        
        self.noWave_1 = WaveClass('noWave')
        
        self.noWaveCIC_1 = WaveClass('noWaveCIC')
        self.noWaveCIC_2 = WaveClass('noWaveCIC')
        self.noWaveCIC_3 = WaveClass('noWaveCIC')
        
        self.regular_1 = WaveClass('regular')
        
        self.regularCIC_1 = WaveClass('regularCIC')
        
        self.irregular_1 = WaveClass('irregular')
        self.irregular_2 = WaveClass('irregular')
        self.irregular_3 = WaveClass('irregular')
        self.irregular_4 = WaveClass('irregular')
        self.irregular_5 = WaveClass('irregular')
        
        self.spectrumImport_1 = WaveClass('spectrumImport')
        self.spectrumImport_2 = WaveClass('spectrumImport')
        
        self.etaImport_1 = WaveClass('etaImport')
        self.etaImport_2 = WaveClass('etaImport')
        
    def tearDown(self):
        print("tearDown\n")
    
    def test_setWaveProps(self):
        """
        Test SetWaveProps
        """
        print("SetWaveProps")
        
        wDepth = 'infinite'
        self.noWave_1.setWaveProps(wDepth)
        self.assertEqual(self.noWave_1.H, 0)
        
        # wDepth = 'infinite'
        self.noWaveCIC_1.T = 'NOT DEFINED'
        self.noWaveCIC_1.bemFreq = [5.19999512307279,0.0199999977946844]
        self.noWaveCIC_1.setWaveProps(wDepth)
        self.assertEqual(self.noWaveCIC_1.H, 0)
        self.assertEqual(self.noWaveCIC_1.w, 0.0199999977946844)
        self.assertIsNone(np.testing.assert_allclose(self.noWaveCIC_1.T, [314.159300000000]))
        
        # wDepth = 'infinite'
        self.noWaveCIC_2.T = 8
        self.noWaveCIC_2.bemFreq = [5.19999512307279,0.0199999977946844]
        self.noWaveCIC_2.setWaveProps(wDepth)
        self.assertEqual(self.noWaveCIC_2.H, 0)
        self.assertIsNone(np.testing.assert_allclose(self.noWaveCIC_2.w, [0.785398163397448]))
        
        # wDepth = 'infinite'
        self.noWaveCIC_3.w = 0.785398163397448
        self.noWaveCIC_3.bemFreq = [5.19999512307279,0.0199999977946844]
        self.noWaveCIC_3.setWaveProps(wDepth)
        self.assertEqual(self.noWaveCIC_3.H, 0)
        self.assertIsNone(np.testing.assert_allclose(self.noWaveCIC_3.T, [8]))
        
        # wDepth = 'infinite'
        self.regular_1.setWaveProps(wDepth)
        self.assertEqual(self.regular_1.deepWaterWave, 1)
        self.assertEqual(self.regular_1.waterDepth, 200)
        
        # wDepth = 'infinite'
        self.regularCIC_1.setWaveProps(wDepth)
        self.assertEqual(self.regularCIC_1.deepWaterWave, 1)
        self.assertEqual(self.regularCIC_1.waterDepth, 200)
        
        # wDepth = 'infinite'
        self.irregular_1.setWaveProps(wDepth)
        self.assertEqual(self.irregular_1.deepWaterWave, 1)
        self.assertEqual(self.irregular_1.waterDepth, 200)
        
        # wDepth = 'infinite'
        self.spectrumImport_1.setWaveProps(wDepth)
        self.assertEqual(self.spectrumImport_1.deepWaterWave, 1)
        self.assertEqual(self.spectrumImport_1.waterDepth, 200)
        self.assertEqual(self.spectrumImport_1.H, 0)
        self.assertEqual(self.spectrumImport_1.T, 0)
        self.assertEqual(self.spectrumImport_1.freqDisc, 'Imported')
        self.assertEqual(self.spectrumImport_1.spectrumType, 'spectrumImport')
        
    def test_setWavePhase(self):
        """
        Test SetWavePhase
        for single wave direction it creates single array of phase values
        for multiple wave direction it creates multiple arrays of phase values
        """
        
        print("SetWavePhase")
        self.irregular_1.waveDir = [0]
        self.irregular_1.phaseSeed = 1
        self.irregular_1.freqDisc = 'EqualEnergy'
        self.irregular_1.numFreq = 4
        self.irregular_1.setWavePhase()
        result = [[2.62022653271779,4.52593227359735,0.000718638171852741,1.89961157824218]]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.phase, result))
        
        self.irregular_2.waveDir = [0]
        self.irregular_2.phaseSeed = 1
        self.irregular_2.freqDisc = 'Traditional'
        self.irregular_2.numFreq = 4
        self.irregular_2.setWavePhase()
        result = [[2.62022653271779,4.52593227359735,0.000718638171852741,1.89961157824218]]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_2.phase, result))
        
        self.irregular_4.waveDir = [0,90]
        self.irregular_4.phaseSeed = 1
        self.irregular_4.freqDisc = 'EqualEnergy'
        self.irregular_4.numFreq = 4
        self.irregular_4.setWavePhase()
        result = [[2.62022653271779,0.000718638171852741,0.922094456924136,1.17030742344035],[4.52593227359735,1.89961157824218,0.580180501936920,2.17122208289517]]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_4.phase, result))
        
        self.spectrumImport_1.waveDir = [0]
        self.spectrumImport_1.phaseSeed = 1
        self.spectrumImport_1.freqDisc = 'Imported'
        self.spectrumImport_1.bemFreq = [5.19999512307279,0.0199999977946844]
        self.spectrumImport_1.spectrumDataFile = "./testData/spectrumImport_1_test/spectrumData.txt"
        self.spectrumImport_1.setWavePhase()
        result = np.array([np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/phase.txt")))])
        self.assertIsNone(np.testing.assert_allclose(self.spectrumImport_1.phase, result))
        
        self.spectrumImport_2.waveDir = [0]
        self.spectrumImport_2.phaseSeed = 1
        self.spectrumImport_2.freqDisc = 'Imported'
        self.spectrumImport_2.bemFreq = [5.19999512307279,0.0199999977946844]
        self.spectrumImport_2.spectrumDataFile = "./testData/spectrumImport_1_test/spectrumData.mat"
        self.spectrumImport_2.setWavePhase()
        result = np.array([np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/phase.txt")))])
        self.assertIsNone(np.testing.assert_allclose(self.spectrumImport_2.phase, result))
        
    def test_waveNumber(self):
        """
        Test waveNumber
        """
        print("waveNumber")
        g = 9.81
        self.noWave_1.w = np.array([0.785398163397448])
        self.noWave_1.deepWaterWave = 1
        self.noWave_1.waveNumber(g)
        result =  [0.0628797426165224]
        self.assertIsNone(np.testing.assert_allclose(self.noWave_1.k, result))
        
        # g = 9.81
        self.regular_1.w = np.array([0.785398163397448])
        self.regular_1.deepWaterWave = 1
        self.regular_1.waveNumber(g)
        result =  [0.0628797426165224]
        self.assertIsNone(np.testing.assert_allclose(self.regular_1.k, result))
        
        # g = 9.81
        self.irregular_1.w = np.array([0.525743641856087,0.541511547017433,0.551933697209493,0.559973049643924,0.566624163384781])
        self.irregular_1.deepWaterWave = 1
        self.irregular_1.waveNumber(g)
        result =  [0.0281759813406831,0.0298914123907455,0.0310530893083935,0.0319643033973004,0.0327281286984203]  
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.k, result))
        
        # g = 9.81
        self.irregular_2.w = np.array([0.525743641856087,0.541511547017433,0.551933697209493,0.559973049643924,0.566624163384781])
        self.irregular_2.deepWaterWave = 0
        self.irregular_2.waterDepth = 150
        self.irregular_2.waveNumber(g)
        result =  [0.0281879608383217,0.0298990180718354,0.0310586687570460,0.0319686743818314,0.0327316884302928]  
        self.assertIsNone(np.testing.assert_allclose(self.irregular_2.k, result))
        
        # g = 9.81
        self.spectrumImport_1.w = np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/w.txt")))
        self.spectrumImport_1.deepWaterWave = 1
        self.spectrumImport_1.waveNumber(g)
        result =  np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/k.txt"))) 
        self.assertIsNone(np.testing.assert_allclose(self.spectrumImport_1.k, result))
        
        
    def test_irregWaveSpectrum(self):
        """
        Test irregWaveSpectrum
        """
    
        print("irregWaveSpectrum")
        g = 9.81
        rho = 1000
        self.irregular_1.w = np.array([0.525743641856087,0.541511547017433,0.551933697209493,0.559973049643924,0.566624163384781,0.572373957973840,0.577471073177114,0.582081268838611,0.586297784870588,0.590203501195047,0.593850217763243,0.597279374536177,0.600522051484601,0.603609328579267,0.606551565810425])
        self.irregular_1.T = 8
        self.irregular_1.H = 2.5
        self.irregular_1.spectrumType = 'BS'
        self.irregular_1.freqDisc = 'EqualEnergy'
        self.irregular_1.deepWaterWave = 1
        self.irregular_1.numFreq =4
        self.irregular_1.dw = np.array([0.505743644061402,0.0157679051613465,0.0104221501920596,0.00803935243443155,0.00665111374085714,0.00574979458905867,0.00509711520327372,0.00461019566149745,0.00421651603197637,0.00390571632445969,0.00364671656819582,0.00342915677293409,0.00324267694842406,0.00308727709466572,0.00294223723115805])
        self.irregular_1.irregWaveSpectrum(g,rho)
        result1 = [0.566624163384781,0.582081268838611,0.590203501195047,0.600522051484601]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.w, result1))
        result2 = [0.0408805215286940,0.0154571054538299,0.00812223235643605,0.0103185502895540]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.dw, result2))
        result3 = [0.126744401973820,0.177305360376672,0.206784788616044,0.246458761526240]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.S, result3))
        result4 = [0.253488803947639,0.354610720753343,0.413569577232088,0.492917523052480]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.A, result4))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.Pw, [16836.8900561635]))
        
        # g = 9.81
        # rho = 1000
        self.irregular_2.w = np.array([0.525743641856087,0.541511547017433,0.551933697209493,0.559973049643924])
        self.irregular_2.T = 8
        self.irregular_2.H = 2.5
        self.irregular_2.spectrumType = 'PM'
        self.irregular_2.freqDisc = 'Traditonal'
        self.irregular_2.deepWaterWave = 1
        self.irregular_2.dw = np.array([0.505743644061402,0.0157679051613465,0.0104221501920596,0.00803935243443155])
        self.irregular_2.irregWaveSpectrum(g,rho)
        result1 = [0.0383934968538238,0.0662991701526572,0.0904660461442089,0.112249475236517]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_2.S, result1))
        result2 = [0.0767869937076475,0.132598340305314,0.180932092288418,0.224498950473034]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_2.A, result2))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_2.Pw, [12753.4644033446]))
        
        # g = 9.81
        # rho = 1000
        self.irregular_3.w = np.array([0.525743641856087,0.541511547017433,0.551933697209493,0.559973049643924])
        self.irregular_3.T = 8
        self.irregular_3.H = 2.5
        self.irregular_3.spectrumType = 'JS'
        self.irregular_3.freqDisc = 'Traditonal'
        self.irregular_3.deepWaterWave = 0
        self.irregular_3.waterDepth = 100
        self.irregular_3.dw = np.array([0.505743644061402,0.0157679051613465,0.0104221501920596,0.00803935243443155])
        self.irregular_3.irregWaveSpectrum(g,rho)
        result1 = [6.10293411905465,10.5392404574227,14.3820944077894,17.8473515284263]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_3.S, result1))
        result2 = [12.2058682381093,21.0784809148455,28.7641888155788,35.6947030568527]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_3.A, result2))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_3.Pw, [2089159.87712874]))
        
        # g = 9.81
        # rho = 1000
        self.spectrumImport_1.freqDisc = 'Imported'
        self.spectrumImport_1.w = np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/w.txt")))
        self.spectrumImport_1.spectrumDataFile = "./testData/spectrumImport_1_test/spectrumData.txt"
        self.spectrumImport_1.H = 0
        self.spectrumImport_1.T = 0
        self.spectrumImport_1.spectrumType = 'spectrumImport'
        self.spectrumImport_1.bemFreq = [5.19999512307279,0.0199999977946844]
        self.spectrumImport_1.deepWaterWave = 1
        self.spectrumImport_1.dw = np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/dw.txt")))
        self.spectrumImport_1.phase = np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/phase.txt")))
        self.spectrumImport_1.irregWaveSpectrum(g,rho)
        result1 = np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/S.txt")))
        self.assertIsNone(np.testing.assert_allclose(self.spectrumImport_1.S, result1))
        result2 = np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/A.txt")))
        self.assertIsNone(np.testing.assert_allclose(self.spectrumImport_1.A, result2))
        self.assertIsNone(np.testing.assert_allclose(self.spectrumImport_1.Pw, [65833.3074598130]))
        
    def test_waveElevNowave(self):
        """
        #Test waveElevNowave
        """
        print("waveElevNowave")
        dt = 0.1
        maxIt = 4
        self.noWave_1.waveElevNowave(maxIt,dt)
        result = [[0,0.1,0.2,0.3,0.4],[0,0,0,0,0]]
        self.assertIsNone(np.testing.assert_allclose(self.noWave_1.waveAmpTime, result))    
        
        
    def test_waveElevReg(self):
        """
        #Test waveElevReg
        """
        print("waveElevReg")
        rampTime = 0.2
        dt = 0.1
        maxIt = 4
        self.regular_1.waveDir = [0]
        self.regular_1.A = np.array([1.25])
        self.regular_1.waveSpread = [1]
        self.regular_1.w = np.array([0.7853981633974481])
        self.regular_1.k = np.array([0.0628797426165224])
        self.regular_1.wavegauge1loc = [0,0]
        self.regular_1.wavegauge2loc = [0,0]
        self.regular_1.wavegauge3loc = [0,0]
        self.regular_1.waveElevReg(rampTime, dt, maxIt)
        result = [[0,0.1,0.2,0.3,0.4],[0,0.623073333583205,1.23461042574392,1.21546240049710,1.18882064536894]]
        self.assertIsNone(np.testing.assert_allclose(self.regular_1.waveAmpTime, result))
        
    def test_wavePowerReg(self):
        """
        Test wavePowerReg
        """
        print("wavePowerReg")
        g  = 9.81
        rho = 1000
        self.regular_1.deepWaterWave = 1
        self.regular_1.A = np.array([1.25])
        self.regular_1.T = np.array([8])
        self.regular_1.wavePowerReg(g, rho)
        result = [47863.9094340186]
        self.assertIsNone(np.testing.assert_allclose(self.regular_1.Pw, result))
        
        # g  = 9.81
        # rho = 1000
        self.regularCIC_1.deepWaterWave = 0
        self.regularCIC_1.A = np.array([1.25])
        self.regularCIC_1.T = np.array([8])
        self.regularCIC_1.k = np.array([0.0628797426165224])
        self.regularCIC_1.waterDepth = np.array([200])
        self.regularCIC_1.wavePowerReg(g, rho)
        result = [47872.2259960703]
        self.assertIsNone(np.testing.assert_allclose(self.regularCIC_1.Pw, result))
        
    def test_waveElevIrreg(self):
        """
        Test waveElevIrreg
        """
        
        print("waveElevIrreg")
        rampTime = 0.2
        dt = 0.1
        maxIt = 4
        df = np.array([0.0408805215286940,0.0154571054538299,0.00812223235643605,0.0103185502895540])
        self.irregular_1.waveDir = [0]
        self.irregular_1.A = np.array([0.253488803947639,0.354610720753343,0.413569577232088,0.492917523052480])
        self.irregular_1.waveSpread = [1]
        self.irregular_1.w = np.array([0.566624163384781,0.582081268838611,0.590203501195047,0.600522051484601])
        self.irregular_1.phase = np.array([[2.62022653271779,4.52593227359735,0.000718638171852741,1.89961157824218]])
        self.irregular_1.k = np.array([0.0327281286984203,0.0345380839482943,0.0355086822449431,0.0367611349968679])
        self.irregular_1.wavegauge1loc = [0,0]
        self.irregular_1.wavegauge2loc = [0,0]
        self.irregular_1.wavegauge3loc = [0,0]
        self.irregular_1.waveElevIrreg(rampTime, dt, maxIt, df)
        result = [[0,0.1,0.2,0.3,0.4],[0,-0.0348281635986628,-0.0720229579155233,-0.0741609426395505,-0.0760625424978483]]
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.waveAmpTime, result))
        
        rampTime = 100
        dt = 0.1
        maxIt = 2000
        df = np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/dw.txt")))
        self.spectrumImport_1.waveDir = [0]
        self.spectrumImport_1.A = np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/A.txt")))
        self.spectrumImport_1.waveSpread = [1]
        self.spectrumImport_1.w = np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/w.txt")))
        self.spectrumImport_1.phase = np.array([np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/phase.txt")))])
        self.spectrumImport_1.k = np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/k.txt")))
        self.spectrumImport_1.wavegauge1loc = [0,0]
        self.spectrumImport_1.wavegauge2loc = [0,0]
        self.spectrumImport_1.wavegauge3loc = [0,0]
        self.spectrumImport_1.waveElevIrreg(rampTime, dt, maxIt, df)
        result = np.conj(np.transpose(np.loadtxt("./testData/spectrumImport_1_test/waveAmpTime.txt")))
        self.assertIsNone(np.testing.assert_allclose(self.spectrumImport_1.waveAmpTime, result))
        
        
        
    def test_waveElevUser(self):
        print("waveElevUser")
        
        rampTime = 100
        dt = 0.1
        maxIt = 4000
        endTime = 400
        t = np.arange(0,endTime+dt,dt) 
        data = np.conj(np.transpose(np.loadtxt("./testData/etaImport_1_test/etaData.txt")))
        self.etaImport_1.waveElevUser(rampTime,dt,maxIt,data,t)
        result = np.conj(np.transpose(np.loadtxt("./testData/etaImport_1_test/waveAmpTime.txt")))
        self.assertIsNone(np.testing.assert_allclose(self.etaImport_1.waveAmpTime, result))
        
        
        
    def test_waveSetup(self):
        """
        commented repeated code but left it as commented so that we can see the values for each tests.
        since only the max and the min value of beam frequency is used in the actual code,
        I changed beam frequency to an array with two value for faster computation
        However, actual beam frequency can be used to test for confirmation.

        """
        print("waveSetup")
        
        # wType = 'noWave'
        # bemFreq = np.array(np.loadtxt('regular_1_test/beamFreq.txt'))
        bemFreq = [5.19999512307279,0.0199999977946844]
        wDepth = 'infinite'
        rampTime = 100
        dt = 0.1
        maxIt = 2000
        g = 9.81
        rho = 1000
        endTime = 200
        self.noWave_1.T = np.array([8])
        self.noWave_1.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)
        result1=  np.conj(np.transpose(np.array(np.loadtxt('./testData/noWave_1_test/waveAmpTime.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.noWave_1.waveAmpTime, result1))
        result2 =  [0.785398163397448]
        self.assertIsNone(np.testing.assert_allclose(self.noWave_1.w, result2))
        
        
        # wType = 'regular'
        # bemFreq = np.array(np.loadtxt('regular_1_test/beamFreq.txt'))
        # bemFreq = [5.19999512307279,0.0199999977946844]
        # wDepth = 'infinite'
        # rampTime = 100
        # dt = 0.1
        # maxIt = 2000
        # g = 9.81
        # rho = 1000
        # endTime = 200
        self.regular_1.T = np.array([8])
        self.regular_1.H = np.array([2.5])
        self.regular_1.wavegauge1loc = [5,5]
        self.regular_1.wavegauge2loc = [10,0]
        self.regular_1.wavegauge3loc = [0,-10]
        self.regular_1.waveDir = [0]
        self.regular_1.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)
        result1 =  np.conj(np.transpose(np.array(np.loadtxt('./testData/regular_1_test/waveAmpTime.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.regular_1.waveAmpTime, result1))
        result2 =  np.conj(np.transpose(np.array(np.loadtxt('./testData/regular_1_test/waveAmpTime1.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.regular_1.waveAmpTime1, result2))
        result3 =  np.conj(np.transpose(np.array(np.loadtxt('./testData/regular_1_test/waveAmpTime2.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.regular_1.waveAmpTime2, result3))
        result4 =  np.conj(np.transpose(np.array(np.loadtxt('./testData/regular_1_test/waveAmpTime3.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.regular_1.waveAmpTime3, result4))
        result5 =  [0.785398163397448]
        self.assertIsNone(np.testing.assert_allclose(self.regular_1.w, result5))
        
        
        # wType = 'regularCIC'
        # bemFreq = np.array(np.loadtxt('regularCIC_1_test/beamFreq.txt'))
        # bemFreq = [5.19999512307279,0.0199999977946844]
        # wDepth = 'infinite'
        # rampTime = 100
        # dt = 0.1
        # maxIt = 2000
        # g = 9.81
        # rho = 1000
        # endTime = 200
        self.regularCIC_1.T = np.array([8])
        self.regularCIC_1.H = np.array([2.5])
        self.regularCIC_1.wavegauge1loc = [0,0]
        self.regularCIC_1.wavegauge2loc = [0,0]
        self.regularCIC_1.wavegauge3loc = [0,0]
        self.regularCIC_1.waveDir = [0]
        self.regularCIC_1.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)
        result1=  np.conj(np.transpose(np.array(np.loadtxt('./testData/regularCIC_1_test/waveAmpTime.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.regularCIC_1.waveAmpTime, result1))
        result2 =  [0.785398163397448]
        self.assertIsNone(np.testing.assert_allclose(self.regularCIC_1.w, result2))
               
        # bemFreq = np.array(np.loadtxt('irregular_1_test/bemFreq.txt')[0]
        # bemFreq = [5.19999512307279,0.0199999977946844]
        # wDepth = 'infinite'
        # rampTime = 100
        # dt = 0.1
        # maxIt = 2000
        # g = 9.81
        # rho = 1000
        # endTime = 200
        self.irregular_1.freqRange = []
        self.irregular_1.freqDisc = 'EqualEnergy'
        self.irregular_1.numFreq = 500
        self.irregular_1.phaseSeed = 1
        self.irregular_1.T = 8
        self.irregular_1.H = 2.5
        self.irregular_1.spectrumType = 'BS'
        self.irregular_1.waveDir = [0]
        self.irregular_1.waveSpread = [1]
        self.irregular_1.wavegauge1loc = [5,5]
        self.irregular_1.wavegauge2loc = [10,0]
        self.irregular_1.wavegauge3loc = [0,-10]
        self.irregular_1.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)
        result1 = np.conj(np.transpose(np.array(np.loadtxt('./testData/irregular_1_test/waveAmpTime.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.waveAmpTime, result1))
        result2 = np.conj(np.transpose(np.array(np.loadtxt('./testData/irregular_1_test/waveAmpTime1.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.waveAmpTime1, result2))
        result3 = np.conj(np.transpose(np.array(np.loadtxt('./testData/irregular_1_test/waveAmpTime2.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.waveAmpTime2, result3))
        result4 = np.conj(np.transpose(np.array(np.loadtxt('./testData/irregular_1_test/waveAmpTime3.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_1.waveAmpTime3, result4))
        
        bemFreq = [5.19999512307279,0.0199999977946844]
        wDepth = 150
        rampTime = 100
        dt = 0.1
        maxIt = 2000
        g = 9.81
        rho = 1000
        endTime = 200
        self.irregular_2.freqRange = [0.03,5]
        self.irregular_2.freqDisc = 'EqualEnergy'
        self.irregular_2.numFreq = 500
        self.irregular_2.phaseSeed = 1
        self.irregular_2.T = 8
        self.irregular_2.H = 2.5
        self.irregular_2.spectrumType = 'BS'
        self.irregular_2.waveDir = [0]
        self.irregular_2.waveSpread = [1]
        self.irregular_2.wavegauge1loc = [0,0]
        self.irregular_2.wavegauge2loc = [0,0]
        self.irregular_2.wavegauge3loc = [0,0]
        self.irregular_2.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)
        result = np.conj(np.transpose(np.array(np.loadtxt('./testData/irregular_2_test/waveAmpTime.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_2.waveAmpTime, result))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_2.Pw, 128939.416506839))
        

        # bemFreq = [5.19999512307279,0.0199999977946844]
        wDepth = 'infinite'
        # rampTime = 100
        # dt = 0.1
        # maxIt = 2000
        # g = 9.81
        # rho = 1000
        # endTime = 200
        self.irregular_3.freqRange = [0.01,6]
        self.irregular_3.freqDisc = 'EqualEnergy'
        self.irregular_3.numFreq = 500
        self.irregular_3.phaseSeed = 1
        self.irregular_3.T = 8
        self.irregular_3.H = 2.5
        self.irregular_3.spectrumType = 'BS'
        self.irregular_3.waveDir = [0]
        self.irregular_3.waveSpread = [1]
        self.irregular_3.wavegauge1loc = [0,0]
        self.irregular_3.wavegauge2loc = [0,0]
        self.irregular_3.wavegauge3loc = [0,0]
        self.irregular_3.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)
        result = np.conj(np.transpose(np.array(np.loadtxt('./testData/irregular_3_test/waveAmpTime.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_3.waveAmpTime, result))
        
        # bemFreq = np.array(np.loadtxt('irregular_1_test/bemFreq.txt')[0]
        # bemFreq = [5.19999512307279,0.0199999977946844]
        # wDepth = 'infinite'
        # rampTime = 100
        # dt = 0.1
        # maxIt = 2000
        # g = 9.81
        # rho = 1000
        # endTime = 200
        self.irregular_4.freqRange = []
        self.irregular_4.freqDisc = 'EqualEnergy'
        self.irregular_4.numFreq = 500
        self.irregular_4.phaseSeed = 1
        self.irregular_4.T = 8
        self.irregular_4.H = 2.5
        self.irregular_4.spectrumType = 'BS'
        self.irregular_4.waveDir = [0,90]
        self.irregular_4.waveSpread = [0,90]
        self.irregular_4.wavegauge1loc = [0,0]
        self.irregular_4.wavegauge2loc = [0,0]
        self.irregular_4.wavegauge3loc = [0,0]
        self.irregular_4.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)
        result = np.conj(np.transpose(np.array(np.loadtxt('./testData/irregular_4_test/waveAmpTime.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_4.waveAmpTime, result))
        
        # bemFreq = [5.19999512307279,0.0199999977946844]
        # wDepth = 'infinite'
        # rampTime = 100
        # dt = 0.1
        # maxIt = 2000
        # g = 9.81
        # rho = 1000
        # endTime = 200
        self.irregular_5.freqRange = []
        self.irregular_5.freqDisc = 'EqualEnergy'
        self.irregular_5.numFreq = 500
        self.irregular_5.phaseSeed = 1
        self.irregular_5.T = 8
        self.irregular_5.H = 2.5
        self.irregular_5.spectrumType = 'BS'
        self.irregular_5.waveDir = [30,60]
        self.irregular_5.waveSpread = [30,60]
        self.irregular_5.wavegauge1loc = [0,0]
        self.irregular_5.wavegauge2loc = [0,0]
        self.irregular_5.wavegauge3loc = [0,0]
        self.irregular_5.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)
        result = np.conj(np.transpose(np.array(np.loadtxt('./testData/irregular_5_test/waveAmpTime.txt'))))
        self.assertIsNone(np.testing.assert_allclose(self.irregular_5.waveAmpTime, result))
        
        # wType = 'spectrumImport'
        # bemFreq = [5.19999512307279,0.0199999977946844]
        # wDepth = 'infinite'
        # rampTime = 100
        # dt = 0.1
        # maxIt = 2000
        # g = 9.81
        # rho = 1000
        endTime = 400
        self.spectrumImport_1.freqRange = []
        self.spectrumImport_1.spectrumDataFile = "./testData/spectrumImport_1_test/spectrumData.txt"
        self.spectrumImport_1.phaseSeed = 1
        self.spectrumImport_1.waveDir = [0]
        self.spectrumImport_1.waveSpread = [1]
        self.spectrumImport_1.wavegauge1loc = [0,0]
        self.spectrumImport_1.wavegauge2loc = [0,0]
        self.spectrumImport_1.wavegauge3loc = [0,0]
        self.spectrumImport_1.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)
        result =  np.conj(np.transpose(np.loadtxt('./testData/spectrumImport_1_test/waveAmpTime.txt')))
        self.assertIsNone(np.testing.assert_allclose(self.spectrumImport_1.waveAmpTime, result))
        
        # bemFreq = [5.19999512307279,0.0199999977946844]
        # wDepth = 'infinite'
        # rampTime = 100
        # dt = 0.1
        maxIt = 4000
        # g = 9.81
        # rho = 1000
        # endTime = 400
        self.etaImport_1.etaDataFile = './testData/etaImport_1_test/etaData.txt'
        self.etaImport_1.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)
        result1 = np.conj(np.transpose(np.loadtxt("./testData/etaImport_1_test/waveAmpTime.txt")))
        self.assertIsNone(np.testing.assert_allclose(self.etaImport_1.waveAmpTime, result1))
        result2 = np.conj(np.transpose(np.loadtxt("./testData/etaImport_1_test/waveAmpTime1.txt")))
        self.assertIsNone(np.testing.assert_allclose(self.etaImport_1.waveAmpTime1, result2))
        
        # bemFreq = [5.19999512307279,0.0199999977946844]
        # wDepth = 'infinite'
        # rampTime = 100
        # dt = 0.1
        # maxIt = 4000
        # g = 9.81
        # rho = 1000
        # endTime = 400
        self.etaImport_2.etaDataFile = './testData/etaImport_1_test/etaData.mat'
        self.etaImport_2.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)
        result1 = np.conj(np.transpose(np.loadtxt("./testData/etaImport_1_test/waveAmpTime.txt")))
        self.assertIsNone(np.testing.assert_allclose(self.etaImport_2.waveAmpTime, result1))
        result2 = np.conj(np.transpose(np.loadtxt("./testData/etaImport_1_test/waveAmpTime1.txt")))
        self.assertIsNone(np.testing.assert_allclose(self.etaImport_2.waveAmpTime1, result2))
        
if __name__ == '__main__':
    unittest.main()