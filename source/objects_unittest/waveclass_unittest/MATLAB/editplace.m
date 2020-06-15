clear all,close all,clc

waves = waveClass('irregular')
bemFreq = [5.19999512307279,0.0199999977946844];
wDepth = 'infinite';
rampTime = 100;
dt = 0.1;
maxIt = 2000;
g = 9.81;
rho = 1000;
endTime = 200;
waves.freqRange = [];
waves.freqDisc = 'EqualEnergy';
waves.numFreq = 500;
waves.phaseSeed = 1;
waves.T = 8;
waves.H = 2.5;
waves.spectrumType = 'BS';
waves.waveDir = [0];
waves.waveSpread = [1];
% waves.wavegauge1loc = [0,0];
% waves.wavegauge2loc = [0,0];
% waves.wavegauge3loc = [0,0];
tic;
waves.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime);
toc;
% waves.plotEta(rampTime)
% hold on
% waves.plotEta()
% hold on
% waves.plotSpectrum()
