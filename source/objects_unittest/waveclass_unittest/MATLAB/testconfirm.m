clc;clear all; close all;
wType = 'etaImport';
bemFreq = [5.19999512307279,0.0199999977946844];
wDepth = 'infinite';
rampTime = 100;
dt = 0.1;
maxIt = 4000;
g = 9.81;
rho = 1000;
endTime = 400;
freqRange = [];
numFreq = 0;
phaseSeed = 0;
waveDir = 0;
waveSpread = 1;
wavegauge1loc = [0,0];
wavegauge2loc = [0,0];
wavegauge3loc = [0,0];
wF = [];
dw = [];
waveAmpTime = [];
waveAmpTime1 = [];
waveAmpTime2 = [];
waveAmpTime3 = [];
etaDataFile = 'etaData.mat';

data = importdata(etaDataFile) ;    % Import time-series
t = [0:dt:endTime]';      % WEC-Sim simulation time [s]
waveAmpTime = zeros(maxIt+1,2);
maxRampIT=round(rampTime/dt);
data_t = data(:,1)';                    % Data Time [s]
data_x = data(:,2)';                    % Wave Surface Elevation [m]
waveAmpTime(:,1) = t;
waveAmpTime(:,2) = interp1(data_t,data_x,t);
if rampTime~=0
    for i=1:maxRampIT
        waveAmpTime(i,2) = waveAmpTime(i,2)*(1+cos(pi+pi*(i-1)/maxRampIT))/2;
    end
end
waveAmpTime1        = zeros(maxIt+1,2);
waveAmpTime1(:,1)   = [0:maxIt]*dt;
waveAmpTime2        = zeros(maxIt+1,2);
waveAmpTime2(:,1)   = [0:maxIt]*dt;
waveAmpTime3        = zeros(maxIt+1,2);
waveAmpTime3(:,1)   = [0:maxIt]*dt;