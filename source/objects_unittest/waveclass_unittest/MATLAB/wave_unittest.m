clear all,close all,clc
%test setWavePhase
% waves = waveClass('irregular')
% waves.phaseSeed = 1
% waves.freqDisc = 'EqualEnergy'
% waves.numFreq = 4
% waves.waveDir = [0;90]
% waves.setWavePhase()
% waves.phase

%test waveNumber
% waves = waveClass('irregular')
% g = 9.81
% waves.w = [0.566624163384781;
% 0.582081268838611;
% 0.590203501195047;
% 0.600522051484601];
% waves.deepWaterWave = 1;
% waves.waterDepth = 200;
% waves.waveNumber(g);
% waves.k

% waves = waveClass('irregular')
% g = 9.81
% waves.w = [0.525743641856087;
%     0.541511547017433;
%     0.551933697209493;
%     0.559973049643924;
%     0.566624163384781];
% waves.deepWaterWave = 0;
% waves.waterDepth = 150;
% waves.waveNumber(g);
% waves.k

%test irregWaveSpectrum
% waves = waveClass('irregular')
% g = 9.81
% rho = 1000
% waves.w = [0.525743641856087;
% 0.541511547017433;
% 0.551933697209493;
% 0.559973049643924;
% 0.566624163384781;
% 0.572373957973840;
% 0.577471073177114;
% 0.582081268838611;
% 0.586297784870588;
% 0.590203501195047;
% 0.593850217763243;
% 0.597279374536177;
% 0.600522051484601;
% 0.603609328579267;
% 0.606551565810425]
% waves.T = 8
% waves.H = 2.5
% waves.spectrumType = 'BS'
% waves.freqDisc = 'EqualEnergy'
% waves.deepWaterWave = 1
% waves.numFreq = 4
% waves.dw = [0.505743644061402;
% 0.0157679051613465;
% 0.0104221501920596;
% 0.00803935243443155;
% 0.00665111374085714;
% 0.00574979458905867;
% 0.00509711520327372;
% 0.00461019566149745;
% 0.00421651603197637;
% 0.00390571632445969;
% 0.00364671656819582;
% 0.00342915677293409;
% 0.00324267694842406;
% 0.00308727709466572;
% 0.00294223723115805]
% waves.irregWaveSpectrum(g,rho)
% waves.w
% waves.dw
% waves.S
% waves.A

% waves = waveClass('irregular')
% g = 9.81
% rho = 1000
% waves.w = [0.525743641856087;
% 0.541511547017433;
% 0.551933697209493;
% 0.559973049643924];
% waves.T = 8
% waves.H = 2.5
% waves.spectrumType = 'PM'
% waves.freqDisc = 'Traditional'
% waves.deepWaterWave = 1
% waves.numFreq = 4
% waves.dw = [0.505743644061402;
% 0.0157679051613465;
% 0.0104221501920596;
% 0.00803935243443155];
% waves.irregWaveSpectrum(g,rho)

% waves = waveClass('irregular')
% g = 9.81
% rho = 1000
% waves.w = [0.525743641856087;
% 0.541511547017433;
% 0.551933697209493;
% 0.559973049643924];
% waves.T = 8
% waves.H = 2.5
% waves.spectrumType = 'JS'
% waves.freqDisc = 'Traditional'
% waves.deepWaterWave = 0
% waves.waterDepth = 100
% waves.numFreq = 4
% waves.dw = [0.505743644061402;
% 0.0157679051613465;
% 0.0104221501920596;
% 0.00803935243443155];
% waves.irregWaveSpectrum(g,rho)

%test waveElevReg
% waves = waveClass('regular')
% rampTime = 0.2;
% dt = 0.1;
% maxIt = 4;
% waves.waveDir = 0;
% waves.A = [1.25];
% waves.waveSpread = 1;
% waves.w = [0.785398163397448];
% waves.k = [0.0628797426165224];
% waves.wavegauge1loc = [0,0];
% waves.wavegauge2loc = [0,0];
% waves.wavegauge3loc = [0,0];
% waves.waveElevReg(rampTime,dt,maxIt);
% waves.waveAmpTime;

%test wavePowerReg
% waves = waveClass('regular');
% g = 9.81;
% rho = 1000;
% waves.A = [1.25];
% waves.deepWaterWave = 0;
% waves.T = 8;
% waves.k = [0.0628797426165224];
% waves.waterDepth = 200
% if waves.deepWaterWave == 1
%     % Deepwater Approximation
%     waves.Pw = 1/(8*pi)*rho*g^(2)*(waves.A).^(2).*waves.T;               
% else
%     % Full Wave Power Equation
%     waves.Pw = rho*g*(waves.A).^(2)/4*sqrt(g./waves.k.*tanh(waves.k.*waves.waterDepth))*(1+2*waves.k.*waves.waterDepth./sinh(waves.k.*waves.waterDepth));
% end
% waves.Pw


%test waveElevIrreg
waves = waveClass('irregular')
rampTime = 0.2;
dt = 0.1;
maxIt = 4;
df = [0.0408805215286940;0.0154571054538299;0.00812223235643605;0.0103185502895540];
waves.waveDir = 0;
waves.A = [0.253488803947639;0.354610720753343;0.413569577232088;0.492917523052480];
waves.waveSpread = 1;
waves.w = [0.566624163384781;0.582081268838611;0.590203501195047;0.600522051484601];
waves.phase = [2.62022653271779;4.52593227359735;0.000718638171852741;1.89961157824218];
waves.k = [0.0327281286984203;0.0345380839482943;0.0355086822449431;0.0367611349968679];
waves.wavegauge1loc = [0,0];
waves.wavegauge2loc = [0,0];
waves.wavegauge3loc = [0,0];
waves.waveElevIrreg(rampTime,dt,maxIt,df);
waves.waveAmpTime;

%test waveSetup
% waves = waveClass('noWave');
% bemFreq = [5.19999512307279,0.0199999977946844];
% wDepth = 'infinite';
% rampTime = 100;
% dt = 0.1;
% maxIt = 2000;
% g = 9.81;
% rho = 1000;
% endTime = 200;
% waves.T = 8;
% waves.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)

% waves = waveClass('noWaveCIC');
% bemFreq = [5.19999512307279,0.0199999977946844];
% wDepth = 'infinite';
% rampTime = 100;
% dt = 0.1;
% maxIt = 2000;
% g = 9.81;
% rho = 1000;
% endTime = 200;
% waves.waveSetup(bemFreq,wDepth,rampTime,dt,maxIt,g, rho, endTime)


