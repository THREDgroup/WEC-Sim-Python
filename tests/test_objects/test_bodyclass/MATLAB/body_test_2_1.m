%% Simulation Data
simu = simulationClass();                   % Initialize Simulation Class
simu.simMechanicsFile = 'RM3.slx';          % Specify Simulink Model File
simu.rampTime = 100;                        % Wave Ramp Time [s]
simu.endTime=200;                           % Simulation End Time [s]
simu.dt = 0.1;                              % Simulation Time-Step [s]
simu.explorer = 'off';
simu.CITime = 20;

%% Wave Information
% Irregular Waves using BS Spectrum with State Space Calculation
waves = waveClass('irregular');             % Initialize Wave Class and Specify Type
waves.H = 2.5;                              % Significant Wave Height [m]
waves.T = 8;                                % Peak Period [s]
waves.spectrumType = 'BS';                  % Specify Wave Spectrum Type
waves.phaseSeed = 1;

%% Body Data
% Float
body(1) = bodyClass('rm3.h5');
body(1).geometryFile = 'float.stl'; 
body(1).mass = 'equilibrium';
body(1).momOfInertia = [20907301 21306090.66 37085481.11];   

% Spar/Plate
body(2) = bodyClass('rm3.h5');
body(2).geometryFile = 'plate.stl';
body(2).mass = 'equilibrium';
body(2).momOfInertia = [94419614.57 94407091.24 28542224.82];

body(1).bodyNumber = 1;
body(1).readH5File();
simu.setupSim;
waves.waveSetup(body(1).hydroData.simulation_parameters.w, body(1).hydroData.simulation_parameters.water_depth, simu.rampTime, simu.dt, simu.maxIt, simu.g, simu.rho,  simu.endTime);
kk = 1;
body(kk).hydroForcePre(waves.w,waves.waveDir,simu.CIkt,simu.CTTime,waves.numFreq,simu.dt,...
    simu.rho,simu.g,waves.type,waves.waveAmpTime,kk,simu.numWecBodies,simu.ssCalc,simu.nlHydro,simu.b2b);
