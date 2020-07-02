clc;clear all;close all;
%% Simulation Data
simu = simulationClass();               % Initialize Simulation Class
simu.simMechanicsFile = 'OSWEC/OSWEC.slx';    % Specify Simulink Model File
simu.mode = 'normal';                   % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
simu.explorer='on';                     % Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                     % Simulation Start Time [s]
simu.rampTime = 100;                    % Wave Ramp Time [s]
simu.endTime=400;                       % Simulation End Time [s]        
simu.solver = 'ode4';                   % simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu.dt = 0.1;                          % Simulation Time-Step [s]
simu.CITime = 30;                       % Specify CI Time [s]

%% Wave Information
% % noWaveCIC, no waves with radiation CIC  
% waves = waveClass('noWaveCIC');       % Initialize Wave Class and Specify Type  

% % Regular Waves 
% waves = waveClass('regular');           % Initialize Wave Class and Specify Type                                 
% waves.H = 2.5;                          % Wave Height [m]
% waves.T = 8;                            % Wave Period [s]

% Irregular Waves using PM Spectrum with Directionality 
waves = waveClass('irregular');         % Initialize Wave Class and Specify Type
waves.H = 2.5;                          % Significant Wave Height [m]
waves.T = 8;                            % Peak Period [s]
waves.spectrumType = 'PM';              % Specify Spectrum Type
waves.waveDir = [0,30,90];              % Wave Directionality [deg]
waves.waveSpread = [0.1,0.2,0.7];       % Wave Directional Spreading [%}

% % Irregular Waves with imported spectrum
% waves = waveClass('spectrumImport');        % Create the Wave Variable and Specify Type
% waves.spectrumDataFile = 'spectrumData.mat';  %Name of User-Defined Spectrum File [:,2] = [f, Sf]

% % Waves with imported wave elevation time-history  
% waves = waveClass('etaImport');         % Create the Wave Variable and Specify Type
% waves.etaDataFile = 'etaData.mat'; % Name of User-Defined Time-Series File [:,2] = [time, eta]


%% Body Data
% Flap
body(1) = bodyClass('OSWEC/oswec.h5');      % Initialize bodyClass for Flap
body(1).geometryFile = 'OSWEC/flap.stl';     % Geometry File
body(1).mass = 127000;                          % User-Defined mass [kg]
body(1).momOfInertia = [1.85e6 1.85e6 1.85e6];  % Moment of Inertia [kg-m^2]

% Base
body(2) = bodyClass('OSWEC/oswec.h5');      % Initialize bodyClass for Base
body(2).geometryFile = 'OSWEC/base.stl';     % Geometry File
body(2).mass = 'fixed';                         % Creates Fixed Body

body(1).bodyNumber = 1;
body(1).readH5File();
simu.setupSim;
waves.waveSetup(body(1).hydroData.simulation_parameters.w, body(1).hydroData.simulation_parameters.water_depth, simu.rampTime, simu.dt, simu.maxIt, simu.g, simu.rho,  simu.endTime);
kk = 1;
% body(kk).hydroForcePre(waves.w,waves.waveDir,simu.CIkt,simu.CTTime,waves.numFreq,simu.dt,...
%     simu.rho,simu.g,waves.type,waves.waveAmpTime,kk,simu.numWecBodies,simu.ssCalc,simu.nlHydro,simu.b2b);
rho = simu.rho;
g = simu.g;
numFreq = waves.numFreq;
waveDir = waves.waveDir;
w = waves.w;
body(1).irrExcitation(w,numFreq,waveDir,rho,g);
B2B = simu.b2b;
CTTime = simu.CTTime;
ssCalc = simu.ssCalc;
nDOF = body(1).dof;
if B2B == 1
    LDOF = body(1).bodyTotal*6;
else
    LDOF = body(1).dof;
end
% Convolution integral formulation
if B2B == 1
    body(1).hydroForce.fAddedMass=body(1).hydroData.hydro_coeffs.added_mass.inf_freq .*rho;
else
    body(1).hydroForce.fAddedMass=body(1).hydroData.hydro_coeffs.added_mass.inf_freq(:,body(1).dof_start:body(1).dof_end) .*rho;
end
% Radition IRF
body(1).hydroForce.fDamping=zeros(nDOF,LDOF);
irfk = body(1).hydroData.hydro_coeffs.radiation_damping.impulse_response_fun.K  .*rho;
irft = body(1).hydroData.hydro_coeffs.radiation_damping.impulse_response_fun.t;
%body(1).hydroForce.irkb=zeros(CIkt,6,lenJ);
if B2B == 1
    for ii=1:nDOF
        for jj=1:LDOF
            body(1).hydroForce.irkb(:,ii,jj) = interp1(irft,squeeze(irfk(ii,jj,:)),CTTime,'spline');
        end
    end
else
    for ii=1:nDOF
        for jj=1:LDOF
            jjj = body(1).dof_start-1+jj;
            body(1).hydroForce.irkb(:,ii,jj) = interp1(irft,squeeze(irfk(ii,jjj,:)),CTTime,'spline');
        end
    end
end
% State Space Formulation
if ssCalc == 1
    if B2B == 1
        for ii = 1:nDOF
            for jj = 1:LDOF
                arraySize = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.it(ii,jj);
                if ii == 1 && jj == 1 % Begin construction of combined state, input, and output matrices
                    Af(1:arraySize,1:arraySize) = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.A.all(ii,jj,1:arraySize,1:arraySize);
                    Bf(1:arraySize,jj)        = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.B.all(ii,jj,1:arraySize,1);
                    Cf(ii,1:arraySize)          = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.C.all(ii,jj,1,1:arraySize);
                else
                    Af(size(Af,1)+1:size(Af,1)+arraySize,size(Af,2)+1:size(Af,2)+arraySize)     = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.A.all(ii,jj,1:arraySize,1:arraySize);
                    Bf(size(Bf,1)+1:size(Bf,1)+arraySize,jj) = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.B.all(ii,jj,1:arraySize,1);
                    Cf(ii,size(Cf,2)+1:size(Cf,2)+arraySize)   = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.C.all(ii,jj,1,1:arraySize);
                end
            end
        end
        body(1).hydroForce.ssRadf.D = zeros(nDOF,LDOF);
    else
        for ii = 1:nDOF
            for jj = body(1).dof_start:body(1).dof_end
                jInd = jj-body(1).dof_start+1;
                arraySize = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.it(ii,jj);
                if ii == 1 && jInd == 1 % Begin construction of combined state, input, and output matrices
                    Af(1:arraySize,1:arraySize) = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.A.all(ii,jj,1:arraySize,1:arraySize);
                    Bf(1:arraySize,jInd)        = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.B.all(ii,jj,1:arraySize,1);
                    Cf(ii,1:arraySize)          = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.C.all(ii,jj,1,1:arraySize);
                else
                    Af(size(Af,1)+1:size(Af,1)+arraySize,size(Af,2)+1:size(Af,2)+arraySize) = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.A.all(ii,jj,1:arraySize,1:arraySize);
                    Bf(size(Bf,1)+1:size(Bf,1)+arraySize,jInd) = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.B.all(ii,jj,1:arraySize,1);
                    Cf(ii,size(Cf,2)+1:size(Cf,2)+arraySize)   = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.C.all(ii,jj,1,1:arraySize);
                end
            end
        end
        body(1).hydroForce.ssRadf.D = zeros(nDOF,nDOF);
    end
    body(1).hydroForce.ssRadf.A = Af;
    body(1).hydroForce.ssRadf.B = Bf;
    body(1).hydroForce.ssRadf.C = Cf .*rho;
end
