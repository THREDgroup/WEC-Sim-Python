clear all
simu = simulationClass();               % Initialize simulationClass
simu.simMechanicsFile = '.slx';         % Simulink Model File
simu.startTime = 477;                     % Simulation Start Time [s]
simu.rampTime = 1750;                   	% Wave Ramp Time [s]
simu.endTime=6920;                       % Simulation End Time [s]
simu.dt = 0.5; 
simu.setupSim()

simu_CTTime = simu.CTTime;
simu_time = simu.time;

% mkdir('simu_1')
% 
% save simu_1/simulation_CTTime.mat  simu_CTTime
% save simu_1/simulation_time.mat  simu_time 
x1 = 477:0.3:6920;
x2 = 477:0.4:6920;
x3 = 477:0.6:6920;
x4 = 477:0.7:6920;
x5 = 477:0.9:6920;
x6 = 78:0.3:6922;