clc,clear all; close all;
%% Simulation Data
simu = simulationClass();               
simu.simMechanicsFile = 'rm3/RM3.slx';      
simu.solver = 'ode4';                   
simu.explorer='off';                                       
simu.startTime = 0;                     
simu.rampTime = 100;                       
simu.endTime=400;                       
simu.dt = 0.1; 							
simu.b2b = 0;                   	% Turn B2B interactions 'off'   

%% Wave Information 
% Regular Waves  
waves = waveClass('regularCIC');       % Regular CIC    
waves.H = 2.5;                          
waves.T = 8;                            

%% Body Data
% Float
body(1) = bodyClass('rm3/rm3.h5');      
body(1).geometryFile = 'rm3/float.stl';   
body(1).mass = 'equilibrium';                   
body(1).momOfInertia = [20907301 21306090.66 37085481.11];       

% Spar/Plate
body(2) = bodyClass('rm3/rm3.h5'); 
body(2).mass = 'equilibrium';                   
body(2).momOfInertia = [94419614.57 94407091.24 28542224.82];
body(2).geometryFile = 'rm3/plate.stl'; 

%%
body(1).bodyNumber = 1;
body(1).bodyTotal = 2;
body(1).readH5File();
body(2).bodyTotal = 2;
body(2).bodyNumber = 2;
body(2).readH5File();
simu.setupSim;
waves.waveSetup(body(1).hydroData.simulation_parameters.w, body(1).hydroData.simulation_parameters.water_depth, simu.rampTime, simu.dt, simu.maxIt, simu.g, simu.rho,  simu.endTime);

kk = 1;
body(kk).hydroForcePre(waves.w,waves.waveDir,simu.CIkt,simu.CTTime,waves.numFreq,simu.dt,...
    simu.rho,simu.g,waves.type,waves.waveAmpTime,kk,simu.numWecBodies,simu.ssCalc,simu.nlHydro,simu.b2b);
kk = 2;
body(kk).hydroForcePre(waves.w,waves.waveDir,simu.CIkt,simu.CTTime,waves.numFreq,simu.dt,...
    simu.rho,simu.g,waves.type,waves.waveAmpTime,kk,simu.numWecBodies,simu.ssCalc,simu.nlHydro,simu.b2b);

%%
body1_linearHydroRestCoef = body(1).hydroForce.linearHydroRestCoef 
body1_visDrag = body(1).hydroForce.visDrag
body1_linearDamping = body(1).hydroForce.linearDamping
body1_userDefinedFe = body(1).hydroForce.userDefinedFe
body1_re = body(1).hydroForce.fExt.re
body1_im = body(1).hydroForce.fExt.im
body1_md = body(1).hydroForce.fExt.md
body1_fAddedMass = body(1).hydroForce.fAddedMass
body1_fDamping = body(1).hydroForce.fDamping
body1_irkb = body(1).hydroForce.irkb
% body1_totDOF = body(1).hydroForce.totDOF


body2_linearHydroRestCoef = body(2).hydroForce.linearHydroRestCoef 
body2_visDrag = body(2).hydroForce.visDrag
body2_linearDamping = body(2).hydroForce.linearDamping
body2_userDefinedFe = body(2).hydroForce.userDefinedFe
body2_re = body(2).hydroForce.fExt.re
body2_im = body(2).hydroForce.fExt.im
body2_md = body(2).hydroForce.fExt.md
body2_fAddedMass = body(2).hydroForce.fAddedMass
body2_fDamping = body(2).hydroForce.fDamping
body2_irkb = body(2).hydroForce.irkb
% body2_totDOF = body(2).hydroForce.totDOF

CTTime = simu.CTTime
w = waves.w
waveAmpTime = waves.waveAmpTime
% 
% mkdir body_6_test
% save body_6_test/body1_linearHydroRestCoef.mat  body1_linearHydroRestCoef
% save body_6_test/body1_visDrag.mat  body1_visDrag 
% save body_6_test/body1_linearDamping.mat body1_linearDamping
% save body_6_test/body1_userDefinedFe.mat  body1_userDefinedFe 
% save body_6_test/body1_re.mat  body1_re
% save body_6_test/body1_im.mat  body1_im 
% save body_6_test/body1_md.mat  body1_md
% save body_6_test/body1_fAddedMass.mat body1_fAddedMass
% save body_6_test/body1_fDamping.mat  body1_fDamping
% save body_6_test/body1_irkb.mat body1_irkb
% % save body_6_test/body1_totDOF.mat  body1_totDOF
% 
% save body_6_test/body2_linearHydroRestCoef.mat  body2_linearHydroRestCoef
% save body_6_test/body2_visDrag.mat  body2_visDrag 
% save body_6_test/body2_linearDamping.mat body2_linearDamping
% save body_6_test/body2_userDefinedFe.mat  body2_userDefinedFe 
% save body_6_test/body2_re.mat  body2_re
% save body_6_test/body2_im.mat  body2_im 
% save body_6_test/body2_md.mat  body2_md 
% save body_6_test/body2_fAddedMass.mat body2_fAddedMass
% save body_6_test/body2_fDamping.mat  body2_fDamping
% save body_6_test/body2_irkb.mat body2_irkb
% % save body_6_test/body2_totDOF.mat  body2_totDOF
% 
% save body_6_test/CTTime.mat CTTime
% save body_6_test/w.mat w
% save body_6_test/waveAmpTime.mat waveAmpTime


