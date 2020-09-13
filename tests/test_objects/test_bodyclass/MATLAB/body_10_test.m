clc,clear all; close all;
%% Simulation Data
simu = simulationClass();               
% simu.simMechanicsFile = 'ellipsoid/ode4/Regular/ellipsoid.slx';   
simu.simMechanicsFile = 'ellipsoid/ellipsoid.slx';   
simu.mode = 'normal';                   % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
simu.explorer='off';                     % Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                     
simu.rampTime = 50;                        
simu.endTime=150;                       
simu.dt = 0.05;                         
simu.rho=1025;                                                                                                           
simu.nlHydro = 2;                       % Non-linear hydro on/off

%% Wave Information 
% Regular Waves 
waves = waveClass('regular');                 
waves.H = 4;                            
waves.T = 6;                

%% Body Data
body(1) = bodyClass('ellipsoid/ellipsoid.h5');
body(1).mass = 'equilibrium';           
body(1).momOfInertia = ...              
    [1.375264e6 1.375264e6 1.341721e6];      
body(1).geometryFile = 'ellipsoid/elipsoid.stl' ;    
body(1).viscDrag.cd=[1 0 1 0 1 0];
body(1).viscDrag.characteristicArea=[25 0 pi*5^2 0 pi*5^5 0];


%%


body(1).bodyNumber = 1;
% body(1).bodyTotal = 1;
body(1).readH5File();

simu.setupSim;
waves.waveSetup(body(1).hydroData.simulation_parameters.w, body(1).hydroData.simulation_parameters.water_depth, simu.rampTime, simu.dt, simu.maxIt, simu.g, simu.rho,  simu.endTime);



kk = 1;
body(kk).bodyGeo(body(kk).geometryFile)
body(kk).hydroForcePre(waves.w,waves.waveDir,simu.CIkt,simu.CTTime,waves.numFreq,simu.dt,...
    simu.rho,simu.g,waves.type,waves.waveAmpTime,kk,simu.numWecBodies,simu.ssCalc,simu.nlHydro,simu.b2b);


%%
body10_norm = body(1).bodyGeometry.norm
body10_area = body(1).bodyGeometry.area
body10_center = body(1).bodyGeometry.center
body10_vertex = body(1).bodyGeometry.vertex
body10_face = body(1).bodyGeometry.face
body10_numVertex = body(1).bodyGeometry.numVertex
body10_numFace = body(1).bodyGeometry.numFace
body10_linearHydroRestCoef = body(1).hydroForce.linearHydroRestCoef 
body10_visDrag = body(1).hydroForce.visDrag
body10_linearDamping = body(1).hydroForce.linearDamping
body10_userDefinedFe = body(1).hydroForce.userDefinedFe
body10_re = body(1).hydroForce.fExt.re
body10_im = body(1).hydroForce.fExt.im
body10_md = body(1).hydroForce.fExt.md
body10_fAddedMass = body(1).hydroForce.fAddedMass
body10_fDamping = body(1).hydroForce.fDamping
% body1_irkb = body(1).hydroForce.irkb
% body1_A = body(1).hydroForce.ssRadf.A
% body1_B = body(1).hydroForce.ssRadf.B
% body1_C = body(1).hydroForce.ssRadf.C
% body1_D = body(1).hydroForce.ssRadf.D
body10_totDOF = body(1).hydroForce.totDOF

CTTime = simu.CTTime
w = waves.w
waveAmpTime = waves.waveAmpTime
% 
mkdir body_10_test
save body_10_test/body1_norm.mat body10_norm
save body_10_test/body1_area.mat body10_area
save body_10_test/body1_center.mat body10_center
save body_10_test/body1_vertex.mat body10_vertex
save body_10_test/body1_face.mat body10_face
save body_10_test/body1_linearHydroRestCoef.mat  body10_linearHydroRestCoef
save body_10_test/body1_visDrag.mat  body10_visDrag 
save body_10_test/body1_linearDamping.mat body10_linearDamping
save body_10_test/body1_userDefinedFe.mat  body10_userDefinedFe 
save body_10_test/body1_re.mat  body10_re
save body_10_test/body1_im.mat  body10_im 
save body_10_test/body1_md.mat  body10_md
save body_10_test/body1_fAddedMass.mat body10_fAddedMass
save body_10_test/body1_fDamping.mat  body10_fDamping
% save body_10_test/body1_irkb.mat body1_irkb
% save body_10_test/body1_A.mat body1_A
% save body_10_test/body1_B.mat body1_B
% save body_10_test/body1_C.mat body1_C
% save body_10_test/body1_D.mat body1_D
save body_10_test/body1_totDOF.mat  body10_totDOF
% 
% save body_10_test/body2_linearHydroRestCoef.mat  body2_linearHydroRestCoef
% save body_10_test/body2_visDrag.mat  body2_visDrag 
% save body_10_test/body2_linearDamping.mat body2_linearDamping
% save body_10_test/body2_userDefinedFe.mat  body2_userDefinedFe 
% save body_10_test/body2_re.mat  body2_re
% save body_10_test/body2_im.mat  body2_im 
% save body_10_test/body2_md.mat  body2_md 
% save body_10_test/body2_fAddedMass.mat body2_fAddedMass
% save body_10_test/body2_fDamping.mat  body2_fDamping
% save body_10_test/body2_irkb.mat body2_irkb
% save body_10_test/body2_A.mat body2_A
% save body_10_test/body2_B.mat body2_B
% save body_10_test/body2_C.mat body2_C
% save body_10_test/body2_D.mat body2_D
% % save body_10_test/body2_totDOF.mat  body2_totDOF
% 
save body_10_test/CTTime.mat CTTime
save body_10_test/w.mat w
save body_10_test/waveAmpTime.mat waveAmpTime
