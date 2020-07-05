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
simu.b2b = 0;                       % Turn B2B interactions 'off'   
simu.ssCalc = 1;

%% Wave Information 
% Regular Waves  
waves = waveClass('regularCIC');    % Regular CIC           
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
% kk = 1;
% body(kk).hydroForcePre(waves.w,waves.waveDir,simu.CIkt,simu.CTTime,waves.numFreq,simu.dt,...
%     simu.rho,simu.g,waves.type,waves.waveAmpTime,kk,simu.numWecBodies,simu.ssCalc,simu.nlHydro,simu.b2b);
% kk = 2;
% body(kk).hydroForcePre(waves.w,waves.waveDir,simu.CIkt,simu.CTTime,waves.numFreq,simu.dt,...
%     simu.rho,simu.g,waves.type,waves.waveAmpTime,kk,simu.numWecBodies,simu.ssCalc,simu.nlHydro,simu.b2b);

B2B = simu.b2b;
rho = simu.rho
ssCalc = simu.ssCalc
CIkt = simu.CIkt
CTTime = simu.CTTime

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
                    arraySize
                    Af(1:arraySize,1:arraySize) = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.A.all(ii,jj,1:arraySize,1:arraySize)
                    Bf(1:arraySize,jInd)        = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.B.all(ii,jj,1:arraySize,1);
                    Cf(ii,1:arraySize)          = body(1).hydroData.hydro_coeffs.radiation_damping.state_space.C.all(ii,jj,1,1:arraySize);
                    
                else
                    AAA = size(Af,1);
                    AAA1 = size(Af,2);
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
