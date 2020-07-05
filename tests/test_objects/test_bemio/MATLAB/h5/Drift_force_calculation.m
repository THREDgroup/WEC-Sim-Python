% AIM: Calculation of the drift force on the body thanks to the "Far Field
% Formula" (Maruo,1960 and Newmann,1967), with implementation into NEMOH
% based on the PhD thesis of G.Delhommeau, "Les problemes de diffraction-radiation et de resistance de vagues: Etude theoriaue et resolution numerique par la methode des singularites", 1987.
% 
% INPUTS :
% - RAO : Displacement RAO of our body incident angle of waves= dir. Matrix (nb frequencies*num_DoF*nBodies)
% - w : frequency vector (rad/s)
% - depth : Water depth (m)
% - dir : direction of the incident wave field (degrees)
% - plotflag: 1= plot the figure, 0= no plot.
% - ampl_wave : waves amplitude (m).
% - nBodies: Number of Bodies.
% - num_DoF : Number  of degrees of freedom
% 
% OUTPUTS :
% - Fdrift_x : Matrix (1*nb frequencies)
% - Fdrift_y : Matrix (1*nb frequencies)

function [Fdrift_x,Fdrift_y]=Drift_force_calculation(depth,w,RAO,ampl_wave,dir,nBodies,num_DoF,plotflag)
fid1=fopen('ID.dat');
line=fgetl(fid1);
rep=fscanf(fid1,'%s',1);
fclose(fid1); 
nw=length(w);

%================= READING KOCHIN FILES ===================%
clear Kochin_BVP x theta i H
for p=1:nw
    for i=1:(num_DoF*nBodies+1)
        clear Kochin
x=(num_DoF*nBodies+1)*(p-1)+i; 
        switch  numel(num2str(x))  
            case 1
                filename=['\Results\Kochin.    ',num2str(x),'.dat'];                  
            case 2
                filename=['\Results\Kochin.   ',num2str(x),'.dat'];
            case 3
                filename=['\Results\Kochin.  ',num2str(x),'.dat'];
            case 4
                filename=['\Results\Kochin. ',num2str(x),'.dat'];
            case 5
                filename=['\Results\Kochin.',num2str(x),'.dat'];
        end
        fid=fopen([rep,filename],'r');
        Kochin= fscanf(fid,'%f');

        for ntheta=1:size(Kochin,1)/3
         theta(ntheta)= Kochin(3*(ntheta-1)+1);
         Kochin_BVP(ntheta,1,x)= Kochin(3*(ntheta-1)+2);
         Kochin_BVP(ntheta,2,x)= Kochin(3*(ntheta-1)+3);        
        end
        status=fclose(fid);
    end
end
%======================END======================%

%------Initialisation-----
first_constant= zeros(1,nw);
second_constant= zeros(1,nw);
Fdrift_x=zeros(1,nw);
Fdrift_y=zeros(1,nw);
H=zeros(ntheta,2,nw);
rho=1025;
%--------------- CALCULATION-----------------------%
Kochin_BVP_complex(:,:)=Kochin_BVP(:,1,:).*exp(1i*Kochin_BVP(:,2,:));
for p=1:nw  
    clear m0 k0
    m0=wave_number(w(p)/(2*pi),depth);
    k0=(w(p)^2)/9.81;
first_constant(p)=-2*pi*ampl_wave*rho*w(p);
second_constant(p)=-(8*pi*rho*m0*(k0*depth)^2)/(depth*(m0^2*depth^2-k0^2*depth^2+k0*depth));

local1=zeros(ntheta,7);
local1(:,1)=ampl_wave*Kochin_BVP_complex(:,(num_DoF*nBodies+1)*(p-1)+1)*exp(1i*pi/2); % Diffraction term

for i=2:(num_DoF*nBodies+1)% sum of the radiation terms times velocities RAOs
x=(num_DoF*nBodies+1)*(p-1)+i;
local1(:,i)=ampl_wave*(1i*w(p))*(RAO(p,i-1)*exp(-1i*pi/2))*(Kochin_BVP_complex(:,x)); 
end
H(:,p)=sum(local1,2); % H= Kochin function per frequency
H_real(:,p)=real(H(:,p));
H_imag(:,p)=imag(H(:,p));
end
dir = pi/180*dir; % conversion radians to degrees
ind_beta=find(abs(theta-dir)==min(abs(theta-dir))); % ind_beta used for determining the appropriate angle of H(dir)
for p=1:nw
  Fdrift_x(p)=first_constant(p)*cos(dir)*imag(H(ind_beta,p)) + second_constant(p)*imag(trapz(theta,H_real(:,p).*imag(conj(H(:,p))).*cos(theta')));
  Fdrift_y(p)=first_constant(p)*sin(dir)*imag(H(ind_beta,p)) + second_constant(p)*imag(trapz(theta,H_real(:,p).*imag(conj(H(:,p))).*sin(theta')));
end

%-----------------PLOT-------------------------------%
if plotflag==1
figure
plot(2*pi./w,Fdrift_x,'g-o',2*pi./w,Fdrift_y,'r-o','LineWidth',2.0)
xlabel('T (s)') 
grid ON
ylabel('(N)' )
legend(' F Drift X',' F Drift Y')
title('Horizontal Drift Forces')
end
end