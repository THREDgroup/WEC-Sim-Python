clc
clear all
close all
n=3;			%No. of points required for describing the shape of the rotating curve (generatrix) of the surface of revolution (submerged part only). 
Radius=5			%Radius of the cylinder
Draft=10			%Height of the submerged part
r=[Radius Radius 0];		%Array of radial coordinates
z=[0 -Draft -Draft];		%Array of vertical coordinates


[Mass,Inertia,KH,XB,YB,ZB]=axiMesh(r,z,n)	%Call the function axiMesh.m

M=Inertia; 		%Mass Matrix

%save('Mesh outputs','KH','M')
