function [dr,r,d2r,d3r]=CRKprof_circarc2(t,crk_prams);

t=t';
%%x\in[a,b]:
% alf=crk_prams{1};
% k=pi*alf;

th = pi*(t+1)/2;
rho = sin(th).^3 + cos(th).^3;

x=rho.*cos(th); 
y=rho.*sin(th); 

r=[x;y];
dr=[0 0];%[0 -k;k 0]*r;
d2r=[0 0];%-k^2*r;
d3r=[0 0];%-k^2*dr;