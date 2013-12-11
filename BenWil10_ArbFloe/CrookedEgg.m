function [dr,r,d2r,d3r]=CrookedEgg(t,~)

t=t';

th = pi*(t+1)/2; % - Tim: this begins the curve at (x,y)=(1,0) (t=0) & doesn't work (not symm)
% th = pi*t/2 + 3*pi/4; % - Tim: this begins the curve at (x,y)=(0,0) (t=0) & does work (symm)
rho = sin(th).^3 + cos(th).^3;

x=rho.*cos(th)- 0.35;  %
y=rho.*sin(th)- 0.35;  %

fac = pi/2;

rho_dth = (3/2)*sin(2*th).*(sin(th)-cos(th));
rho_dth2 = (3/2)*( sin(2*th).*(sin(th)+cos(th)) ...
    ...
    + 2*cos(2*th).*(sin(th)-cos(th)));
rho_dth3 = (3/2)*(5*sin(2*th).*(cos(th)-sin(th)) ...
    ...
    + 4*cos(2*th).*(cos(th)+sin(th)));

r=[x;y];
dr=[cos(th).*rho_dth - sin(th).*rho; sin(th).*rho_dth + cos(th).*rho];
d2r=[cos(th).*rho_dth2 - 2*sin(th).*rho_dth - cos(th).*rho;...
    sin(th).*rho_dth2 + 2*cos(th).*rho_dth - sin(th).*rho];
d3r=[cos(th).*rho_dth3 - 3*sin(th).*rho_dth2 - 3*cos(th).*rho_dth + sin(th).*rho;...
    sin(th).*rho_dth3 + 3*cos(th).*rho_dth2 - 3*sin(th).*rho_dth - cos(th).*rho];

dr = fac*dr;
d2r = fac^2*d2r;
d3r = fac^3*d3r;
