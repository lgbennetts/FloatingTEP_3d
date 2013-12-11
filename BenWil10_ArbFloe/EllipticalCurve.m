function [dr,r,d2r,d3r]=EllipticalCurve(t,~)

t=t';
%%x\in[a,b]:
alf=1; %crk_prams{1};
k=pi*alf;

x=cos(k*t);%+pi/4
y=sin(k*t);%+pi/4

r=[x;y];
dr=[0 -k;k 0]*r;
d2r=-k^2*r;
d3r=-k^2*dr;

return