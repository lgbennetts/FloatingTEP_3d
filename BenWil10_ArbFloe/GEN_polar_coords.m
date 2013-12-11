function [r,theta]=GEN_polar_coords(x,y)

z=x+1i*y;
r=abs(z);
theta=angle(z);