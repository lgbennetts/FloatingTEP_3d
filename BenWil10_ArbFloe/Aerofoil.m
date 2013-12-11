function [dr,r,d2r,d3r]=Aerofoil(t,crk_prams);

k=.5;

if nargin==0
  np=100;
  t=-1+(0:np)'*2/np;
  %%
  phi=.8*(pi/2);
else
  phi=crk_prams{1}*pi/2;
end
s=k*(t+1);
n0=length(t);

%% want to scale perimeter to approx 1;
n1=200;
s=[s;(0:n1)'/n1];

if phi<=pi/4;
  a=1;
  b=tan(phi);
else
  a=cot(phi);
  b=1;
end

%%
x=a*s.*(s-1);
dx=a*k*(2*s-1);
d2x=2*a*k^2+0*s;
d3x=0*s;
%%
y=b*s.*(s-1).^2;
dy=k*b*(s-1).^2 + 2*k*b*s.*(s-1);
d2y=2*k^2*b*(s-1)+2*k^2*b*(2*s-1);
d3y=6*k^3*b+0*s;
%%
r=[x,y]';
dr=[dx,dy]';
d2r=[d2x,d2y]';
d3r=[d3x,d3y]';

%% SCALE SO PERIMETER IS APPROX 1:
L=GEN_arclength(r(:,n0+1:n0+n1));
r=r(:,1:n0)/L;
dr=dr(:,1:n0)/L;
d2r=d2r(:,1:n0)/L;
d3r=d3r(:,1:n0)/L;

if 0
  plot(r(1,:),r(2,:),'k');
end