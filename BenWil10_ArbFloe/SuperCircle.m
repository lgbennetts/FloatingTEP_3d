% function [dr,r,d2r,d3r]=SuperCircle(t,crk_prams)
%
% Description: Supercircle of form |x|^n + |y|^n=1
%              http://en.wikipedia.org/wiki/Squircle
%
% L Bennetts Nov 13 / Adelaide

function [dr,r,d2r,d3r]=SuperCircle(t,crk_prams)

t=reshape(t,1,length(t));

alf=1; %crk_prams{1};
n  =crk_prams{1};
k  =pi*alf;

if n<=2
 print('warning: Supercircle assumes n > 2')
 r=[0*t;0*t]; dr=r; d2r=r; d3r=r;
 return
end

c2n=sign(cos(k*t)).*(cos(k*t).^2).^(1/n);
s2n=sign(sin(k*t)).*(sin(k*t).^2).^(1/n);

x=c2n;
y=s2n;

r=[x;y];

x=-(2*k/n)*sin(k*t).*c2n.*(cos(k*t).^(-1));
y= (2*k/n)*cos(k*t).*s2n.*(sin(k*t).^(-1));

dr=[x;y];

x= (2*k/n)*k*(-1+2/n)*(sin(k*t).^2).*c2n.*(cos(k*t).^(-2)) ...
   -(2*k/n)*k*cos(k*t).*c2n.*(cos(k*t).^(-1));
y= (2*k/n)*k*(-1+2/n)*(cos(k*t).^2).*s2n.*(sin(k*t).^(-2)) ...
   -(2*k/n)*k*sin(k*t).*s2n.*(sin(k*t).^(-1));

d2r=[x;y];

d3r=[0*x;0*y];

return