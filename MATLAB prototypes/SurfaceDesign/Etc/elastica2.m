function [x,y,dx,dy,kappa,s,theta,dkappa]=elastica2(k,s0,l,scale,phi,tt)
% Calculates points etc, on elastica without inflections (high energy)
% k: modulus
% s0: start parameter
% l: length/scaling
% scale: scaling
% phi: angle
% tt: parameter values in [0,1]
%
% x,y: curve
% dx,dy: tangent
% kappa: curvature
% s: arclength starting at 0
s = s0+l*tt;
K = 4*ellipke(k^2);

[sn,cn,dn] = ellipj(s/k,k^2);


am = atan2(sn,cn);    % -pi< am <= pi

I  = find(am<0);
am(I) = am(I) + 2*pi; % now: 0<= am < 2*pi
n  = floor(s/K);
am = am + 2*pi*n;

E = ellipticE(am,k^2);

x0 = scale*(2/k*E + (1-2/k^2)*s);
y0 = scale*2/k*(1-dn);
s = scale*l*tt;
x = cos(phi)*x0 - sin(phi)*y0;
y = sin(phi)*x0 + cos(phi)*y0;

% normalized derivative
dx0 = 1 - 2*sn.^2; 
dy0 = 2*sn.*cn;
dx = cos(phi)*dx0 - sin(phi)*dy0;
dy = sin(phi)*dx0 + cos(phi)*dy0;
theta = atan2(dy,dx);
theta = unwrap(theta);%+2*n*pi;
kappa = 2/k*dn/scale;
if nargout > 7
    dkappa = -2*sn.*cn/scale^2;
end
if l<0
    dx    = -dx;
    dy    = -dy;
    theta = -theta;
    kappa = -kappa;
    s     = -s;
end
end