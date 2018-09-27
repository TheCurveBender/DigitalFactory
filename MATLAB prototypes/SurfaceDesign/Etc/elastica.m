function [x,y,dx,dy,kappa,s,theta,dkappa]=elastica(e,tt)
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
if e.inflexion
    [x,y,dx,dy,kappa,s,theta,dkappa]=elastica1(e.k,e.s0,e.l,e.scale,e.phi,tt);
else
    [x,y,dx,dy,kappa,s,theta,dkappa]=elastica2(e.k,e.s0,e.l,e.scale,e.phi,tt);
end
x = x+e.x0;
y = y+e.y0;
end