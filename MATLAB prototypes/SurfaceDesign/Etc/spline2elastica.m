function [ela,L] =spline2elastica(kn,cx,cy)
%
% Input:
% kn: knots of spline
% cx: x-coordinates of control points
% cy: y-coordinates of control points
%
% Output:
% ela: struct with parameters for elastica
% err: residuals
%
deg = numel(kn)-numel(cx)-1;
kn  = kn - kn(deg+1); % start at 0
kn  = kn/kn(end-deg); % end at 1

% The splines and their derivatives
xsp   = spmak(kn,cx);
ysp   = spmak(kn,cy);
dxsp  = fnder(xsp);
dysp  = fnder(ysp);
d2xsp = fnder(dxsp);
d2ysp = fnder(dysp);

% use Gaussian quadrature on each knot interval
np      = 6*deg-3;    % Number of Gauss points
[gt,gw] = get_gauss_points_and_weights(kn(deg+1:end-3),np);
% The values in the Gauss points
xg    = fnval(xsp,gt);
yg    = fnval(ysp,gt);
dxg   = fnval(dxsp,gt);
dyg   = fnval(dysp,gt);
d2xg  = fnval(d2xsp,gt);
d2yg  = fnval(d2ysp,gt);
% speed, tangent, and curvature
speed  = sqrt(dxg.^2 + dyg.^2);
tx     = dxg./speed; 
ty     = dyg./speed;
det12  = (dxg.*d2yg - dyg.*d2xg);
kappa  = det12./speed.^3;
% weights for integration with respect to arc-length
gsw   = gw.*speed;  % ds = speed*dt
% find s(t) in Gauss points
gs    = find_arclength(dxsp,dysp,np);
[ela,L] = find_elastica(xg,yg,tx,ty,kappa,gs,gsw);
end