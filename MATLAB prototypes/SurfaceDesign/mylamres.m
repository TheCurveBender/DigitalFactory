function [lambdares,lambda1,lambda2,detA]= mylamres(p0,p1,p2,p3,ng)
if nargin < 5
    ng = 550;
end

cx = [p0(1),    p1(1),   p2(1) ,  p3(1)]; %x coordinates of the 4 control points, 
cy = [p0(2)  ,  p1(2) ,  p2(2) ,  p3(2)]; %y
kn = [0,0,0,0,1,1,1,1];

deg = 3;
[gt,gw] = legpts(ng,[0 1]);
gt = gt';

xsp   = spmak(kn,cx);
ysp   = spmak(kn,cy);
dxsp  = fnder(xsp);
dysp  = fnder(ysp);
d2xsp = fnder(dxsp);
d2ysp = fnder(dysp);
%   Evaluation at the Gauss points
xg    = fnval(xsp,gt);
yg    = fnval(ysp,gt);
dxg   = fnval(dxsp,gt);
dyg   = fnval(dysp,gt);
d2xg  = fnval(d2xsp,gt);
d2yg  = fnval(d2ysp,gt);

%Speed, tangent and curvature

speed = sqrt(dxg.^2 + dyg.^2);
kappa = (dxg.*d2yg - dyg.*d2xg)./speed.^3;

%   Gauss weights for computing the arc length
gsw   = gw.*speed; 

L = sum(gsw);

A = zeros(3,3);

A(1,1) = sum(yg.^2.*gsw);
A(1,2) = -sum(xg.*yg.*gsw);
A(1,3) = -sum(yg.*gsw);
A(2,1) = A(1,2);
A(2,2) = sum(xg.^2.*gsw);
A(2,3) = sum(xg.*gsw);
A(3,1) = A(1,3);
A(3,2) = A(2,3);
A(3,3) = L;
detA = det(A);

b = zeros(3,1);
b(1) = -sum(kappa.*yg.*gsw);
b(2) = sum(kappa.*xg.*gsw);
b(3) = sum(kappa.*gsw);

lambda = A\b;
lambda1 = lambda(1);
lambda2 = lambda(2);
alpha   = lambda(3);

lambdares = sqrt(...
           (kappa + lambda1*yg - lambda2*xg - alpha).^2*gsw'/...
           (kappa.^2*gsw'));
end
 