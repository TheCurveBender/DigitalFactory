function[par,el,lambdares,errtan,L,p00,p01,t1,t2] = FirstGuessD_el(p0,p1,p2,p3,n_curves,isplot)
%For Bezier curves. Uses a first guess and then scales and translates
%to match the end points. Rotates the points p1 and p2 so that p1-p0 and
%p3-p2 match the tangents of the first guess elastic curves
global curvemode
if nargin<5
    isplot=0;
end
if nargin<4
    p3=[1,0];
end
if nargin<3
    p0=[0,0];
end
kn = [0,0,0,0,1,1,1,1];
xc=[p0(1),p1(1),p2(1),p3(1)];
yc=[p0(2),p1(2),p2(2),p3(2)];
[ela,L]= spline2elastica(kn,xc,yc);


if ela.inflexion
    k=ela.k;
else
    k=1/(ela.k);
end
s0=ela.s0;
l=ela.l;
scale=ela.scale;
phi=ela.phi;
x0=ela.x0;
y0=ela.y0;
par=[k,s0,l,scale,phi,x0,y0];
residues=ela.residue;
lambdares=residues(1);

[XE,YE,errtan,p1e,p2e]=matchends(p0,p1,p2,p3,ela,par,isplot,n_curves);
el=[XE;YE];
if curvemode==1  %change the control polygon to match the end-tangents of the elastica
    p1=p1e;
    p2=p2e;
end
[p00x,p00y,dx,dy,kappa,s,theta,dkappa]=elastica(ela ,0);    
p00 = [p00x; p00y];
t1 = [dx dy];
[p01x,p01y,dx,dy,kappa,s,theta,dkappa]=elastica(ela ,1);    
p01 = [p01x; p01y];
t2 = [dx dy];
end