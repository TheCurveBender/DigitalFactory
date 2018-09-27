function[par,lambdares,errtan,p1,p2,L,ela] = FirstGuessD(p1,p2,p0,p3,isplot,n_curves)
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

%yc = yc+[0 0 0 0.0001]
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
[~,~,errtan]=matchends(p0,p1,p2,p3,ela,par,isplot,n_curves);  
lambdares
end