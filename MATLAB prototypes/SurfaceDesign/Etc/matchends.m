function [XE,YE, errtan,p1e,p2e]=matchends(p0,p1,p2,p3,ela,par,isplot,n_curves)
global curvemode
if nargin<7
    isplot=0;
end
residues=ela.residue;  lambdares=residues(1);
[xe,ye]=  eplot(par,n_curves-1);
p0el=[xe(1),ye(1)];
p3el=[xe(end),ye(end)];
LE=p3el-p0el;
LB=p3-p0;
S=norm(LB)/norm(LE);
costh=dot(LE,LB)/(norm(LE)*norm(LB));
sinth=det([LE;LB])/(norm(LE)*norm(LB));
Rot=[costh, -sinth; sinth, costh];
Tx=ela.Tx;
Ty=ela.Ty;
tangentsE=Rot*[Tx(1), Tx(end); Ty(1), Ty(end)];
V1=tangentsE(1:2,1:1)';
V2=tangentsE(1:2,2:2)';
p1e=p0+norm(p1-p0)*V1;
p2e=p3-norm(p3-p2)*V2;
errth1=acos(dot( V1 ,p1-p0)/(norm(V1)*norm(p1-p0)))*2/pi;
errth2=acos(dot(V2,  p3-p2)/(norm(V2)*norm(p3-p2)))*2/pi;
errtan=0.5*(errth1+errth2);
matched1=S*Rot*[xe; ye];
matched=p0'-matched1(:,1:1)+matched1;

XE=matched(1:1,:);
YE=matched(2:2,:);
if isplot 
    colorL=[min(10*errth1,1), 0.8*(1-min(10*errth1,1)), 0];
    colorR=[min(10*errth2,1), 0.8*(1-min(10*errth2,1)), 0];
    cs=min(4*lambdares,1);
    curvecolor=cs*[0.5,0.5,0.5]+(1-cs)*[0 0.9 0];
    polyplot(p0,p1,p2,p3,colorL,colorR);
    plot(XE,YE,'--','color',[0.75 0.75 0.75],'linewidth',1.7); axis equal
    t=linspace(0,1,51);
    curve = kron((1-t).^3,p0') + kron(3*(1-t).^2.*t,p1') + kron(3*(1-t).*t.^2,p2') + kron(t.^3,p3');
    x=curve(1,:); y=curve(2,:);
end
end