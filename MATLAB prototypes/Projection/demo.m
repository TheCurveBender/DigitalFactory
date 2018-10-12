% Demo of the function Matlab function 'BezProj.m', which demonstrates the
% projection algorithm in the article "Bezier curves that are close to
% elastica".  

%The function BezProj is designed to take as input Bezier curves that satisfy certain angle
%constraints on the end-tangents (described in the article).  Therefore, this 
%script creates a set of Bezier curves each given by 4 control points
%randomly placed in the unit disc.  It then selects those that satisfy the
%angle constraints.  Finally it applys the function 'BezProj' to 10 of
%these curves, with two different threshold levels for the lambda residual
%(0.3 and 0.1).

%The angle constraints can be changed here:
%%%%%%%%%% parameters that determine the angle constraints
ab2=pi/3; %set to 0 to get completely general curves
ab=0.4*pi;  % set to pi to get completely general curves
%%%%%%%%%%%%%

%% create some random points p0 p1, p2, p3 in the unit disc
X=cell(4,1); Y=cell(4,1); mm=zeros(4,1);
for i=1:4
 x=2*rand(6,1)-1;
 y=2*rand(6,1)-1;
 I=find(x.^2+y.^2<=1);
 X{i} = x(I);  Y{i}=y(I);
 mm(i)=numel(I);
end
m=min(mm);
x=zeros(m,4); y=zeros(m,4);
for i=1:4
    Xt=X{i}; Yt=Y{i};
    X{i}=Xt(1:m); Y{i}=Yt(1:m);
end
 
a1=[kron(ones(m,1),[X{1},Y{1}]), kron([X{2},Y{2}],ones(m,1))];
a2=[kron(ones(m,1),a1), kron([X{3},Y{3}],ones(m^2,1))];
pp=[kron(ones(m,1),a2), kron([X{4},Y{4}],ones(m^3,1))];
J=randperm(m^4)';
p=pp(J,:);

%find the points that satisfy the angle constraints
 x0=p(:,1:1); y0=p(:,2:2); x1=p(:,3:3); y1=p(:,4:4);
 x2=p(:,5:5); y2=p(:,6:6);x3=p(:,7:7); y3=p(:,8:8);
 
%coordinate change that takes (0,0) -> (x0,y0)
x1=x1-x0; y1=y1-y0;
x2=x2-x0; y2=y2-y0;
x3=x3-x0; y3=y3-y0;
% coordinate change that takes (1,0) -> (x3,y3)
r=sqrt(x3.^2+y3.^2); 
th = atan2(y3,x3);
nx1=(1./r).*(cos(th).*x1 + sin(th).*y1);
ny1=(1./r).*(-sin(th).*x1+cos(th).*y1);
nx2=(1./r).*(cos(th).*x2 + sin(th).*y2);
ny2=(1./r).*(-sin(th).*x2+cos(th).*y2);

x1=nx1; y1=ny1; x2=nx2; y2=ny2;


beta1=mod(2*pi+atan2(y1,-x1),2*pi);  %the clockwise angle from -ve y axis to L0
beta2=mod(2*pi+atan2(y2,x2-1),2*pi);   %the counterclockwise angle from x axis to L2
ang1= beta1>ab2&beta1<2*pi-ab2& beta2>ab2&beta2<2*pi-ab2;
ang2= abs(beta1-beta2)<ab; %satisfy the ang 1-2 constraints!
J=find(ang1&ang2);

p=p(J,:);
for j=1:12
    projgendemo2(p(j:j,:))
end


function projgendemo2(x)
thresh=[0.2]; %0.2 threshold for lambda-residual
p00=[x(1);x(2)]; p10=[x(3);x(4)]; p20=[x(5);x(6)]; p30=[x(7);x(8)];

lr=zeros(1,2);
t=linspace(0,1,101);
figure()
hold on
incurve = kron((1-t).^3,p00) + kron(3*(1-t).^2.*t,p10) + kron(3*(1-t).*t.^2,p20) + kron(t.^3,p30);
lr0= lamres64(x(3),x(4),x(5),x(6),x(1),x(2),x(7),x(8));
plot(incurve(1,:),incurve(2,:),'color',[0.7 0.7 0.7],'lineWidth',1,'LineStyle','--')
cols=[0.4 0.4 0.8; 0 0.7 0;0.7 0 0];

for j=1:1
    col= cols(j:j,:);
    data=BezProj(x,thresh(j),2,0);
    lr(j)=data(1);
    x=data(2:9);
    p0=[x(1);x(2)]; p1=[x(3);x(4)]; p2=[x(5);x(6)]; p3=[x(7);x(8)];
    curve = kron((1-t).^3,p0) + kron(3*(1-t).^2.*t,p1) + kron(3*(1-t).*t.^2,p2) + kron(t.^3,p3);
    hold on
    plot(curve(1,:),curve(2,:),'color',col,'lineWidth',1)
    axis equal
end

legend({ [ 'Input Curve, $e_\lambda=$',num2str(lr0,2)],...
    [ 'Projection, $e_\lambda=$',num2str(lr(1),2)]...
    },'Interpreter','latex','FontSize',14);
polyplot(p00,p10,p20,p30,[0.7 0.7 0.7],'--')
polyplot(p0,p1,p2,p3,col,'-')
end
function polyplot(p0,p1,p2,p3,col,style)
t=linspace(0,1,2);
l1=t.*p0+(1-t).*p1;
l2=t.*p1+(1-t).*p2;
l3=t.*p2+(1-t).*p3;
lw=1;
plot(l1(1,:),l1(2,:),'color',col,'lineWidth',lw,'LineStyle',style);
plot(l2(1,:),l2(2,:),'color',col,'lineWidth',lw,'LineStyle',style);
plot(l3(1,:),l3(2,:),'color',col,'lineWidth',lw,'LineStyle',style);
plot(p1(1),p1(2), 'o','color',col)
plot(p2(1),p2(2), 'o','color',col)
plot(p0(1),p0(2), 'o','color',col,'MarkerFaceColor',col)
plot(p3(1),p3(2), 'o','color',col,'MarkerFaceColor',col)

end