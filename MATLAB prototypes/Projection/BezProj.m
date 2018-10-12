function returndata = BezProj(x,thresh,n,isplot)
%This is a demo of the feed-back based projection algorithm described in the article
%"Bezier curves that are close to elastica" by David Brander, Jakob Andreas Bærentzen, Ann-Sofie Fisker, Jens Gravesen
%The input is an arbitrary Bezier curve, and the output is a Bezier curve
%with the same endpoints and end-tangents, but different length.
%The desired closeness to an elastic curve is set by "thresh", a number
%between 0 and 1.  Default is 0.2, and anything less than 0.3 should be
%visually very close to an elastic curve.  

%Note: the algorithm is designed to guarantee a good output for Bezier
%curves that satisfy certain angle constraints on the end-tangents
%(described in the paper).  For curves that deviate too far from these, the
%result may not be close to an elastic curve.  This is measured by the
%first return data,  'lr', the lambda residual.   Roughly speaking, 
%good curves have lr<0.3 and bad curves have lr>0.5

%Author: David Brander,   email:  dbra@dtu.dk

%Input: 
%                    x= [x0,y0, x1,y1, x2,y2, x3,y3] 
%for the Bezier curve with control polygon p0=[x0,y0], p1=[x1,y1], p2=[x2,y2], p3=[x3 y3]

% Other input are optional: 
% thresh is between 0 and 1, the threshold value for the lambda residual which terminates the search
%n should usually be 1 or 2.
%isplot=1 to plot the solution in a figure

%Output: returns a row-vector of length 9
%             returndata=[lr, xproj]    
%where lr is the lambda residual of the projected Bezier curve,
%xproj=  [x0,y0, x1,y1, x2,y2, x3,y3]  is the 4 control polygon vertices of
%the projected Bezier curve

%Example:
%      BezProj([0,0,-1,1,0.5,1,1,0],0.4);

if nargin<4
    isplot=1; %plot the solution
end
if nargin<3
   n=2;  % number of passes
end
if nargin<2
   thresh=0.2;  %target value for lambda residual
end
minLL=0.27;  %minimal edge-length for inflectional curves
speed=0;  %speed=0  -> 21 tvalues, more accurate
                %speed=1 -> 11 tvalues, faster, usually the same result. 
p0=[x(1);x(2)]; p1=[x(3);x(4)]; 
p2=[x(5);x(6)]; p3=[x(7);x(8)];
 if isplot
 figure()
 bezplot4(p0,p1,p2,p3,'--',[0.5 0.5 0.5]);  %
 hold on
xlima=xlim; ylima=ylim;
end

% put the polygon into standard form:
 tran=p0;
    p3b=p3-tran;
    theta=atan2(p3b(2),p3b(1));
    A=(1/norm(p3b))*[cos(theta), sin(theta); -sin(theta), cos(theta)];
    p1=A*(p1-tran);
    p2=A*(p2-tran);
 %%%%%%%%%
    x1=p1(1); y1=p1(2);
    x2=p2(1); y2=p2(2);
 
lr0=lamres64(x1,y1,x2,y2);
  [phi1,phi2] = calcphi(x1,y1,x2,y2);
  inf=~(sign(phi1-pi)==sign(phi2-pi)); %detects whether inflectional or not
  data=[inf,lr0,0,x1,y1,x2,y2];
  iter1=0; iter2=0; 
  
if lr0>thresh||min(x1.^2+y1.^2, (x2-1).^2+y2.^2)<minLL^2 %if lr0 less than thresh we get our original data
  [x1,y1,x2,y2,inf]=initializecurve(x1,y1,x2,y2,minLL);%remove crossings, etc.
end

  if inf
    [data,~]= projInf([x1,y1,x2,y2],thresh,n,minLL,speed);
    inf=data(1);
    iter1=data(3);
  end
  if ~inf
     [data,~,~,~] = projNI(data(4:7),thresh,n,minLL,speed);
     iter2=data(3);
  end 
   p1=A^(-1)*data(4:5)'+tran;
   p2=A^(-1)*data(6:7)'+tran;
   returndata=[data(2:2), p0',p1', p2',p3'];
  if isplot
     bezplot4(p0,p1,p2,p3,'-',[0 0.5 0]);  %
     xlimb=xlim; ylimb=ylim;
     Xlim=[min(xlima(1),xlimb(1)), max(xlima(2),xlimb(2))];
     Ylim=[min(ylima(1),ylimb(1)), max(ylima(2),ylimb(2))];
     lrplotheight=0.1*(Ylim(2)-Ylim(1));
     Ylim=[Ylim(1), Ylim(2)+1.1*lrplotheight];
    xlim(Xlim); ylim(Ylim);
    title(['lres input: ',num2str(lr0,2),'. lres out: ',num2str(data(2),2),...
        '. Itns:',num2str(iter1),',  ',num2str(iter2)],...
      'fontsize',10,'interpreter','latex');
  end
end

function   [x1,y1,x2,y2,inf]=initializecurve(x1,y1,x2,y2,minLL)
      beta1=atan2(y1,x1);
 beta2=atan2(y2,x2-1);  %the angles of the legs to the x-axis
 angdiff=sign(beta2).*(beta1-beta2);
     LMin=max(0.4*(1+6*angdiff/pi),minLL); %0.27
      LMax=max(1.2*(1+5.0*angdiff/pi), 0.58);
         L0=sqrt(x1.^2+y1.^2);
         L2=sqrt((x2-1).^2+y2.^2);
         u1=x1/L0; v1=y1/L0; 
         u2=(x2-1)/L2; v2=y2/L2;
        L0=max(L0, LMin); 
         L2=max(L2,  LMin);
         %inflectional type is determined by setting the leg-length to minLL
         L0test=min(L0,LMax);  L2test=min(L2,LMax);
        [phi1,phi2] = calcphi(L0test*u1,L0test*v1,1+L2test*u2,L2test*v2);
        inf=~(sign(phi1-pi)==sign(phi2-pi));
          x1=L0*u1;
          y1=L0*v1;
          x2=1+L2*u2; 
          y2=L2*v2; 
   if (x1-1)^2+y1^2<0.1  %avoid coincidence of control points
      x1=0.75*x1; y1=0.75*y1;
  end
  if x2^2+y2^2<0.1
      x2=1+0.75*(x2-1); y2=0.75*y2;
  end
  if y1*(y1*x2-x1*y2)<=0&&-y2*(y2*x1-y2-y1*x2+y1)<=0 %remove intersection
    X=x1*y2/(x1*y2+y1*(1-x2));
    Y=y1*y2/(x1*y2+y1*(1-x2));
    x1=0.9*X; y1=0.9*Y;
    x2=1+0.9*(X-1);  y2=0.9*Y;
  end
 end


        
function [phi1,phi2] = calcphi(x1,y1,x2,y2)
u1=[x1';y1';zeros(1,length(x1))];
u2=[x1'-x2'; y1'-y2';zeros(1,length(x1))];
x2s=x2-1;
u3=[x2s'; y2';zeros(1,length(y2))];
th1=acos(   dot(u1,u2)./sqrt(   dot(u1,u1).*dot(u2,u2) ))';
th2=acos(dot(u3,-u2)./sqrt(dot(u2,u2).*dot(u3,u3)))';
cross1=cross(u1,u2);
s1=sign(cross1(3:3,:))';
cross2=cross(u3,u2);
s2=sign(cross2(3:3,:))';
phi1=s1.*th1+(1-s1)*pi;
phi2=s2.*th2+(1-s2)*pi;
end


function bezplot4(p0,p1,p2,p3,style,col)
%Plot's a cubic Bezier curve
%iniput bez=[x0 x1 x2 x3; y0 y1 y2 y3]
%ouput [x,y] are the xy values at n poins along the curve
if nargin<6
    col=[0 0 0];
end
if nargin<5
    style='-';
end
n=51;

t=linspace(0,1,n);
curve = kron((1-t).^3,p0) + kron(3*(1-t).^2.*t,p1) + kron(3*(1-t).*t.^2,p2) + kron(t.^3,p3);
l1=t.*p0+(1-t).*p1;
l2=t.*p1+(1-t).*p2;
l3=t.*p2+(1-t).*p3;
lw=2;
x=curve(1,:);
y=curve(2,:);
plot(x,y,'lineWidth',1.5,'LineStyle',style,'color',col)
hold on
plot(l1(1,:),l1(2,:),'LineStyle',style,'color',[0.4 0.4 0.4],'lineWidth',0.4*lw);
plot(l2(1,:),l2(2,:),'LineStyle',style,'color',[0.4 0.4 0.4],'lineWidth',0.4*lw);
plot(l3(1,:),l3(2,:),'LineStyle',style,'color',[0.4 0.4 0.4],'lineWidth',0.4*lw);
plot(p1(1),p1(2), 'or')
plot(p2(1),p2(2), 'og')
plot(p0(1),p0(2), 'ok','MarkerFaceColor',[0 0 0])
plot(p3(1),p3(2), 'ok', 'MarkerFaceColor', [0 0 0])
axis equal

end




%%%%%%%%%%%% Projection function for non-inflectional curves%%
function [data,lrlist, NIL1,NIL2]  = projNI(x,thresh,n,minLL,speed)
%projection for non-inflectional curves
if nargin<4
    minLL=0.25;
    speed=1;
end
m=21;
if speed==1
    m=11;
end
if nargin<3
    n=2;   %number of passes
end
if nargin<2
   thresh=0.2;
end
x1=x(1); y1=x(2); x2=x(3); y2=x(4);

 [~,lr,inf]=checkNI(x1,y1,x2,y2,1,1,thresh);
lrlist=lr;
iter=0;
 [NIL1,NIL2]=getoptimal(x1,y1,x2,y2);
if lr>thresh
    [x1,y1,x2,y2,u1,u2,v1,v2,L0,L2]=adjustlegs(x1,y1,x2,y2,minLL);
    [OK,lr,inf]=checkNI(x1,y1,x2,y2,NIL1,NIL2,thresh);
    lrlist=[lrlist, lr];
    iter=1;
end

for i=1:n
       if lr>thresh
           iter=1+i;
           t=linspace(0,1,m);  lrs=[lr,ones(1,m-1)]; infs=[inf,zeros(1,m-1)];
           for j=2:m
                if ~OK
                    x1=((1-t(j))*L0+t(j)*NIL1)*u1;
                    y1=((1-t(j))*L0+t(j)*NIL1)*v1;
                    x2=1+((1-t(j))*L2+t(j)*NIL2)*u2;
                    y2=((1-t(j))*L2+t(j)*NIL2)*v2;
                    [OK,lr,inf]=checkNI(x1,y1,x2,y2,NIL1,NIL2,thresh);
                    lrs(j)=lr;  infs(j)=inf;
                else
                     j=j-1;
                     break
                end
           end
           % we choose the result with lowest value of lr
           [lr,ind]=min(lrs(1:j));
           lrlist=[lrlist, lrs(1:ind)];
           inf=infs(ind);
           L0=(1-t(ind))*L0+t(ind)*NIL1;
           L2=(1-t(ind))*L2+t(ind)*NIL2;
           x1=L0*u1;    y1=L0*v1;  x2=1+L2*u2;   y2=L2*v2;  
       end
end
       
data=[inf,lr,iter,x1,y1,x2,y2];

end

function  [NIL1,NIL2]=getoptimal(x1,y1,x2,y2)
      beta1=atan2(y1,x1);
     beta2=atan2(y2,x2-1); %the angles of the legs to the x-axis
       psi1=atan2(y1,-x1);  %the clockwise angle from -ve y axis to L0
       psi2= atan2(y2,-x2+1); %the clockwise angle from x axis to L2
       NIL1=max(0.12,f(sign(y1)*beta1, sign(y1)*beta2));
       NIL2=max(0.12, f(sign(y2)*psi2, sign(y2)*psi1));
end

 function   [x1,y1,x2,y2,u1,u2,v1,v2,L0,L2]=adjustlegs(x1,y1,x2,y2,minLL)
      LMin=minLL;
         L0=sqrt(x1.^2+y1.^2);
         L2=sqrt((x2-1).^2+y2.^2);
         u1=x1/L0; v1=y1/L0; 
         u2=(x2-1)/L2; v2=y2/L2;
        L0=max(L0-0.01, LMin-0.001); %to stop inflections
          L2=max(L2-0.01,  LMin-0.001);
          x1=L0*u1;
          y1=L0*v1;
          x2=1+L2*u2; 
          y2=L2*v2; 
 end
        
 
    function [OK,lr,inf]=checkNI(x1,y1,x2,y2,NIL1,NIL2,thresh)
    lr=lamres64(x1,y1,x2,y2);
      L0=sqrt(x1.^2+y1.^2);
      L2=sqrt((x2-1).^2+y2.^2);
      mycross=y1.*(y1.*x2-x1.*y2)<=0&-y2.*(y2.*x1-y2-y1.*x2+y1)<=0;
      [phi1,phi2] = calcphi(x1,y1,x2,y2);
      inf=~(sign(phi1-pi)==sign(phi2-pi));
       if nargout > 2
           OK=((~mycross)& lr<thresh)|...
           (abs(L2-NIL2)<=0.001&abs(L0-NIL1)<=0.001);
       end
    end
    

function R=f(x,y)
p0 =   1.475e-05;
       p1 =       10.39 ;
       p2 =      -10.48 ;
       p3 =      0.4574 ;
       p4 =       1.711 ;
       p5 =      -2.535 ;
       p6 =       2.772 ;
       p7 =    -0.08504 ;
       p8 =     -0.9109 ;
       p9 =     -0.2957 ;
       p10 =     -0.6606 ;
 R=p0*exp(x*p1+p2*y)+x*p3*exp(x*p4+p5*y)+...
     y*p6*exp(x*p7+p8*y)+x*y*p9*exp(x*p10);

end


%%%%%%%%%%%% Projection function for inflectional curves%%
function [data,lrlist]  = projInf(x,thresh,n,minLL,speed)
%projection for inflectional curves
if nargin<4
    minLL=0.27;
    speed=1;
end
m=21;
if speed==1
    m=11;
end
if nargin<3
   n=2; %number of passes
end
if nargin<2
   thresh=0.2;
end
x1=x(1); y1=x(2); x2=x(3); y2=x(4);

 beta1=atan2(y1,x1);
 beta2=atan2(y2,x2-1);  %the angles of the legs to the x-axis
 angdiff=sign(beta2).*(beta1-beta2);
 LMin=max(0.4*(1+6*angdiff/pi),minLL); 
 LMax=max(1.2*(1+5.0*angdiff/pi), 0.58);
 [~,lr,inf]=checkInf(x1,y1,x2,y2,LMin,thresh);
lrlist=lr; iter=0;
if lr>thresh
        L0=sqrt(x1.^2+y1.^2);
         L2=sqrt((x2-1).^2+y2.^2);
         u1=x1/L0; v1=y1/L0; 
         u2=(x2-1)/L2; v2=y2/L2;
         L0=max(L0,minLL);
         L2=max(L2,minLL);
         L0=min(L0, LMax);
         L2=min(L2, LMax);  
          x1=L0*u1;
          y1=L0*v1;
          x2=1+L2*u2; 
          y2=L2*v2; %adjust curve 
          [OK,lr,inf]=checkInf(x1,y1,x2,y2,LMin,thresh);
          lrlist=[lrlist, lr];
         iter=1;
end

for i=1:n
   if lr>thresh %check if adjusted curve need to be projected
       iter=1+i;
       t=linspace(0,1,m);  lrs=[lr,ones(1,m-1)]; infs=[inf,zeros(1,m-1)];
       for j=2:m
           if ~OK %continue until ningj or proj z or e tol
             x1=((1-t(j))*L0+t(j)*LMin)*u1;
             y1=((1-t(j))*L0+t(j)*LMin)*v1;
             x2=1+((1-t(j))*L2+t(j)*LMin)*u2; 
             y2=((1-t(j))*L2+t(j)*LMin)*v2;  
            [OK,lr,inf]=checkInf(x1,y1,x2,y2,LMin,thresh);
            lrs(j)=lr;  infs(j)=inf;
           else
               j=j-1;
               break
           end
       end
        % we choose the result with lowest value of lr
       [lr,ind]=min(lrs(1:j));
        lrlist=[lrlist, lrs(1:j)];
        inf=infs(ind);  
        L0=(1-t(ind))*L0+t(ind)*LMin;
        L2=(1-t(ind))*L2+t(ind)*LMin;
        x1=L0*u1;    y1=L0*v1;  x2=1+L2*u2;   y2=L2*v2;  
   end
end
data=[inf,lr,iter,x1,y1,x2,y2];
end

    function [OK,lr,inf]=checkInf(x1,y1,x2,y2,Min,thresh)
    lr=lamres64(x1,y1,x2,y2);
      L0=sqrt(x1.^2+y1.^2);
      L2=sqrt((x2-1).^2+y2.^2);
      mycross=y1.*(y1.*x2-x1.*y2)<=0&-y2.*(y2.*x1-y2-y1.*x2+y1)<=0;
      [phi1,phi2] = calcphi(x1,y1,x2,y2);
      inf=~(sign(phi1-pi)==sign(phi2-pi));
       if nargout > 2
           OK=((~mycross)& lr<thresh)|...
           (L0<=Min+0.001&L2<=Min+0.001)| ...
           ~inf;
       end
    end
    