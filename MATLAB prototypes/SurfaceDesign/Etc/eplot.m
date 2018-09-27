function [X1,X2]=eplot(p,n)
if nargin<2
    n=100;
end
%computes the ielastic curve with  parameters 
% p= [k, s0, L, S, phi, x0, y0]
%s0 is the initial point on a stadard representation of the elastic curve
%with parameter k.  L is the length os the segment is the standard rep.
%S is a scaling factor. phi is a rotation in the plane. (x0,y0) is a
%translation in the plane.

k=p(1); s0=p(2); L=p(3); S=p(4); phi=p(5); 
x0=p(6); y0=p(7);
t=s0:L/n:(s0+L);

options=odeset('RelTol',0.001,'AbsTol',0.001);
v0=0;
if k<=1
    
    if s0 ~= 0
      t1=[0 s0];
    [~,sol1]=ode45(@(x,y)ellipj(x,k^2)^2,t1,0,options); 
    v0=sol1(end);
    end
  [~,sol] = ode45(@(x,y)ellipj(x,k^2)^2,t,v0,options); 
  [~,cn,~] =  ellipj(t,k^2);
   ex=t-(2*k^2*sol)';
else
     if s0 ~= 0
     t1=[0 s0];
    [~,sol1]=ode45(@(x,y)ellipj(x,(1/k)^2)^2,k*t1,0,options); 
         v0=sol1(end);
     end
  [~,sol] = ode45(@(x,y)ellipj(x,(1/k)^2)^2,k*t,v0,options); 
  [~,~,dn1] =  ellipj(k*t,(1/k)^2); 
  cn=dn1;
   ex=t-(2*(1/k)*sol)';  
end
ey=2*k*(1-cn);
X1=x0 + S*(cos(phi)*ex-sin(phi)*ey);
X2=y0 +S*(sin(phi)*ex+cos(phi)*ey);

end