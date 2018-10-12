 function [xptsout,yptsout] = invTransformCrv_only1pt(xpts,ypts,m,S,th,t0,t1)

    p0 = [xpts'; ypts'];
    if m
        p0 = [-xpts'; ypts'];
    end
    %scale
    
    S = 1/S;
    p0 = S*p0;

    R = [cos(-th) -sin(-th);sin(-th) cos(-th)];
    
    p0 = R*p0;

    
    %translate
    t = [t0;t1];
    p0 = p0 +t;

    
    %transform Bezier curve
    xptsout = p0(1,:)';
    yptsout = p0(2,:)';

 end