 function [p0,p1,p2,p3,Database,th,t0,t1,S,m,angError] = TransformCrv(p0,p1,p2,p3,point)
 
 global G11 G12 angerrorsG12 angerrors
 
 if point == 1
     Database = G11;
     angError = angerrors;
 else
     Database = G12; 
     angError = angerrorsG12;
 end
    %transform Bezier curve
    p0 = p0'; p1 = p1'; p2 = p2'; p3 = p3';
    m = false;
    
    %translate
    p0 = p0 - p3;
    p1 = p1-p3;
    p2 = p2-p3;
    t0 = p3(1);
    t1 = p3(2);
    p3 = p3-p3;


    %rotate
    th = -(atan2(p2(2),p2(1))-pi/2); % -acos(p2(2)/norm(p2,2))
    %th = atan2(p2(2),p2(1))
    R = [cos(th) -sin(th);sin(th) cos(th)];

    p0 = R*p0;
    p1 = R*p1;
    p2 = R*p2;
    
    %scale
    S = 1/norm(p2,2);
    p0 = S*p0;
    p1 = S*p1;
    p2 = S*p2;
    p3 = S*p3;
    
    %mirror maybe?

    if ((p1(1)>0) && point == 1) || ((p0(1)>0) && point == 2) 
        p0 = [-p0(1); p0(2)];
        p1 = [-p1(1);p1(2)];
        m = true;
    end
    
    
    p0 = p0'; p1 = p1'; p2 = p2'; p3 = p3';
 end