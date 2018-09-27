function [colin] = colinear4(p0,p1,p2,p3,tol)
%Check if 4 points are co-linear
if nargin < 5
    tol = 1;
end

val0 = abs(colinear3(p0,p1,p2));
val1 = abs(colinear3(p0,p1,p3));
val2 = abs(colinear3(p0,p2,p3));
val3 = abs(colinear3(p1,p2,p3));

val = val0+val1+val2+val3;
colin = false;
if val < tol
    colin = true;
end
end