function [val] = colinear3(p0,p1,p2)
%Check if 3 points are co-linear
val = det([1 p0(1) p0(2);1 p1(1) p1(2); 1 p2(1) p2(2)]);
end