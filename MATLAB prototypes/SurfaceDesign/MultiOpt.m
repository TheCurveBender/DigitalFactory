classdef MultiOpt
   properties
    nBlocks
    bW
    bH
    bL
    nICurves
    yvals
    Manual
    Plot
    nCurves 
    max_curvature_con
    ExtLength
    nTotal
    PlotExt
    Prev
    prevtol
   end
   methods
   function obj = MultiOpt(varargin)
       if mod(nargin,2) ~= 0
           error('Number of arguments must be even.')
       end
       %default values
       obj.nBlocks = 1;
       obj.bW = 5.5; obj.bH =5.; obj.bL = 5.5; obj.nICurves = 3; 
       obj.yvals = linspace(0,5.5,3);
       obj.Manual = false;
       obj.Plot = true; obj.nCurves = 8; obj.max_curvature_con =1;
       obj.ExtLength = 8;
       obj.nTotal = 150;
       obj.PlotExt = true;
       obj.Prev = false;
       
       %obj.prevtol = inf;

       %set values
       nVarargs = length(varargin);
       for k = 1:(nVarargs/2)
           if isnumeric(varargin{2*k})
               val = num2str(varargin{2*k});
           else
               val = varargin{2*k};
           end
          eval(strcat('obj.',varargin{2*k-1},'=',val));
       end
       
   end
end
end