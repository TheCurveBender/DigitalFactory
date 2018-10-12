function h=nrbplot_render(nrb,subd,varargin) 


% Function Name: 
%  
%   nrbplot - Plot a NURBS curve or surface. 
%  
% Calling Sequence: 
%  
%   h=nrbplot(nrb,subd,properties/values) 
%  
%  
% Parameters: 
%  
%  
%   nrb		: NURBS curve or surface, see nrbmak. It can be a structure
%             array or cell array of nurbs curves or surfaces
%  
%   subd	: array of number of subdivisions or cell array of evaluation
%             points along U and V directions. For ex: [nu, nv] or {vecu,vecv}
%   h       : Handle to the graphics object
%   poperties/values: Properties values pairs for surface plot
%
% Examples:
%
% The two statements below are equivalent and plot a nurbs surface with
% 
% nurbs=nrbtestsrf;
% h=nrbplot(nurbs,[50,50],'FaceColor','r','EdgeColor','k','FaceLighting'...
% ,'gouraud');
% light
% 
% h=nrbplot(nurbs,{linspace(0,1,50),linspace(0,1,50)},'FaceColor','r','Edge
% Color','k','FaceLighting','gouraud');
% light
% D.M. Spink 
% Copyright (c) 2000. 
 
 
nargs = nargin; 
if nargs < 2 
  error('Need a NURBS to plot and the number of subdivisions!'); 
elseif (nargs-2)/2~=floor(nargs-2)/2
    error('Need properties/value pairs');
end 
 
for i=1:length(nrb)
 
    if iscell(nrb)
        nurbs=nrb{i};
    else
        nurbs=nrb(i);
    end
% plot the curve or surface 
if iscell(nurbs.knots) 
  % plot a NURBS surface 
  if length(subd)~=2
      error('surface nurbs need evaluation points along U and V directions');
  end
  
  if iscell(subd)
     if max(subd{1})>1 || max(subd{2})>1 || min(subd{1})<0 || min(subd{2})<0
         error('values must be between 0 and 1');
     end
     evp=subd;  
  else    
     evp={linspace(0,1,subd(1)+1),linspace(0,1,subd(2)+1)};
  end
     p = nrbeval(nurbs,evp);
     h(i) = surf(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)),varargin{:});
else 
    
  % plot a NURBS curve 
  if isempty(subd)
      error('nurbs curves need evaluation points along U direction');
  end
  if iscell(subd)
     if max(subd{1})>1 ||  min(subd{1})<0 
         error('values must be between 0 and 1');
     end
     evp=subd{1};  
     
  else     
     evp=linspace(0,1,subd(1)+1);         
  end
  p = nrbeval(nurbs,evp); 
 
  if any(nurbs.coefs(3,:)) 
    % 3D curve 
    h(i)=plot3(p(1,:),p(2,:),p(3,:),varargin{:});  
    
  else 
    % 2D curve 
    h(i)=plot(p(1,:),p(2,:),varargin{:}); 
  end 
end 

end
axis equal; 
 

