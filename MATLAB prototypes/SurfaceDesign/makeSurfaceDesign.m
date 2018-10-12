
function [curves,obj,lambdaresi_curves] = makeSurfaceDesign(opt)

%The database of lambda residuals:
s = load('databaseLambdares.mat','G11','G12','angerrors','angerrorsG12');
n = fieldnames(s);
for k=1:length(n)
    eval(sprintf('global %s; %s=s.%s;',n{k},n{k},n{k}));
end

%Initialization
[curves,kn] = initCrvs(opt);

%Design part

%Create the nurbs surface
surffig = figure('Position',[100 300 450 450]);
set(gcf,'renderer','opengl');
obj = PlotNurbsDesign(curves,opt,kn);
light

%make the interpol curves on the surface
CurveHandles = PlotCurvesOnDesign(curves,opt);

lambdaresi_curves = []; %lambda residuals of design curves

%Select curve
finished = 0;
while finished == 0
    disp('Choose curve to modify surface.')
    [obj,CurveHandles,curves,lambdaresi_curves] = CurveToModify(opt,obj,CurveHandles,curves,kn,surffig);
    if numel(curves) == 0
        finished = 1;
        curves = curves_bef;
    else
        curves_bef = curves;
    end
end

end