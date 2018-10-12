function [obj] = PlotNurbsDesign(curves,opt,kn)

curves_input = cell(1,opt.nICurves);
yvals = opt.yvals;
cntrl = zeros(3,2*(opt.nBlocks+1));

for j = 1:opt.nICurves
     cntrl(1,:) =   curves(j,1:2:end);
     cntrl(2,:) = repmat(yvals(j),1,2*(opt.nBlocks+1));
     cntrl(3,:) =   curves(j,2:2:end);
     nurbs =   nrbmak(cntrl,kn);
     curves_input{j} = nurbs;
end

hold on
surface = nrbloft(curves_input);
obj = nrbplot_render(surface,[50 50],'EdgeColor','none','FaceColor',[0.8 0.8 0.8]);
obj.HitTest = 'off';
axis equal
%colormap gray
view(-36,22) 
end