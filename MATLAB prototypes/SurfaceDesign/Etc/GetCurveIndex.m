function [index] = GetCurveIndex(opt,obj,CurveHandles,curves,kn)


old_colors = zeros(3,opt.nICurves*opt.nBlocks);
for j = 1:opt.nICurves*opt.nBlocks
    old_colors(:,j) = CurveHandles(j).Color;
end

disp('Click anywhere on the surface, then hit return to select the curve.')
pause
[pout]  = select3d(obj);

if numel(pout) == 0
    index = inf;
else
    yvals= opt.yvals; %linspace(0,opt.bL,opt.nICurves);
    [~, index] = min(abs(yvals-pout(2)));
    index_long = 1+opt.nBlocks*(index-1);
    for j = 1:opt.nICurves*opt.nBlocks
        if index_long <=  j && j <index_long+opt.nBlocks
            CurveHandles(j).Color = [1 1 0];
        else
            CurveHandles(j).Color = old_colors(:,j);
        end
    end
end

end
