function [listOfHandles,lambdaresi_curves] = PlotCurvesOnDesign(curves,opt)
yvals = opt.yvals;
listOfHandles = [];
lambdaresi_curves = [];
for j=1:opt.nICurves
    for i = 1:opt.nBlocks
        if opt.nBlocks == 1
            p0=curves(j:j,1:2);
            p1=curves(j:j,3:4);
            p2=curves(j:j,5:6);
            p3=curves(j:j,7:8);
        else
            if i == 1
                p0=curves(j:j,1:2);
                p1=curves(j:j,3:4);
                p2=curves(j:j,5:6);
                p3 = 0.5.*(p2+curves(j:j,7:8));
            elseif i == opt.nBlocks
                p0=0.5*(curves(j:j,end-5:end-4)+curves(j:j,end-7:end-6));
                p1=curves(j:j,end-5:end-4);
                p2=curves(j:j,end-3:end-2);
                p3 = curves(j:j,end-1:end);
            else
                p0 = p3;
                p1 = curves(j:j,7+(i-2)*4:8+(i-2)*4);
                p2 = curves(j:j,7+(i-2)*4+2:8+(i-2)*4+2);
                p3 = 0.5*(p2+curves(j:j,7+(i-2)*4+4:8+(i-2)*4+4));
            end
        end

        [~,elas2D,lambdares,errtan,~,~,~,~,~]  = FirstGuessD_el(p0,p1,p2,p3,100,0);
        lambdaresi_curves = [lambdaresi_curves  lambdares];
        %Plot IG given by elas2D
         if lambdares>0.35|| errtan>0.1
          colorE=[1 0 0];
        end
        if lambdares<0.2 && errtan<0.01
          colorE=[0 1 0];
        end
                
        if colinear4(p0,p1,p2,p3)
            colorE = [1 0 0];
        end
        XVALS = elas2D(1,:); YVALS = repmat(yvals(j),1,numel(elas2D)/2); ZVALS = elas2D(2,:);
        hold on
        FIG = plot3(XVALS,YVALS,ZVALS,'Color',colorE);
        listOfHandles = [listOfHandles FIG];
    end
end

end