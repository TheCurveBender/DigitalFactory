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

        %Compute lambda residual of curve

        lambdares = mylamres(p0,p1,p2,p3,19);
        cx=[p0(1),p1(1),p2(1),p3(1)];
        cy=[p0(2),p1(2),p2(2),p3(2)];
        spl = spmak([0,0,0,0,1,1,1,1],[cx;cy]);
        xy = fnval(spl,linspace(0,1,100));
        lambdaresi_curves = [lambdaresi_curves  lambdares];
       
         if lambdares>0.35
          colorE=[1 0 0];
         else
          colorE=[0 1 0];
        end
                
        XVALS = xy(1,:); YVALS = repmat(yvals(j),1,numel(xy)/2); ZVALS = xy(2,:);
        hold on
        FIG = plot3(XVALS,YVALS,ZVALS,'Color',colorE);
        listOfHandles = [listOfHandles FIG];
    end
end

end