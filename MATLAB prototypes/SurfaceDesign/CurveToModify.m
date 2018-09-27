function [obj,CurveHandles,curves,lambdaresi_curves] = CurveToModify(opt,obj,CurveHandles,curves,kn,surffig)
global active coefs csize R polys searchrad Mold boxcon Rx1 Rx2 Ry
global G11 G12 angerrorsG12 angerrors

lambdaresi_curves =[];
boxcon =  true;
searchrad = 0.1;
csize=3; %circle of control point size
R=opt.bH+3;
Ry = opt.bH;
if mod(opt.nBlocks,2) == 0
    Rx1 = opt.bW*(opt.nBlocks/2);
    Rx2 = Rx1;
else
    Rx1 = opt.bW*((opt.nBlocks-1)/2);
    Rx2 =opt.bW*((opt.nBlocks+1)/2);
end

figure(surffig)
[index] = GetCurveIndex(opt,obj,CurveHandles,curves,kn);

if index == inf
    curves = [];
else
    %% Open new window with the curve
    
    ModWindow = figure('Position',[600 300 450 450]);
    set(gcf,'renderer','opengl')
    curve_coefx = curves(index,1:2:end);
    curve_coefy = curves(index,2:2:end);
    coefs = [curve_coefx ; curve_coefy]; crv = spmak(kn,coefs);
    hold on
    points = fnplt(crv);
    handleModCrv = plot(points(1,:),points(2,:),'b');
    axis equal
    
    uknots = unique(kn);
    n = numel(uknots);
    for j = 1:n-2
        crv =  fnrfn(crv,uknots(j+1));
    end
    
    %Create the polygons
    polys=zeros(4,2,n-1);
    M = crv.coefs;  xM = M(1,:);yM = M(2,:); Mold = M;
    active = 1;

    redraw_fig2(opt.nBlocks,opt.bL,opt.bH); %this sets the axis and plot the control polygons

    for j = 1:n-1
        polys(:,:,j) = [xM(1+3*(j-1):1+3*(j-1)+3);yM(1+3*(j-1):1+3*(j-1)+3)]';
    end
    actv = coefs(:,active);
    hold on
    plot(actv(1),actv(2),'.','color',[1 0 0],'MarkerSize',9*csize) %The active point gets the color red
    drawnow
    plot_data(opt.nBlocks)
    for j=1:opt.nBlocks
        compute_curve(polys(:,:,j),100); %compute the inital guesses
    end
    set_axes2    
    drawnow

    %Interactive part
    nit = 0;
    finished=0;
    while finished == 0
        figure(ModWindow)
        [x,y,button]=ginput(1);
        
        if button==3
            nit = nit+1;
            for j=1:numel(coefs(1,:))
                if norm([coefs(1,j),coefs(2,j)]-[x,y])<0.4
                    active=j;
                end
            end
            clf; redraw_fig2(opt.nBlocks,opt.bL,opt.bH) ;fnplt(crv,'b');
            set_axes2

            drawnow
            plot_data(opt.nBlocks);
            drawnow
 
            for j=1:opt.nBlocks
                compute_curve(polys(:,:,j),100);
            end
            set_axes2;

        elseif button==1
            clf
            if opt.Manual
                    valxyz = input('value: ');
                    x = valxyz(1);
                    y = valxyz(3);
            end
            if opt.nBlocks == 1
                if active == 1
                    coefs(1,active) = 0;
                elseif active == 4
                    coefs(1,active) = opt.bW;
                else
                    coefs(1,active) = x;
                end
            elseif opt.nBlocks == 2
                if active == 1
                    coefs(1,active) = -opt.bW;
                elseif active == 6
                    coefs(1,active) = opt.bW;
                else
                    coefs(1,active) = x;
                end
            else
                if active == 1
                    if mod(opt.nBlocks,2) == 0
                        coefs(1,active) = -(opt.bW*opt.nBlocks)/2;
                    else
                        coefs(1,active) = -(opt.bW*(opt.nBlocks-1))/2;
                    end
                elseif active == 2*(opt.nBlocks+1)
                      if mod(opt.nBlocks,2) == 0
                        coefs(1,active) = (opt.bW*opt.nBlocks)/2;
                    else
                        coefs(1,active) = (opt.bW*(opt.nBlocks+1))/2;
                    end
                else
                    coefs(1,active) = x;
                end
            end
            coefs(2,active) = y;
            if active == 1 || active == 4+(opt.nBlocks-2)+(opt.nBlocks-1)+1
                Categ = 1;
            elseif active-1 == 1 || active == 4+(opt.nBlocks-2)+(opt.nBlocks-1)
                Categ = 2;
            else
                if mod(active,2) == 0
                    Categ = 4;
                else
                    Categ = 3;
                end
            end
            if (Categ == 3 || Categ ==4 ) %boxcon

                if Categ == 3
                    midtpt_y = (y+coefs(2,active+1))/2;
  
         
                    xval = Mold(1,active+((active-1)/2));
                    yval = midtpt_y;
                    coefs(1,active+1) = xval+(xval-x);
                    coefs(2,active+1) = yval+(yval-y);
                else
                    midtpt_y = (y+coefs(2,active-1))/2;
                    xval = Mold(1,active+(active-4)/2);
                    yval = midtpt_y;
                    coefs(1,active-1) = xval+(xval-x);
                    coefs(2,active-1) = yval+(yval-y);
                end
            end
            crv=spmak(kn,coefs);
            fnplt(crv,'b')
            uknots = unique(kn);
            n = numel(uknots);
            polys=zeros(4,2,n-1);
            C = [coefs(1,:)',coefs(2,:)'].';
            P = C(:)';
            for j = 1:n-1 %2 segm
                out1 = fnbrk(crv,[uknots(j) uknots(j+1)]);
                koeff = out1.coefs;
                polys(:,:,j)=[koeff(1,:); koeff(2,:)]';
            end
            redraw_fig2(opt.nBlocks,opt.bL,opt.bH);
            set_axes2;
            drawnow
            plot_data(opt.nBlocks);
            set_axes2;
            drawnow
            for j=1:opt.nBlocks
                compute_curve(polys(:,:,j),100);
            end
            set_axes2

            %update spline surface!!
            delete(obj)
            delete(CurveHandles)
            curves(index,:) = coefs(:)';
            
  
            %update surface and curves
            figure(surffig);
            obj = PlotNurbsDesign(curves,opt,kn);
            nit = nit+1;
        
            [CurveHandles,lambdaresi_curves] = PlotCurvesOnDesign(curves,opt);
        else
            if nit == 0
                delete(CurveHandles)
                %update surface and curves
                figure(surffig);
                [CurveHandles,lambdaresi_curves] = PlotCurvesOnDesign(curves,opt);
                curves = [];
            end
            figure(ModWindow);

            close;
            finished = 1;
        end
        
    end

end

end