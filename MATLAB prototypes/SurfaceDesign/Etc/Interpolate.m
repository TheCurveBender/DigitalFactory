function [IGlist, spllist, normallist,eSurf,eCurv,X_vec,Y_vec,Z_vec,lambda_residuals] = Interpolate(curves,opt) 
%Obtain the initial guesses and visualize the spline surface that is close
%to being production-ready.
%spline surface is interpolated using cubic interpolation with not-a-knot
%condition

n_total = opt.nCurves*(opt.nICurves-1)+opt.nICurves;
lambda_residuals = []; %lres of spline foliation
yvals = opt.yvals;
%Make start to be a line

t=linspace(0,1,n_total);
eCurv = cell(2,opt.nBlocks);


if opt.nBlocks > 1
    mp_vec = cell(opt.nBlocks-1,2);
    pp_vec = cell(opt.nBlocks-1,2);
    derpp_vec  = cell(opt.nBlocks-1,3);
    for j = 1:opt.nBlocks-1
            mpx = 0.5*(curves(:,4*j+1)+curves(:,4*j+3));
            mpy = 0.5*(curves(:,4*j+2)+curves(:,4*j+4));
            ppx = spline(1:opt.nICurves,mpx);
            ppy = spline(1:opt.nICurves,mpy);
            ppz = spline(1:opt.nICurves,1:opt.nICurves);
            derppx = ppder(ppx);
            derppy = ppder(ppy);
            derppz = ppder(ppz);
            mp_vec{j,1} = mpx;
            mp_vec{j,2} = mpy;
            pp_vec{j,1} = ppx;
            pp_vec{j,2} = ppy;
            derpp_vec{j,1} = derppx;
            derpp_vec{j,2} = derppy;
            derpp_vec{j,3} = derppz;
    end
end

%close all;

%Interpolation

%clf
for j=1:2*(opt.nBlocks+1)*2
     x=1:opt.nICurves;
     y=curves(:,j:j); 
     xx=linspace(1,opt.nICurves,n_total);
     curvesI(:,j:j)=spline(x,y,xx);
end


yvalsext = spline(1:opt.nICurves,yvals,linspace(1,opt.nICurves,n_total));


%close all


X_vec = cell(1,opt.nBlocks); Y_vec = X_vec; Z_vec = X_vec;
colors_vec = cell(1,opt.nBlocks);
endpts_vec = cell(1,opt.nBlocks+1);
el_vec = cell(1,opt.nBlocks);

for j = 1:opt.nBlocks
    X_vec{1,j} = zeros(n_total,n_total);
    Y_vec{1,j} = zeros(n_total,n_total);
    Z_vec{1,j} = zeros(n_total,n_total);
    colors_vec{1,j} = zeros(n_total,3);
    el_vec{1,j} = zeros(2*opt.nICurves,n_total);
end


X1=zeros(n_total,n_total); Y1=X1; Z1=X1; X2=X1; Y2=X2; Z2=X2;X3=X2;Y3=Y2; Z3 = Z2;
colors1=zeros(n_total,3); colors2=colors1;colors3=colors2;colors1(:,2:2)=0.6*ones(n_total,1); colors2(:,2:2)=colors1(:,2:2);colors3(:,2:2)=colors2(:,2:2);
beziercurv1 =[];beziercurv2 = [];beziercurv3=[];endpts1 = [];endpts2 = [];endpts3 = []; endpts4 = [];

splines = [];igs = [];epoints = [];IGepoints = [];Btans = []; IGtans = [];IGlens = [];Blens = [];normals = []; 
el1=zeros(2*opt.nICurves,n_total);el2=zeros(2*opt.nICurves,n_total);el3=zeros(2*opt.nICurves,n_total);

xx=linspace(1,opt.nICurves,n_total);


for j=1:n_total
        splines = [splines curvesI(j,1:end)]; 
        for i = 1:opt.nBlocks
        
        if i == 1 & opt.nBlocks > 1
            p0=curvesI(j:j,1:2);
            p1=curvesI(j:j,3:4);
            p2=curvesI(j:j,5:6);
            p3=curvesI(j:j,7:8); 
            mp = 0.5.*(p2+p3);
            tval = [ppval(derpp_vec{i,1},xx(j)), ppval(derpp_vec{i,2},xx(j)),ppval(derpp_vec{i,3},xx(j))];

            v2 = 3*(mp-p2);
            n = cross([v2(1),v2(2),0],tval);
            n = n/norm(n,2);

            b0 = p0;
            b1 = p1;
            b2 = p2;
            b3 = mp;

        end
        if i > 1 & i < opt.nBlocks
            
            b0 = 0.5*(curvesI(j,4*i-3:4*i-2)+curvesI(j,4*i-1:4*i)); %mp
            b1 = curvesI(j,4*i-1:4*i);
            b2 = curvesI(j,4*i+1:4*i+2);
            b3 = 0.5*(b2+curvesI(j,4*(i+1)-1:4*(i+1)));%mp1
            tval = [ppval(derpp_vec{i,1},xx(j)) ppval(derpp_vec{i,2},xx(j)),ppval(derpp_vec{i,3},xx(j))];
            v2 = 3*(b3-curvesI(j,4*(i+1)-1:4*(i+1)));
            n1 =  cross([v2(1),0,v2(2)],tval);
            n = n1/norm(n1,2);

        end
        if i == opt.nBlocks
            if opt.nBlocks >1
                b0 = 0.5*(curvesI(j,4*i-3:4*i-2)+curvesI(j,4*i-1:4*i));
            else
                b0 = curvesI(j,4*i-3:4*i-2);
            end
            b1 = curvesI(j,4*i-1:4*i);
            b2 = curvesI(j,4*i+1:4*i+2);
            b3 = curvesI(j,4*(i+1)-1:4*(i+1));
            n = [];
        end


        [de_elastiske_kurver,elvec,lambdares1,errtan1,L1,p00,p01,t0,t1]  = FirstGuessD_el(b0,b1,b2,b3,n_total,0);

        el_vec{1,i}(2*j-1:2*j,:) = elvec;
        lambda_residuals = [lambda_residuals lambdares1];
        temp = kron((1-t).^3,b0') + kron(3*(1-t).^2.*t,b1') + kron(3*(1-t).*t.^2,b2') + kron(t.^3,b3');
        X_vec{1,i}(j:j,:) = temp(1:1,:);
        Y_vec{1,i}(:,j:j) = yvalsext; %opt.bL*t';
        Z_vec{1,i}(j:j,:) = temp(2:2,:);
        if lambdares1>0.35|| errtan1>0.1
            colors_vec{1,i}(j,1) = 1;
        end
        if lambdares1<0.2 && errtan1<0.01
            colors_vec{1,i}(j,2)=1;
        end
        if colinear4(b0,b1,b2,b3)
            colors_vec{1,i}(j,:) = [1 0 0];
        end
        endpts_vec{1,i} = [endpts_vec{1,i} b0];
        endpts_vec{1,i+1} = [endpts_vec{1,i+1} b3];

        igs = [igs  de_elastiske_kurver];
        normals = [normals n];
        end   
end

IGlist =reshape(igs,[n_total,opt.nBlocks*7]);
spllist = reshape(splines,[n_total,numel(splines)/n_total]);
normallist = reshape(normals,[n_total,(opt.nBlocks-1)*3]);



eSurf = [];

for j = 1:opt.nBlocks
    S = [X_vec{1,j}; Y_vec{1,j}; Z_vec{1,j}];
    eSurf = [eSurf S];
    
    eCurv{1,j} =  el_vec{1,j}; eCurv{2,j} = colors_vec{1,j};
    if opt.Plot
        hold on
        splotf(S,2,2,2,2,[0.9,0.9,0.9],[0.9,0.9,0.9])
        hold on
        for i=1:n_total
            plot3(el_vec{1,j}(2*i-1:2*i-1,:),Y_vec{1,j}(i:i,:),  el_vec{1,j}(2*i:2*i,:),'color',colors_vec{1,j}(i:i,:),'LineWidth',0.5)
        end
    end
    hold on
end


drawnow;

end