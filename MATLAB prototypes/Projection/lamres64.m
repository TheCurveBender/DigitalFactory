
    function lambdares= lamres64(x1,y1,x2,y2,x0,y0,x3,y3)
    %calculates the lambda residual for a Bezier curve with control points:
    %p0=[x0 y0],  p1 = [x1,y1],  p2=[x2,y2], p3=[x3 y3]
    %Uses 64 Gauss points for the integration
  if nargin<5
      x0=0; y0=0; x3=1; y3=0;
  end
cx = [x0,    x1,   x2 ,  x3]; %x coordinates of the 4 control points, 
cy = [y0  , y1 ,   y2 ,   y3]; %y
kn = [0,0,0,0,1,1,1,1];

deg = 3;
%   Gauss points, weights
np      = 64; %-3;    % Number of Gauss points
[gt,gw] = get_gauss_points_and_weights(kn(deg+1:end-3),np);

xsp   = spmak(kn,cx);
ysp   = spmak(kn,cy);
dxsp  = fnder(xsp);
dysp  = fnder(ysp);
d2xsp = fnder(dxsp);
d2ysp = fnder(dysp);
%   Evaluation at the Gauss points
xg    = fnval(xsp,gt);
yg    = fnval(ysp,gt);
dxg   = fnval(dxsp,gt);
dyg   = fnval(dysp,gt);
d2xg  = fnval(d2xsp,gt);
d2yg  = fnval(d2ysp,gt);

%Speed, tangent and curvature

speed = sqrt(dxg.^2 + dyg.^2);
kappa = (dxg.*d2yg - dyg.*d2xg)./speed.^3;

%   Gauss weights for computing the arc length
gsw   = gw.*speed; 

L = sum(gsw);

A = zeros(3,3);
A(1,1) = sum(yg.^2.*gsw);
A(1,2) = -sum(xg.*yg.*gsw);
A(1,3) = -sum(yg.*gsw);
A(2,1) = A(1,2);
A(2,2) = sum(xg.^2.*gsw);
A(2,3) = sum(xg.*gsw);
A(3,1) = A(1,3);
A(3,2) = A(2,3);
A(3,3) = L;

b = zeros(3,1);
b(1) = -sum(kappa.*yg.*gsw);
b(2) = sum(kappa.*xg.*gsw);
b(3) = sum(kappa.*gsw);

lambda = A\b;
lambda1 = lambda(1);
lambda2 = lambda(2);
alpha   = lambda(3);

lambdares = sqrt(...
           (kappa + lambda1*yg - lambda2*xg - alpha).^2*gsw'/...
           (kappa.^2*gsw'));
    end

function [gt,gw] = get_gauss_points_and_weights(kn,np)
Gauss = gauss_quad();
gt    = Gauss{np}(:,1); % Gauss points
gw    = Gauss{np}(:,2); % Gauss weights

[B,IX] = sort(-gt);

gt = [B;gt];
gw = [gw(IX);gw];
kn0 = unique(kn);
mean = (kn0(2:end) + kn0(1:end-1))/2;
diff = (kn0(2:end) - kn0(1:end-1))/2;

gw = gw*diff;
gw = reshape(gw,1,numel(gw));
gt = ones(size(gt))*mean + gt*diff;
gt = reshape(gt,1,numel(gt));

end

function Gauss = gauss_quad()
    
Gauss{64}=[0.0243502926634244325089558 0.0486909570091397203833654
         0.0729931217877990394495429 0.0485754674415034269347991
         0.1214628192961205544703765 0.0483447622348029571697695
         0.1696444204239928180373136 0.0479993885964583077281262
         0.2174236437400070841496487 0.0475401657148303086622822
         0.2646871622087674163739642 0.0469681828162100173253263
         0.3113228719902109561575127 0.0462847965813144172959532
         0.3572201583376681159504426 0.0454916279274181444797710
         0.4022701579639916036957668 0.0445905581637565630601347
         0.4463660172534640879849477 0.0435837245293234533768279
         0.4894031457070529574785263 0.0424735151236535890073398
         0.5312794640198945456580139 0.0412625632426235286101563
         0.5718956462026340342838781 0.0399537411327203413866569
         0.6111553551723932502488530 0.0385501531786156291289625
         0.6489654712546573398577612 0.0370551285402400460404151
         0.6852363130542332425635584 0.0354722132568823838106931
         0.7198818501716108268489402 0.0338051618371416093915655
         0.7528199072605318966118638 0.0320579283548515535854675
         0.7839723589433414076102205 0.0302346570724024788679741
         0.8132653151227975597419233 0.0283396726142594832275113
         0.8406292962525803627516915 0.0263774697150546586716918
         0.8659993981540928197607834 0.0243527025687108733381776
         0.8893154459951141058534040 0.0222701738083832541592983
         0.9105221370785028057563807 0.0201348231535302093723403
         0.9295691721319395758214902 0.0179517157756973430850453
         0.9464113748584028160624815 0.0157260304760247193219660
         0.9610087996520537189186141 0.0134630478967186425980608
         0.9733268277899109637418535 0.0111681394601311288185905
         0.9833362538846259569312993 0.0088467598263639477230309
         0.9910133714767443207393824 0.0065044579689783628561174
         0.9963401167719552793469245 0.0041470332605624676352875
         0.9993050417357721394569056 0.0017832807216964329472961];
         
         
end
    