function s = find_arclength(dxsp,dysp,np)
Gauss  = gauss_quad();
gt     = Gauss{np}(:,1); % Gauss points
gw     = Gauss{np}(:,2); % Gauss weights
[B,IX] = sort(-gt);
if mod(np,2) % odd number of Gaussian points
    gt = [B(1:end-1);gt];
    gw = [gw(IX(1:end-1));gw];
else         % even number of Gaussian points
    gt = [B;gt];
    gw = [gw(IX);gw];
end
% Integrating from -1 to gauss points:
for i=1:np
    mean = (gt(i) - 1)/2;
    diff = (gt(i) + 1)/2;
    part{i} = zeros(np,2);
    part{i}(:,1) = ones(size(gt))*mean + gt*diff;
    part{i}(:,2)  = gw*diff;
end
[kn,ord]  = fnbrk(dxsp,'knots','order');
kn0  = unique(kn(ord:end-ord+1));

% All Gauss points
s    = zeros(1,(numel(kn0)-1)*np);
l0   = 0;
for i = 1:numel(kn0)-1
    mean = (kn0(i+1) + kn0(i))/2;
    diff = (kn0(i+1) - kn0(i))/2;
    for j = 1:np
        tmp_t  = ones(size(gt))*mean + diff*part{j}(:,1);
        tmp_w  = diff*part{j}(:,2);
        tmp_dx = fnval(dxsp,tmp_t);
        tmp_dy = fnval(dysp,tmp_t);
        tmp_v  = sqrt(tmp_dx.^2+tmp_dy.^2);
        s(np*(i-1)+j) = l0 + sum(tmp_v.*tmp_w);
    end
    tmp_t = ones(size(gt))*mean + gt*diff;
    tmp_w = gw*diff;
    tmp_dx = fnval(dxsp,tmp_t);
    tmp_dy = fnval(dysp,tmp_t);
    tmp_v  = sqrt(tmp_dx.^2+tmp_dy.^2);
    l0    = l0 + sum(tmp_v.*tmp_w);
end
end