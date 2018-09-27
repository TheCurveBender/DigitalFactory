function [gt,gw] = get_gauss_points_and_weights(kn,np)
Gauss = gauss_quad();
gt    = Gauss{np}(:,1); % Gauss points
gw    = Gauss{np}(:,2); % Gauss weights
[B,IX] = sort(-gt);
if mod(np,2) % odd number of Gaussian points
    gt = [B(1:end-1);gt];
    gw = [gw(IX(1:end-1));gw];
else         % even number of Gaussian points
    gt = [B;gt];
    gw = [gw(IX);gw];
end
% All Gauss points
kn0 = unique(kn);
mean = (kn0(2:end) + kn0(1:end-1))/2;
diff = (kn0(2:end) - kn0(1:end-1))/2;
gw = gw*diff;
gw = reshape(gw,1,numel(gw));
gt = ones(size(gt))*mean + gt*diff;
gt = reshape(gt,1,numel(gt));
end