function [x0,y0,errx,erry,xe,ye,Tx,Ty] = find_translation(x,y,st,ela,gw)
    ela.x0  = 0;
    ela.y0  = 0;
    [xe,ye,Tx,Ty] = elastica(ela,st/ela.scale/abs(ela.l));
    L  = sum(gw);
    x0 = (x-xe)*gw'/L;
    y0 = (y-ye)*gw'/L;
    errx = (x-xe-x0).^2*gw'/L;
    erry = (y-ye-y0).^2*gw'/L;
end