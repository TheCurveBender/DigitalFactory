function [s0,errs0] = find_s0_noinflexionN(lambda,k,u,du,st,gw)
    sqrtlambda = sqrt(lambda);
    L          = sum(gw);
    sn2    = 1 - ( 1 - k*sqrtlambda/2*u).^2;
    sn2    = sn2/k^2;
    I      = find(sn2>1);
    sn2(I) = 1;
    I      = find(sn2<0);
    sn2(I) = 0;
    am     = asin(sqrt(sn2));
    I     = find(du>0);
    ap    = am;
    ap(I) = 0;
    am    = 2*ap-am;
    am    = unwrap(2*am)/2;
    Fam   = ellipticF(am,k^2);
    s0    = (k*Fam + st*sqrtlambda)*gw'/L;
    if nargout > 1
        errs0 = (k*Fam + st*sqrtlambda - s0).^2*gw'/L;
    end
end 