function [s0,errs0] = find_s0_inflexion(lambda,k,u,du,st,gw)   
    sqrtlambda = sqrt(lambda);
    L          = sum(gw);
    cn    = 1-sqrtlambda*u/(2*k);
    I     = find(cn>1);
    cn(I) = 1;
    I     = find(cn<-1);
    cn(I) = -1;
    am    = acos(cn);    
    I     = find(du>0);
    ap    = am;
    ap(I) = 0;
    am    = 2*ap-am;
    am    = unwrap(am);
    Fam   = ellipticF(am,k^2);   
    s0    = (Fam-st*sqrtlambda)*gw'/L;
    if nargout > 1
        errs0 = (Fam - st*sqrtlambda - s0).^2*gw'/L;
    end 
end 