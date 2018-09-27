function [ela,L]= find_elastica(x,y,dx,dy,kappa,st,gw)
    L   = sum(gw);
    KAP = kappa*gw';
    [lambda1,lambda2,alpha,errlambda] = ...
            find_lambda(x,y,kappa,gw,L);
    lambda    = sqrt(lambda1^2+lambda2^2);
    ela.scale = 1/sqrt(lambda);
    ela.phi   = atan2(lambda2,lambda1);
    ela.l     = L/ela.scale;
    u  = (lambda2*x-lambda1*y)/lambda;
    du = (lambda2*dx-lambda1*dy)/lambda;
    dv = (dx*lambda1 + dy*lambda2)/lambda;
    [beta,inflexion,errbeta] = find_beta(u,dv,lambda,alpha,gw,L);
    ela.inflexion = inflexion;
    if ~inflexion
        if KAP<0
            ela.phi = ela.phi+pi;
            ela.l   = -ela.l;
        end
    end
    Dp   = alpha^2-2*lambda*(beta-1);
    sqDp = sqrt(Dp);
    if inflexion
        Du   = 2*sqDp/lambda;
        k    = Du*sqrt(lambda)/4;
        umax = ( -alpha + sqDp )/lambda;
        u    = umax - u;
        [s0,errs0] = find_s0_inflexion(lambda,k,u,du,st,gw);
    else
        Dm   = alpha^2-2*lambda*(beta+1);
        sqDm = sqrt(Dm);
        Du   = (sqDp-sqDm)/lambda;
        k    = 4*sqrt(lambda)*Du/(4+lambda*Du^2);
        if KAP > 0
            umax = ( -alpha + sqDp )/lambda;
            u    = umax - u;
           [s0,errs0] = find_s0_noinflexion(lambda,k,u,du,st,gw);
        else
            umin = ( -alpha - sqDp )/lambda;
            u = u - umin; 
            [s0,errs0] = find_s0_noinflexionN(lambda,k,u,du,st,gw);
        end
    end
    ela.k  = k;
    ela.s0 = s0;
    [x0,y0,errx,erry,xe,ye,Tx,Ty] = find_translation(x,y,st,ela,gw);
    ela.x0 = x0;
    ela.y0 = y0;
    ela.Tx=Tx;
    ela.Ty=Ty;
    ela.xe=xe;
    ela.ye=ye;
    err = [errlambda,errbeta,errs0,errx+erry];
    ela.residue = err;    
end % find_elastica