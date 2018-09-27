function [beta,inflexion,residue]=find_beta(u,dv,lambda,alpha,w,L)

    beta      = (dv - lambda/2*u.^2 - alpha*u)*w'/L;
    valmin    = beta - alpha^2/(2*lambda); % minimum for parabola.
    inflexion = valmin > -1; % true if elastica with inflection points
    if nargout > 2 % we want the residue
        residue = (dv - lambda/2*u.^2 - alpha*u - beta).^2*w';
        residue = residue/L;
    end  
end % 