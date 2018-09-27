function [lambda1,lambda2,alpha,residue] = find_lambda(x,y,kappa,w,L)
%
% x,y: curve = coordinates
% kappa: curvature
% w: weights for integration
% L: length
%


A = zeros(3,3);
A(1,1) = sum(y.^2.*w);
A(1,2) = -sum(x.*y.*w);
A(1,3) = -sum(y.*w);
A(2,1) = A(1,2);
A(2,2) = sum(x.^2.*w);
A(2,3) = sum(x.*w);
A(3,1) = A(1,3);
A(3,2) = A(2,3);
A(3,3) = L;

b = zeros(3,1);
b(1) = -sum(kappa.*y.*w);
b(2) = sum(kappa.*x.*w);
b(3) = sum(kappa.*w);

lambda = A\b;
lambda1 = lambda(1);
lambda2 = lambda(2);
alpha   = lambda(3);
if nargout > 3 % we want the residue
    residue = sqrt(sum((kappa + lambda1*y - lambda2*x - alpha).^2.*w)/...
                   sum(kappa.^2.*w));
end
 