function [f, fkde] = pkde(x, proj, h, alpha, eta, epsilon)
%function [f, fkde] = pkde(x, proj, h, alpha, eta, epsilon)
%
% Returns: 
%	(f): Penalised density on hyperplane at point (x)
%	(fkde): Density on hyperplane at (x)
% Inputs:
%	(x) location of split along the projected data
%	(proj) Univariate projection of dataset (X*v)
%	(h) bandwidth parameter
%	(alpha) range over which minimisers of 1D density are sought
%	(eta) Term in penalty function controlling maximum distance between minimisers of the
%		kde and the penalised density integral (recommended value: 0.01)
%	(epsilon) Term in penalty function controlling smoothness (recommended value: 1)

s = std(proj);
L = exp(-0.5)/(sqrt(2*pi)*(h*h)* eta^epsilon);

if length(x)==1,
	fkde = kdeC(x,proj,h);
else
	fkde = fgt_kde(proj,x,h);
end
f = fkde + L* max([zeros(length(x),1), x-alpha*s, -x-alpha*s], [],2).^(1+epsilon);
