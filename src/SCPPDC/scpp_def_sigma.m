function out = scpp_def_sigma(X)
%Default setting for scaling parameter (sigma) for Gaussian kernel used to estimate kernel/ similarity matrix
%out = scpp_def_sigma(X)
% 
% Returns:
%	(out) Default scaling parameter
% Input:
%	(X) Data matrix 

lambdas = sort(eig(cov(X)), 'descend');
dstar = max( min(20, sum(lambdas>=1)), 1);
out = sqrt(mean(lambdas(1:dstar)))* (4/(3*size(X,1)))^(1/(4+dstar));
