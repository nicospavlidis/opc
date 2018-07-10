function out = scpp_def_sigma(X)
%Default scaling parameter (sigma) for Gaussian kernel used to estimate kernel/ similarity matrix in SCPP
%OUT = SCPP_DEF_SIGMA(X)
%
% Input:
%	(X): Data matrix 
% 
% Output:
%	(OUT): Default scaling parameter

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

lambdas = sort(eig(cov(X)), 'descend');
dstar = max( min(20, sum(lambdas>=1)), 1);
out = sqrt(mean(lambdas(1:dstar)))* (4/(3*size(X,1)))^(1/(4+dstar));
