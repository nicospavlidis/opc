function out = mc_spindex(v, X, params, N)
%Default split_index used to select which cluster MCDC partitions at each iteration
%OUT = MC_SPINDEX(V, X, PARAMS, N)
%
% Inputs:
%	(V): Projection vector
%	(X): Data matrix
%	(PARAMS): Parameters structure used to extract minsize (minimum cluster size)
%	(N): Total number of observations
%
% Output:
%	(OUT): Splitting index

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

[n,d] = size(X);
if N > 2000,
	out = n;
else 
	fval = f_mc(v,X,params.minsize);
	out = 1 - ncfcdf(max(0,n-d-1)*fval/min(n,d+1), min(n,d+1),  max(0,n-d-1), n) ...
		 + 1.e-30*sqrt(n)*fval;
end
