function out = mc_spindex(v, X, params, N)
%Default split index for maximum clusterability projection pursuit
%OUT = MC_SPINDEX(X, FVAL, TRDATA)
%
% Default split_index used to select which cluster MCDC partitions at each iteration
% Inputs:
%	(v) Projection vector
%	(X) Data matrix
%	(params) Parameters structure used to extract minsize (minimum cluster size)
%	(N) Total number of observations


[n,d] = size(X);
if N > 2000,
	out = n;
else 
	fval = f_mc(v,X,params.minsize);
	out = 1 - ncfcdf(max(0,n-d-1)*fval/min(n,d+1), min(n,d+1),  max(0,n-d-1), n) ...
		 + 1.e-30*sqrt(n)*fval;
end
