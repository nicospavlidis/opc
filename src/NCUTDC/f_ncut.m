function [f,bmin,idx] = f_ncut(v, Data, pars)
%Function value, split point, and cluster assignment for minimum normalised cut projection index
%[F,BMIN,IDX] = F_NCUT(V, DATA, PARS)
%
% Inputs:
% 	v: Projection vector
%	Data: Data matrix
%	pars: parameter struct containing 
%		minsize: minimum cluster size and
%		sigma: Scaling parameter for Laplace kernel
%
% Outputs: 
%	f: Value of normalised cut for Data*(v./norm(v,2))
%	bmin: Optimal split point along unit-length vector (v./norm(v,2))
%	idx: Resulting binary assignment 

if isa(pars.sigma,'function_handle')
	sigma = pars.sigma(Data);
else
	sigma = pars.sigma;
end

minN = pars.minsize;

normV = norm(v,2);
proj = Data*( v./(normV*sigma) );
x = sort(proj);
n = length(x);
CP = exp(x(1) - x);
CSPi = cumsum(1./CP);

c = cumsum( CP(n:-1:1) );
S = c((n-1):-1:1) .* CSPi(1:n-1);
updown = cumsum(CP.*CSPi);
downup = cumsum( [c((n-1):-1:1); 0]./CP );
DL = updown+downup;
r = S.*( 1./DL(1:(n-1)) +  1./(DL(n) - DL(1:(n-1))) );
[f,w] = min( r(minN:(n-minN)) );

if nargout > 1,
	bmin = 0.5*(x(w + minN - 1) + x(w+minN))*sigma;
end
if nargout > 2,
	idx = (proj > x(w+minN-1)) + 1;
end
