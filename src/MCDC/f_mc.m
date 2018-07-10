function [f, bmin] = f_mc(v,X,minsize)
%Variance ratio clusterability and split point along unit-length (v)
%[F, BMIN] = F_MC(V,X,MINSIZE)
%
% Inputs:
%	(V): Projection vector
%	(X): N-by-D data matrix
%	(MINSIZE): Minimum cluster size
%
% Output:
%	(F): Variance Ratio Clusterability of projected dataset X*(V./norm(V,2))
%	(BMIN): Optimal split point along unit-length vector (V./norm(V,2))

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

n = size(X,1);
x = X * (v./norm(v,2));
proj = sort(x);

cs = cumsum(proj);
mu = cs(n)/n;

% Explicitly enforcing minimum cluster size contstraint
id = [minsize:(n - minsize)]'; 

bc = (cs(id)./id - mu).^2 .*(id./n) + ((cs(n) - cs(id))./(n-id) - mu).^2 .* ((n-id)./n);
[f, M] = max( bc./(var(proj) - bc));
f = -f;

if nargout > 1,
	bmin = 0.5*(proj(M + minsize) + proj(M + minsize-1));
end
