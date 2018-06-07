function [f, bmin] = f_mc(v,X,minsize)
%Computes variance ratio clusterability and split point along unit-length projection vector
%[F, BMIN] = F_MC(V,X,MINSIZE)
%
% Inputs:
%	(v) Projection vector
%	(X) N-by-D data matrix
%	(minsize): minimum cluster size
%
% Outputs:
%	(f) Variance Ratio Clusterability of projected dataset X*v
% 	(bmin) Optimal split point along unit-length vector (v./norm(v,2))


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
