function [f, df, split] = f_df_mc(v,X,minsize)
%Function value, derivative, and split point along projection vector of variance ratio clusterability projection index
%[F, DF, SPLIT] = F_DF_MC(V,X,PARS)
%
% Returns: 
%	(f) Variance Ratio Clusterability projection index
%       (df) Derivative of (f) w.r.t. projection vector (v)
%	(bmin) split point along (v)
%
% Inputs:
% 	(v) Projection vector
%	(X) N-by-D Data matrix
%	(minsize): minimum cluster size

n = size(X,1);
nv = norm(v,2);

proj = X * (v./nv);
[srt,index] = sort(proj);

cs = cumsum(srt);
id = [minsize:(n - minsize)]';
assert(~isempty(id),'f_df_mc: Error because the cluster being split is too small');

mu = mean(proj);
% BC is the numerator of the objective function (defined as BC in Eq. (10))
bc = (cs(id)./id - mu).^2 .*(id./n) + ((cs(n) - cs(id))./(n-id) - mu).^2 .* ((n-id)./n);

Var = var(proj);
%VRS = bc./(Var - bc);

[f, M] = max( bc./(Var - bc) );
f = -f;

if nargout > 1,
	% split point of one-dimensional sample
	split = 0.5*(srt(M+minsize) + srt(M + minsize-1));
	pos = M+minsize-1;

	m1 = cs(pos)/pos;
	m2 = (cs(n)-cs(pos))/(n-pos);
	%mu = m1*(M/n) + m2*((n-M)/n);

	%beta = dVR in Eq. (7). First construct beta_k from Eq. (9)
	beta = (2/(n-1))*(proj' - mu)*bc(M);

	% next add up alpha_k from Eq. (8)
	ixs = (proj <= split);
	beta(ixs) = ((2./n)*(m1-mu)*Var - beta(ixs))./(Var - bc(M))^2;

	%ixs = find(proj > split);
	beta(~ixs) = ((2./n)*(m2-mu)*Var - beta(~ixs))./(Var - bc(M))^2;

	df = -beta *(X/nv - (X*v)*(v'./nv^3));
end
