function [f, eigenGap] = f_sc(v, X, pars)
%function [f, eigenGap] = f_sc(v, X, pars)
%
% Returns: 
%	(f) Second smallest eigenvalue of Normalised Laplacian
%	(bmin) Optimal split point along unit-length vector (v./norm(v,2))
%	(eigenGap) Difference between 3rd and 2nd smallest eigenvalues
% Inputs:
% 	(v) Projection matrix (vector storing matrix column-wise)
%	(X) Data matrix
%	(pars) parameter struct containing 
%		(sigma) scaling parameter for Gaussian kernel
%		(weights) weights of micro-clusters (if not used should be empty)
%		(beta) (delta) parameters of similarity transform function
%		(omega) penalty term used to ensure orthonormality of (v)

%function [fval,eigenGap] = f_sc(theta,X,weights,kernel, params) 
%
% fval : Second smallest eigenvalue of Normalised Laplacian
% eigenGap : Difference between 3rd and 2nd smallest eigenvalues

dim = size(X,2);  % dimensionality of original space
pdim = size(v,1)/dim; % dimensionality of the projection space

% Reshape v into projection matrix
if size(v,2) == 1 && pdim > 1,
	V = reshape(v,[dim,pdim]);
else
	V = v;
end
V = bsxfun(@rdivide, V, sqrt(sum(V.^2,1)));

% transformed similarities Eq.(19)
%proj = X*V;
proj = sim_transform(X*V, pars.beta, pars.delta, pars.weights);

% Similarity matrix
W = exp( -(squareform(pdist(proj)).^2)./(2*pars.sigma^2));
if ~isempty(pars.weights),
	W = (pars.weights*pars.weights') .* W;
end
u1 = sqrt(sum(W,2));
A = ( (1./u1)*(1./u1') ).*W; 
%L = eye(size(W,1)) - (u1*u1').*W;
%[L, u1] = LaplacianN(W, pars.weights);
%lambda = sort( eigs(L + u1*u1',2, 1.0e-10) );

N = size(X,1);
if ~isempty(pars.weights),
	N = sum(pars.weights);
end

% Unit-length eigenvector of (L) and (A) corresponding to 
% smallest and largest eigenvalue respectively
v1 = u1./sqrt(sum(u1.^2));
lambda = 1 - sort(eigs(A - v1*v1',2, 'lm'), 'descend');

% Second smallest eigenvalue of L
f = lambda(1);

eigenGap =  lambda(2) - lambda(1);

if pdim>1,
	% triu( V(:,1:end-1)'*V(:,2:end) ).^2 = \sum_{i=1}^{n-1} \sum_{j=i+1}^n <v_{i}, v_{j}>^2
	f = f + pars.omega * sum(sum(triu(V'*V,1).^2));
	%f = pars.omega * sum(sum(triu(V'*V,1).^2));
end

end
