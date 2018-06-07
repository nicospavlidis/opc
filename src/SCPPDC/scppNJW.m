function idx = scppNJW(K,v,X,sigma,weights, beta,delta)
%function idx = scppNJW(K,v,X,weights,sigma,params)
%
% Clusters projected data using Ng, Jordan and Weiss (2002) spectral clustering algorithm
%
% Returns:
%	(idx) Cluster assignment vector \in {1,...,K}
%
% Inputs:
%	(K) number of clusters
%	(v) Matrix defining projection subspace
%	(X) Dataset (potentially micro-cluster centers)
%	(weights) Observations per microcluster (empty for no micro-clustering)
%	(sigma) scaling parameter for Gaussian kernel
%	(beta,delta) parameters of similarity transformation function:
%		if empty similarity between projections is based on Euclidean distance

% transformed projections
if nargin==7 & ~isempty(beta) & ~isempty(delta)
	p = sim_transform(X*v, beta, delta, weights);
else
	p = X*v;
end
% Similarity matrix
W = exp( -(squareform(pdist(p)).^2)./(2*sigma^2));

% if micro-clustering has been applied
if nargin >=4 & ~isempty(weights),
	W = (weights*weights') .* W;
end

% zero out diagonal
W(1:size(W,1)+1:end) = 0;

% to avoid division by zero cap degrees smaller than eps
deg = 1./sqrt( max([eps*ones(size(W,1),1), sum(W,2)], [], 2) );

% select largest in magnitude eigenvalues
[U,lambdas,flag] = eigs((deg*deg').*W, K, 'lm');
clear W deg;

% Ng et al. 2002: Normalise rows to unit length
nrm = sqrt(sum(U.^2,2));
empty = find(nrm<sqrt(eps));
U(empty,:) = 1/sqrt(K);
nrm(empty) = 1.0;
U = bsxfun(@rdivide, U, nrm);

%fprintf('Eigen-Decomposition Complete!\n');
% k-means clustering
idx = kmeans(U,K,'EmptyAction','singleton','Replicates',1);

end
