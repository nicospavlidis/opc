function [idx,C,sumd] = dr_kmeans(X,K,e)
%Dimensionality-reduction for k-means algorithm
%[IDX,C,SUMD] = DR_KMEANS(X, K, EPSILON)
%
% IDX = DR_KMEANS(X, K, EPSILON) produces a clustering of the N-by-D data
% matrix (X) into (K) clusters, using the dimensionality reduction k-means
% algorithm (Algorithm 4 in Feldman et al. (2018)). This algorithm first
% reduces the size and dimensionality of the data using the first M right
% singular vectors of the Singular Value Decomposition of X. M is chosen
% internally by the algorithm to ensure that the resulting clustering will be
% an (1+EPSILON) approximation of the clustering in the full dimensional space,
% where EPSILON must be in (0,0.5). (see Corollary 23 in Feldman et al. (2018))
%
% Returns:
%	(IDX): Cluster assignment vector
%	(C): Cluster centroids stored as rows of K-by-D matrix (C)
%	(SUMD): Within-cluster sums of point-to-centroid distances
%
% Inputs:
%	(X): N-by-D data matrix
%	(K): Number of clusters
%	(EPSILON): A real number in the interval (0,0.5) that 
%		determines the approximation error
%
%References: 
%
%D. Feldman, M. Schmidt and C. Sohler. Turning Big Data into Tiny Data:
%Constant-size Coresets for K-means, PCA and Projective Clustering.
%Twenty-fourth Annual ACM-SIAM Symposium on Discrete Algorithms (SODA'13),
%pages 1434--1453, 2013.
%
%D. Feldman, M. Schmidt and C. Sohler. Turning Big Data into Tiny Data:
%Constant-size Coresets for K-means, PCA and Projective Clustering.
%arXiv:1807.04518 [cs.DS], 2018.

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------


if nargin < 3,
	e = 1/3;
elseif e<=0 | e>1/3,
	error('EPSILON must be in (0,1/3].');
end

if nargin < 2 | isempty(K) | K<=0,
	error('K must be a positive integer value.');
end

if nargin < 1 | isempty(X) | size(X,1)<K,
	error('X must have more rows than the number of clusters.');
end

[n,d] = size(X);

% Compute coreset size
m = K + ceil(72*K/(e^2)) - 1;

if m >= min(n,d),
	warning('Choice of K and EPSILON do not allow dimensionality reduction. Using kmeans on full dimensional dataset!');
end

% Compute coresets
[U,S,V] = svds(X,m,'largest')
Am = S*V';

[idx,C] = kmeans(Am,K,'EmptyAction','singleton','Replicates',1);

if nargout == 3,
	sumd = zeros(K,1);
	for i=1:K,
		sumd(i) = norm( bsxfun(@minus, X(idx==i,:), C(i,:)), 'fro')^2;
	end
end
