function [idx,sumD,U,deg] = spclust(X,K,varargin)
%Spectral Clustering algorithm by Ng, Jordan and West (NIPS 2001)
%[IDX,SUMD,DEG] = SPCLUST(X,K,VARARGIN)
%
% Returns:
%	(IDX): Cluster assignment vector
%	(SUMD): Within-cluster sums of point-to-centroid distances after the application
%		of k-means to the K eigenvectors
%	(U):   Eigenvectors of normalised Laplacian
%	(DEG): Vector containing degree of each vertex
%
%
% Inputs:
%	(X): N-by-D Data matrix or N-by-N Distance Matrix, 
%		or 1-by-(N choose 2) vector of pairwise distances produced by pdist(DataMatrix).
%		By default X is assumed to be a data matrix	
%	(K): Number of clusters
%	Optional Input Arguments specified as Name,Value pairs:
%	(s): Scaling parameter used in similarity matrix:
%			A_{ij} = exp( -norm(X(:,i) - X(:,j))^2/(2*s)^2 )
%		If s=0 (which is the default) the local bandwidth selection rule of 
%		Zelnik-Manor and Perona (NIPS 2004) is used
%	(distmat): If TRUE then X is not a data matrix. Depending on its dimensionality X is treated as either 
%		the output of pdist(DataMatrix) (1-by-(N choose 2)) or as a distance matrix (N-by-N)
%		By default: distmat=FALSE
%	(normalise): If TRUE and if X is a data matrix then all variables are scaled in [-1,1]
%		By default normalise=FALSE
%	(nn): Number of Nearest Neighbours to use in the computation of the local bandwidth selection rule
%		Following Zelnik-Manor and Perona (NIPS 2004) the default value is 7.
%
%References:
%A.Y. Ng, M.I. Jordan, and Y. Weiss. On spectral clustering: Analysis and an algorithm.
%Advances in Neural Information Processing Systems 14, pages 849-856. 2001.
%
%L. Zelnik-Manor and P. Perona. Self-tuning spectral clustering.
%Advances in Neural Information Processing Systems, pages 1601--1608, 2004.

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

s = 0;
ds = 0;
nr = 0;
nn = 7;

if (rem(length(varargin),2)==1)
	error('Optional parameters should always go in pairs');
end

for i=1:2:(length(varargin)-1),
	if ~ischar (varargin{i}),
		error (['Unknown type of optional parameter name (parameter names must be strings).']);
	end
	
	switch lower(varargin{i})
	case 's'
		s = double(varargin{i+1});
	case 'distmat'
		ds = double(varargin{i+1});
	case 'normalise'
		nr = double(varargin{i+1});
	case 'nn'
		nn = double(varargin{i+1});
	otherwise
		error(['Unrecognized parameter: ''' varargin{i} '''']);
	end
end

% Normalise if not distance matrix and normalise==TRUE
if nr & ~ds,
	% remove columns with no variation
	active = max(X) - min(X);
	X = X(:, active>0);
	
	% normalise in [-1,1]
	X = bsxfun(@minus,X, mean(X));
	X = X/max(max(abs(X)));
end

if ~ds,
	if s==0,
		% Zelnik Manor local scaling
		A = squareform(pdist(X));
		s = sort(A);
		% select distance to nn-th nearest neighbor
		s = s(nn+1,:);
		s(s<sqrt(eps)) = sqrt(eps);

		A = exp( -(A.^2)./(s'*s) );
		%assert( sum(diag(A) == zeros(n,1))==n);
	else
		A = squareform( exp(-(pdist(X)./(sqrt(2)*s)).^2) );
	end

% X is output of pdist function
elseif size(X,1) == 1,
	if s==0,
		s = sort(squareform(X));
		% select distance to nn-th nearest neighbor
		s = s(nn+1,:);
		s(s<sqrt(eps)) = sqrt(eps);

		A = exp( -(squareform(X).^2)./(s'*s) );
	else
		A = squareform( exp(-(X./(sqrt(2)*s)).^2) );
	end

% X is n \times n distance matrix
elseif (size(X,1) == size(X,2)) & issymmetric(X),
	if s==0,
		s = sort(X);
		% select distance to nn-th nearest neighbor
		s = s(nn+1,:);
		s(s<sqrt(eps)) = sqrt(eps);

		A = exp( -(X.^2)./(s'*s) );
	else
		A = exp(-(X./(sqrt(2)*s)).^2);
	end

else 
	error('spectral_NgJordan: X matrix provided as input is incompatible\n');
end
A(1:size(A,1)+1:end) = 0;

% to avoid division by zero cap degrees smaller than eps
deg = max([eps*ones(size(A,1),1), sum(A,2)], [], 2);
sqd = 1./sqrt(deg);

% select largest in magnitude eigenvalues
[U,lambdas,flag] = eigs((sqd*sqd').*A, K, 'lm');
U = bsxfun(@rdivide, U, sqrt(sum(U.^2)));
clear A;

% Ng et al. 2002: Normalise rows to unit length
Un = U;
nrm = sqrt(sum(Un.^2,2));
empty = find(nrm<sqrt(eps));
Un(empty,:) = 1/sqrt(K);
nrm(empty) = 1.0;
Un = bsxfun(@rdivide, Un, nrm);

%fprintf('Eigen-Decomposition Complete!\n');
% k-means clustering
[idx,~,sumD] = kmeans(Un,K,'EmptyAction','singleton','Replicates',1);
