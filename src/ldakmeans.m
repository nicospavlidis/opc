function [idx,U,fval,C,iter] = ldakmeans(X,K,varargin)
%LDA-K-means algorithm
%[IDX,U,FVAL,C,ITER] = LDAKMEANS(X,K,VARARGIN)
%
% [IDX,U,FVAL,C,ITER] = ldaKmeans(X,K) produces a flat partitioning of the N-by-D
% data matrix (X) into (K) clusters. The algorithm iteratively performs Linear
% Discriminant Analysis and K-means clustering.
%
% [IDX,U,FVAL,C,FLAG] = ldaKmeans(X,K) returns the cluster assignment, (IDX);
% the projection matrix (U) used to perform dimensionality reduction; the value
% of the projection index (FVAL); the cluster centroids (C); and finally the
% iteration at which the algorithm terminated, (ITER). If LDAKMEANS fails
% to converge ITER=0
% 
% [IDX,U,FVAL,C,ITER] = ldaKmeans(X, K, 'PARAM1',val1, 'PARAM2',val2, ...)
% specifies optional parameters in the form of Name,Value pairs. These are only
% used for visualisation
%
%  'maxit' - Number of iterations to perform (default: 50)
%
%  'ftol' - Tolerance level (default: 1.e-5)
%
%  'verb' - Verbosity. Values greater than 0 enable visualisation during execution
%	Enabling this option slows down the algorithm considerably
%	(default: 0)
%
%  'labels' - true cluster labels. Specifying these enables the computation of performance over 
%	successive iterations and a better visualisation of how clusters are split
%
%  'colours' - Matrix containing colour specification for observations in each of the (K) clusters
%
%Reference:
%C. Ding and T. Li. Adaptive dimension reduction using discriminant analysis and k-means clustering. 
%Proceedings of the 24th International Conference on Machine Learning, pages 521-528, 2007.

pars = struct();
pars.maxit = 50;
pars.ftol = 1.e-5;
pars.verb = 0;
pars.labels = [];
pars.colours = [];

pars = myparser(X,K,varargin, pars);
[n, d] = size(X);

% if labels have not been defined, and if colours are less than K 
% (which can happen because palette assumes binary clustering) set default colours
if max(pars.labels)==1 && size(pars.colours,1) < K,
	pars.colours = palette(K,[]);
end

% Centre and compute 1st principal component
mu = mean(X);
X = bsxfun(@minus, X, mu);
U = pcacomp(X,[1:K-1]);

% Total scatter
St = (n-1)*cov(X);
iter = 0;
Up = zeros(d,K-1);
for i = 1:pars.maxit,

	% LDA-Km(1) step in the paper: Solve for optimal cluster assignment
	% given projection subspace
	proj = X*U;
	idx = kmeans(proj,K,'EmptyAction','singleton','Replicates',1);
	

	% Convergence of projection matrix
	if norm(Up-U,Inf)/(d*K-d) < pars.ftol,
		iter = i;
		break;
	else
		Up = U;
	end

	% Perform LDA on original data using k-means clustering as labels
	%  Sw: Within class scatter
	Sw = zeros(d,d);
	for j=1:K,
		Sw = Sw + min(sum(idx==j)-1, 1) * cov(X(idx==j,:));
	end
	% Between class scatter
	Sb = St - Sw;

	if rcond(Sw) < sqrt(eps),
		warning('Condition number of within cluster sum of squares is too low -- Clusters appear to be defined in different subspaces');
		[U,L] = eigs(pinv(Sw)*Sb,K-1);
	else
		[U,L] = eigs(Sb,Sw,K-1,'LA');
	end
	% Normalise eigenvectors to unit length
	%U = bsxfun(@rdivide, U, sqrt(sum(U.^2)));
	
	fval = trace(U'*Sb*U)/trace(U'*Sw*U);
	if pars.verb > 0, 
		plotter(X,U,fval,idx,i, pars.labels, pars.colours); 
	end
end

C = zeros(K, d);
for j=1:K,
	C(j,:) = mean(X(idx==j,:), 1) + mu;
end


function plotter(X, U, fval, idx, iter, labels, colours)

X = X*U;
K = max(idx);

hFig = figure(1);
clf;
set(hFig, 'Position', [0 0 1600 800]);
if max(labels) > 1,

	subplot(1,2,1)
	% ensure reproducibility of colours
	hold on;
	l = unique(labels);
	for i=1:length(l),

		s = find(labels==l(i));
		plot(X(s,1), X(s,2),'bo','MarkerSize',3, 'Color', colours(l(i),:))
		% , 'MarkerFaceColor', colours(l(i),:));
	end
	hold off;
	axis([min(X(:,1)), max(X(:,1)), min(X(:,2)), max(X(:,2))]);
	title('Actual clusters');

	subplot(1,2,2)
	% Set colours for clusters (K might differ from true number of clusters)
	if max(labels) == K,
		colours2 = colours;
		
		T = mycrosstab(labels,idx);
		% number of observations allocated to each (estimated) cluster
		per_clust = sum(T,1);
		
		% available colours
		avail = [1:K];
		[~,order] = sort(per_clust, 'descend');
		for i=1:K,
			% identify the true cluster with maximum number of observations
			% allocated to i-th biggest estimated cluster
			[~,selected] = max( T(avail, order(i)) );

			% set colour of estimated cluster to be same as maximal true cluster
			colours2(order(i),:) = colours(avail(selected), :);

			avail = avail([1:selected-1, selected+1:end]);
		end

	else
		colours2 = palette(K,[]);
	end
	hold on;
	for i=1:K,
		s = find(idx==i);
		plot(X(s,1), X(s,2),'bo','MarkerSize',3, 'Color', colours2(l(i),:));
	end
	str = sprintf('Estimated clusters',fval);
	hold off;

	axis([min(X(:,1)), max(X(:,1)), min(X(:,2)), max(X(:,2))]);
	c = cluster_performance(idx, labels);

	if ~isOctave(),
		title(str);
		suptitle(sprintf('iter: %i fval: %1.3f Pur: %1.2f NMI %1.2f Ad.Rand: %1.2f', ...
			iter, fval, c.Purity, c.NMI, c.AdjRand));
	else
		title(sprintf('iter: %i fval: %1.3f Pur: %1.2f NMI %1.2f', iter, fval, c.Purity, c.NMI));
	end

	drawnow;
else
	%keyboard
	hold on;
	for i=1:K,
		s = find(idx==i);
		plot(X(s,1), X(s,2),'bo','MarkerSize',3, 'Color', colours(i,:));
	end
	str = sprintf('iter %i fval=%1.5f',iter, fval);
	title(str);
	hold off;
	axis([min(X(:,1)), max(X(:,1)), min(X(:,2)), max(X(:,2))]);
	drawnow;
end
