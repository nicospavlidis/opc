function [idx,t] = scppdc(X, K, varargin)
%Spectral Clustering Projection Pursuit Divisive Clustering
%[IDX,T] = SCPPDC(X, K, VARARGIN)
%
% [IDX, T] = SCPPDC(X, K) produces a divisive hierarchical clustering of the
% N-by-D data matrix (X) into (a maximum of) K clusters. This algorithm uses a
% hierarchy of binary partitions each splitting the observations by identifying
% the linear subspace that minimises the second smallest eigenvalue of the
% normalised Laplacian and then splitting by applying spectral clustering on
% the projected data.
% The algorithm can return fewer clusters if no subspace is found that returns 
% a valid clustering.
%
%  [IDX, T] = SCPPDC(X, K) returns the cluster assignment, IDX, and the  binary tree
%  (T) containing the cluster hierarchy
%
%  [IDX, T] = SCPPDC(X, K, 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional parameters
%  in the form of Name,Value pairs. 
%
%  OPTIONAL PARAMETERS:
%  'v0' - Initial projection matrix
%	Function handle: v0(X) returns D-by-S projection matrix
%   	(default: v0 = @(y)(pca(y,'NumComponents',2)) -- first 2 principal components)
%
%  'sigma' - Scale parameter for Gaussian kernel used to estimate similarity matrix
%	Function Handle: sigma(X,params) returns sigma parameter (positive scalar)
%
%  'split_index' - Criterion determining which cluster to split
%	Function Handle: index = split_index(v, X, pars)
%			(v: projection matrix, X:data matrix, pars: parameters structure)
%	Cluster with MAXIMUM INDEX is split at each step of the algorithm
%	Three standard choices of split index can be enabled by settgin 'split_index' to 
%	one of the strings below:
%		+ 'size':    Split largest cluster
%		+ 'fval':    Split cluster whose normalised Laplacian has smallest second eigenvalue
%	(default: split_index = 'size')
%
%  'NumMicroClust' - Number of microclusters
%	(default: 200)
%
%  'minsize' - Minimum cluster size (integer)
%	(default: size(X,1)/(5*K))
%
%  'beta' - Initial/ maximum value of BETA (used in similarity transformation)
%	(default: 3.0)
%
%  'betaMin' - The minimum BETA used in similarity function transformation: 
%	BETA starts from (BETA) and decreaes by 0.2 every until a valid subspace is found.
%	(default: 0.5)
%
%  'maxit' - Number of BFGS iterations to perform
%	(default: 50)
%
%  'ftol' - Stopping criterion for change in objective function value over consecutive iterations
%	(default: 1.e-5)
%
%  'verb' - Verbosity. Values greater than 0 enable visualisation during execution
%	Enabling this option slows down the algorithm considerably
%	(default: 0)
%
%  'labels' - true cluster labels. Specifying these enables the computation of performance over 
%	successive iterations and a better visualisation of how clusters are split
%
%  'colours' - Matrix containing colour specification for observations in different clusters
%	Number of rows must be equal to the number of true clusters (if 'labels' has been specified) or equal to 2.
%
%Reference
%D.P. Hofmeyr, N.G. Pavlidis, and I.A. Eckley. Minimum spectral connectivity projection pursuit. 
%Statistics and Computing, forthcoming.

% Set all parameter values to defaults
pars = schp.def_opar();
pars.minsize = floor(size(X,1)/(5*K));
pars.split_index = 'size';
pars.labels = [];
pars.colours = [];

pars = myparser(X,K,varargin,pars);
[N, dim] = size(X);

nc = max(pars.labels);
% Colours
pars.colours = palette(nc, pars.colours);

% Parameter checking
assert(pars.betaMin <= pars.beta,'Incorrect range of beta');
assert(numel(pars.omega)==1 & pars.omega>=0);
assert(numel(pars.NumMicroClust)==1 & pars.NumMicroClust>10*K);



% User defined parameters debugging
if pars.verb >0,
	fprintf('\nPARAMETER VALUES:\n');
	fprintf('v0: %s\n',func2str(pars.v0)); 
	fprintf('sigma: %s\n', func2str(pars.sigma));
	if isa(pars.split_index,'function_handle'),
		fprintf('split_index: %s\n', func2str(pars.split_index));
	else
		fprintf('split_index: %s\n', pars.split_index);
	end
	fprintf('Number of microclusters: %i\n', pars.NumMicroClust);
	fprintf('beta range: (%1.3f, %1.3f)\n', pars.betaMin, pars.beta);
	fprintf('minsize: %i\n', pars.minsize);
	fprintf('maxit: %i\n', pars.maxit);
	fprintf('ftol: %e\n', pars.ftol);
end


%%%%%%%%% HIERARCHICAL ALGORITHM EFFECTIVELY STARTS HERE

% pass{i}: contains binary allocation of observations allocated to i-th tree node based on separating hyperplane 
pass = {};
% location of each node in t
node_index = [1];
% split_index[i]: contains the value of the split index for the tree node with number node_index[i]
split_index = [];
[optS, pass{1}, split_index(1)] = scpp(X, pars, pars.labels, pars.colours);
if isempty(optS),
	% No valid clusters found: Exit
	idx = ones(N,1);
	t = ctree;
	return;
end

% list{i}: contains indices of observations allocated to node_index[i] tree node
list = {};
list{1} = [1:N]';

% Initialise binary tree
optS.idx = pass{1}.*list{1};
optS.tree_params.nodeid = 1;
optS.tree_params.leaf = 1;
t = ctree(optS, struct('data',X,'method','scpp','params',pars));

fprintf('Clusters: ');
while length(list) < K,
	fprintf('%i ', length(list));

	% identify which cluster to split next
	[criterion, id] = max(split_index);

	if isinf(criterion),
		fprintf('Divisive process terminated because no cluster can be split further\n');
		fprintf('%i Clusters Identified\n', length(list));
		break;
	end
	% Cluster to be split is no longer a leaf
	t.Node{node_index(id)}.tree_params.leaf = 0;

	% Estimate optimal splits for 2 newly created clusters: Required to determine which cluster to split next
	pos = [length(list)+1, id];
	for i=1:2,
		j = pos(i);
		% Observations allocated to each child depend on sign of pass{pid}
		list{j} = list{id}(pass{id} == (-1)^i);

		la = pars.labels(list{j});
		cl = ifelse(nc>1, pars.colours(unique(la),:), pars.colours);

		% Partition data subset
		[optS, pass{j}, split_index(j)] = scpp(X(list{j},:), pars, la, cl);
		optS.idx = pass{j}.*list{j};
		optS.tree_params.nodeid = nnodes(t)+1;
		optS.tree_params.leaf = 1;

		% extend the tree by splitting the node with number: node_index(id)
		t = t.addnode(node_index(id), optS);

		node_index(j) = optS.tree_params.nodeid;
	end
end
fprintf('%i\n', length(list));

% Cluster assignment
idx = tree2clusters(t);
t.cluster = idx;

% Compute performance in terms of Success Ratio and Purity for
% all separating hyperplanes (tree nodes)
if nc > 1,
	t = perf(t,pars.labels);
end
