function [idx,t] = ncutdc(X, K, varargin)
%Minimum Normalised Cut Divisive Clustering
%[IDX,T] = NCUTDC(X, K, VARARGIN)
%
%  [IDX, T] = NCUTDC(X, K) produces a divisive hierarchical clustering of the 
%  N-by-D data matrix X into K clusters. This algorithm uses a hierarchy of
%  binary partitions each splitting the observations with the hyperplane that
%  minimises the normalised cut criterion. The algorithm can return fewer 
%  clusters if no valid hyperplane separators are found.
%
%  [IDX, T] = NCUTDC(X, K) returns the cluster assignment, IDX, and the  binary 
%  tree (T) containing the cluster hierarchy
%
%  [IDX, T] = NCUTDC(X, K, 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional parameters
%  in the form of Name,Value pairs. 
%
%  'v0' - Function handle. v0(X) returns D-by-S matrix of initial projection vectors
%   	(default: v0 = @(X)(pca(X,'NumComponents',1)) -- 1st principal component)
%
%  'sigma' - Function handle. sigma(X,pars) returns scaling parameter (sigma) as function of data matrix
%	(default: sigma = 100*sqrt(l)*N^(-0.2), where l = max(eig(cov(X))); )
%
%  'split_index' - Criterion determining which cluster to split
%	Function Handle: index = split_index(v, X, pars)
%			(v: projection vector, X:data matrix, pars: parameters structure)
%	Cluster with MAXIMUM INDEX is split at each step of the algorithm
%	Two standard choices of split index can be enabled by setting 'split_index' to 
%	one of the strings below:
%		+ 'fval':    Split cluster whose hyperplane achieves the lowest normalised cut value
%		+ 'size':    Split largest cluster
%	(default: split_index = 'fval')
%
%  'minsize' - Minimum cluster size (integer)
%	(default minsize = 1)
%
%  'maxit' - Number of BFGS iterations to perform for each value of alpha (default: 50)
%
%  'ftol' - Stopping criterion for change in objective function value over consecutive iterations
%	(default: 1.e-7)
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
%Reference:
%D.P. Hofmeyr. Clustering by minimum cut hyperplanes. IEEE Transactions on Pattern
%Analysis and Machine Intelligence, 39(8):1547–1560, 2017.

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

% Set default parameters
pars = struct();
pars.v0 = @(x,p)(pcacomp(x,1));
pars.split_index = 'fval';
pars.sigma = @(x,p)(ncut_sigma(x,100));
pars.minsize = 1;
pars.maxit = 50;
pars.ftol = 1.e-7;
pars.verb = 0;
pars.labels = [];
pars.colours = [];

pars = myparser(X,K,varargin,pars);

assert(isa(pars.sigma,'function_handle'), ...
	'scale parameter (sigma) must be set to a function handle');
assert(isa(pars.v0,'function_handle'), ...
	'v0 in divisive hierarchical algorithms must be set to a function handle');

[N, dim] = size(X);

nc = max(pars.labels);
% Colours
pars.colours = palette(nc, pars.colours);

if pars.verb >0,
	fprintf('\nPARAMETER VALUES:\n');
	fprintf('v0: %s\n', func2str(pars.v0));
	fprintf('sigma: %s\n', func2str(pars.sigma));
	if ischar(pars.split_index),
		fprintf('split_index: %s\n', pars.split_index);
	else
		fprintf('split_index: %s\n', func2str(pars.split_index));
	end
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
[opthp, pass{1}, split_index(1)] = ncutpp(X, pars, pars.labels, pars.colours);
if isempty(opthp),
	idx = ones(N,1);
	t = ctree;
	return;
end

% list{i}: contains indices of observations allocated to node_index[i] tree node
list = {};
list{1} = [1:N]';

% Initialise binary tree
opthp.idx = pass{1}.*list{1};
opthp.tree_params.nodeid = 1;
opthp.tree_params.leaf = 1;
t = ctree(opthp, struct('data',X,'method','ncutdc','params',pars));

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
		[opthp, pass{j}, split_index(j)] = ncutpp(X(list{j},:), pars, la, cl);
		opthp.idx = pass{j}.*list{j};
		opthp.tree_params.nodeid = nnodes(t)+1;
		opthp.tree_params.leaf = 1;

		% extend the tree by splitting the node with number: node_index(id)
		t = t.addnode(node_index(id), opthp);

		node_index(j) = opthp.tree_params.nodeid;
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
