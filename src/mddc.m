function [idx,t] = mddc(X, K, varargin)
%Minimum Density Divisive Clustering 
%[IDX,T] = MDDC(X, K, VARARGIN)
%
% [IDX, T] = MDDC(X, K) produces a divisive hierarchical clustering of the
% N-by-D data matrix X into (a maximum of) K clusters. This algorithm uses a
% hierarchy of binary partitions each splitting the observations with the
% hyperplane with minimum density integral.  The algorithm can return fewer
% clusters if no valid hyperplane separators are found.
%
%  [IDX, T] = MDDC(X, K) returns the cluster assignment, IDX, and the  binary tree
%  (T) containing the cluster hierarchy
%
%  [IDX, T] = MDDC(X, K, 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional parameters
%  in the form of Name,Value pairs. 
%
%  OPTIONAL PARAMETERS:
%  'v0' - Initial projection vector(s)
%	Function handle: v0(X,P) returns D-by-S matrix of initial projection vectors
%   	(default: v0 = @(y,p)(pca(y,'NumComponents',1)) -- 1st principal component)
%
%  'bandwidth' - Bandwidth parameter
%	Function Handle: bandwidth(X,pars) returns bandwidth (positive scalar)
%	(default: bandwidth =0.9* sqrt(eigs(cov(X),1)) * size(X,1)^(-0.2))
%
%  'split_index' - Criterion determining which cluster to split
%	Function Handle: index = split_index(v, X, pars)
%			(v: projection vector, X:data matrix, pars: parameters structure)
%	Cluster with MAXIMUM INDEX is split at each step of the algorithm
%	Three standard choices of split index can be enabled by settgin 'split_index' to 
%	one of the strings below:
%		+ 'fval':    Split cluster whose hyperplane achieves the lowest density integral
%		+ 'size':    Split largest cluster
%		+ 'rdepth':  Split cluster with maximum relative depth
%	(default: split_index = 'size')
%
%  'minsize' - Minimum cluster size (integer)
%	(default minsize = 1)
%
%  'alphamin' - The minimum ALPHA over which MDHs are sought: 
%	[mean(X*V) - ALPHA*std(X*V), mean(X*V) + ALPHA*std(X*V)]. 
%	ALPHA starts from (alphamin) and increases by 0.1 every (maxit) iterations until (alphamax) is reached.
%	(default: 0)
%
%  'alphamax' - The maximum ALPHA over which MDHs are sought is 
%	[mean(X*V) - ALPHA*std(X*V), mean(X*V) + ALPHA*std(X*V)]. 
%	ALPHA starts from (alphamin) and increases by 0.1 every (maxit) iterations until (alphamax) is reached.
%	(default: 1)
%
%  'maxit' - Number of BFGS iterations to perform for each value of alpha
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
%Reference:
%N.G. Pavlidis, D.P. Hofmeyr and S.K. Tasoulis. Minimum density hyperplanes.
%Journal of Machine Learning Research, 17(156):1–33, 2016.
%http://jmlr.org/papers/v17/15-307.html.

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if nargin < 2, 
	error('Data matrix and number of clusters need to be specified\n'); 
end;

% Set default parameters
pars = struct();
pars.v0 = @(x,p)(pcacomp(x,1));
pars.bandwidth = @(x,p)(0.9*size(x,1)^(-0.2)*std(x* pcacomp(x,1)));
pars.split_index = 'size';
pars.alphamax = 1.0;
pars.alphamin = 0.0;
pars.minsize = 1;
pars.maxit = 50;
pars.ftol = 1.e-5;
pars.labels = [];
pars.colours = [];
pars.verb = 0;

% We strongly recommend not to modify these
pars.eta = 0.01;
pars.epsilon = 1-1.0e-6;

pars = myparser(X,K,varargin,pars);
[N, dim] = size(X);

assert(pars.alphamin>=0 & pars.alphamax>pars.alphamin);

assert(isa(pars.bandwidth,'function_handle'), ...
	'bandwidth in divisive hierarchical algorithms must be set to a function handle');
assert(isa(pars.v0,'function_handle'), ...
	'v0 in divisive hierarchical algorithms must be set to a function handle');

nc = max(pars.labels);
% Colours
pars.colours = palette(nc, pars.colours);

% Verify parameter settings
if pars.verb >0,
	fprintf('\nPARAMETER VALUES:\n');
	fprintf('v0: %s\n', func2str(pars.v0));
	fprintf('bandwidth: %s\n', func2str(pars.bandwidth));
	if ischar(pars.split_index),
		fprintf('split_index: %s\n', pars.split_index);
	else
		fprintf('split_index: %s\n', func2str(pars.split_index));
	end
	fprintf('alpha range: (%1.3f, %1.3f)\n', pars.alphamin, pars.alphamax);
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
[opthp, pass{1}, split_index(1)] = mdpp(X, pars, pars.labels, pars.colours);
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
t = ctree(opthp, struct('data',X,'method','mddc','params',pars));

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
		[opthp, pass{j}, split_index(j)] = mdpp(X(list{j},:), pars, la, cl);
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
