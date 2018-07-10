function [idx,t] = pddp(X,K, varargin)
%Principal Direction Divisive Partitioning algorithm
%[IDX,T] = PDDP(X, K, VARARGIN) 
%
% [IDX, T] = PDDP(X, K) produces a divisive hierarchical clustering of the
% N-by-D data matrix (X). This algorithm uses a hierarchy of binary partitions
% each splitting the observations by first projecting onto the first principal
% component and then splitting at the mean of the projected data
%
%  [IDX,T] = PDDP(X,K) returns the cluster assignment, IDX, and the  binary 
%  tree (T) containing the cluster hierarchy
%
%  [IDX, T] = PDDP(X, 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional parameters
%  in the form of Name,Value pairs. 
%
%  'split_index' - Criterion determining which cluster to split	next
%	Function Handle: index = split_index(v, X, pars)
%			(v: projection vector, X:data matrix, pars: parameters structure)
%	Cluster with MAXIMUM INDEX is split at each step of the algorithm
%	Two standard choices of split index can be enabled by setting 'split_index' to 
%	one of the strings below:
%		+ 'scatter': Split cluster with largest total scatter value
%		+ 'size':    Split largest cluster
%	(default: split_index = 'scatter')
%
%  'minsize' - Minimum cluster size (integer)
%	(default minsize = 1)
%
%  'labels' - true cluster labels. Only used for performance assessment
%
%Reference:
%D. Boley. Principal Direction Divisive Partitioning. Data Mining and Knowledge Discovery, 2(4):325-344, 1998.

% Set default parameters
pars = struct();
pars.minsize = 1;
pars.split_index = @(v,x,pars)(total_scatter(x));
pars.labels = ones(size(X,1),1);
pars.colours = [];

% PARSING
pars = myparser(X,K,varargin, pars);
[N, dim] = size(X);

%%%%%%%%% HIERARCHICAL ALGORITHM EFFECTIVELY STARTS HERE

% pass{i}: contains binary allocation of observations allocated to i-th tree node based on separating hyperplane 
pass = {};
% location of each node in t
node_index = [1];
% split_index[i]: contains the value of the split index for the tree node with number node_index[i]
split_index = [];
[opthp, pass{1}, split_index(1)] = gpp(X, @pddppp, pars);
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
t = ctree(opthp, struct('data',X,'method','PDDP','params',pars));

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

		% Partition data subset
		[opthp, pass{j}, split_index(j)] = gpp(X(list{j},:), @pddppp, pars);
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
if ~isempty(pars.labels),
	t = perf(t,pars.labels);
end
