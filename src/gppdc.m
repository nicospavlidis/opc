function [idx,t] = gppdc(X, K, pphandle, varargin)
%Generic Projection Pursuit Divisive Clustering
%[IDX,T] = GPPDC(X, K, PPHANDLE, VARARGIN)
%
%  [IDX, T] = GPPDC(X, K, PPHANDLE) produces a divisive hierarchical clustering of the 
%  N-by-D data matrix (X) into (K) clusters, using binary partitions produced by a generic
%  projection pursuit function defined in the function handle (PPHANDLE).
%
%  [IDX, T] = GPPDC(X, K, PPHANDLE, 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional parameters
%  in the form of Name,Value pairs. 
%
%  Inputs:
% 	(X): N-by-D data matrix
%	(K): Number of clusters
%	(PPHANDLE): A handle to a function of the type:
%		[V,FVAL,IDX] = PPHANDLE(X,PARAMS), where,
%			(X): X is the data matrix
%			(PARAMS): Is a structure containining all the parameters of the projection pursuit algorithm
%			Returns:
%			(V): Optimal projection matrix/ vector
%			(FVAL): Value of projection index for (v) [default projection index]
%				(If PP algorithm aims to minimise the projection index then -fval must
%				be returned, since at each step cluster with maximum 'split_index' is partitioned)
%			(IDX): Cluster assignment in {1,2}
%
%  Optional Parameters:
%
%  'param' - Structure containing parameter settings employed by PPHANDLE
%	The contents of this structure are provided as input to (PPHANDLE, see above)
%
%  'split_index' - Criterion determining which cluster to split
%	Function Handle: index = split_index(v, X, pars)
%			(v: projection vector, X:data matrix, pars: parameters structure)
%	Cluster with MAXIMUM INDEX is split at each step of the algorithm
%	Two standard choices of split index can be enabled by settgin 'split_index' to 
%	one of the strings below:
%		+ 'fval':    Split cluster whose hyperplane achieves the lowest density integral
%		+ 'size':    Split largest cluster
%	(default: split_index = 'fval' estimated in PPHANDLE)
%
%  'labels' - true cluster labels. Used only to evaluate quality of binary partitions at the end

if nargin < 3 | ~isa(pphandle,'function_handle'),
	error('MATLAB:gppdc:','pphandle must be a function handle');
end

pars = struct();
pars.minsize = 1;
pars.split_index = 'fval';
pars.labels = ones(size(X,1),1);
pars.colours = [];

pars = myparser(X,K,varargin,pars);
[N, dim] = size(X);

% Projection pursuit function: pphandle
[v,f,idx] = pphandle(X,pars);
if size(v,1) ~= dim,
	error('MATLAB:gppdc: Projection pursuit function returns incompatible projection matrix');
end
if isempty(f) | ~isscalar(f) | isinf(f) | isnan(f),
	error('MATLAB:gppdc: Projection pursuit function incorrect projection index');
end
if isempty(idx) | length(idx) ~= N | min(idx)~=1 | max(idx)~=2,
	error('MATLAB:gppdc: Projection pursuit function returns incompatible cluster assignment');
end

%%%%%%%%% HIERARCHICAL ALGORITHM EFFECTIVELY STARTS HERE

% pass{i}: contains binary allocation of observations allocated to i-th tree node based on separating hyperplane 
pass = {};
% location of each node in t
node_index = [1];
% split_index[i]: contains the value of the split index for the tree node with number node_index[i]
split_index = [];

% Generic Projection Pursuit
[optS, pass{1}, split_index(1)] = gpp(X, pphandle, pars);
if isempty(optS),
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

method = regexprep(regexp(func2str(pphandle),'\([a-zA-Z]+\(','match'), '(','');
t = ctree(optS, struct('data',X,'method',strcat('gdc','_',method{1}),'params',pars));

fprintf('Clusters: ');
while length(list) < K,
	fprintf('%i ', length(list));

	% identify which cluster to split next
	[criterion, id] = max(split_index);

	if isinf(criterion),
		fprintf('Divisive process terminated because no cluster can be split further\n');
		fprintf('Clusters Identified: ');
		break;
	end
	% Cluster to be split is no longer a leaf
	t.Node{node_index(id)}.tree_params.leaf = 0;

	% Estimate optimal splits for 2 newly created clusters: Required to determine which cluster to split next
	pos = [length(list)+1, id];
	for i=1:2,
		j = pos(i);
		% Observations allocated to each child depend on sign of pass{pid}
		list{j} = list{id}(pass{id} == (2*i-3));

		% Partition data subset
		[optS, pass{j}, split_index(j)] = gpp(X(list{j},:), pphandle, pars);

		% if PP method fails to identify valid separating hyperplane
		if isinf(optS.fval),
			% Use projection matrix from parent to enable visualisation
			optS = gsep(t.Node{node_index(id)}.v, -inf, [1:length(list{j})]', pphandle, pars);
		end
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
if ~isempty(pars.labels),
	t = perf(t,pars.labels);
end
