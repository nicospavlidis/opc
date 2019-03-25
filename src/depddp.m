function [idx,t] = depddp(X, K, varargin)
%density-enhanced Principal Direction Divisive Partitioning algorithm
%[IDX,T] = DEPDDP(X, VARARGIN)
%
% [IDX, T] = DEPDDP(X, K) produces a divisive hierarchical clustering of the
% N-by-D data matrix (X) into (K) clusters. This algorithm uses a hierarchy of
% binary partitions each splitting the observations by first projecting onto
% the first principal component and then identifying the lowest local minimum
% of the 1D KDE constructed from the projected data.
%
%  [IDX,T] = dePDDP(X,K) returns the cluster assignment, (IDX), and the  binary 
%  tree (T) containing the cluster hierarchy. If K==[] the number of clusters is estimated
%
%  [IDX,T] = dePDDP(X, K, 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional parameters
%  in the form of Name,Value pairs. 
%
%  'bandwidth' - Bandwidth parameter
%	Function Handle: bandwidth(X,pars) returns bandwidth (positive scalar)
%	(default: bandwidth =0.9* std(projections) * N^(-0.2))
%
%  'split_index' - Criterion determining which cluster to split	next (only relevant if K is specified)
%	Function Handle: index = split_index(v, X, pars)
%			(v: projection vector, X:data matrix, pars: parameters structure)
%	Cluster with MAXIMUM INDEX is split at each step of the algorithm
%	Two standard choices of split index can be enabled by setting 'split_index' to 
%	one of the strings below:
%		+ 'size':    Split largest cluster
%		+ 'fval':    Split cluster whose local minimum on the 1D KDE is the lowest
%	(default: split_index = 'fval')
%
%  'minsize' - Minimum cluster size (integer)
%	(default minsize = 1)
%
%  'labels' - true cluster labels. Only used for performance assessment.
%
%Reference:
%S.K. Tasoulis, D.K. Tasoulis and V.P. Plagianakos. Enhancing principal direction divisive clustering.
%Pattern Recognition, 43(10):3391-3411, 2010.

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if nargin==1,
	K = size(X,1);
end

% Set default parameters
pars = struct();
pars.bandwidth = @(x,p)(0.9*size(x,1)^(-0.2)*std(x* pcacomp(x,1)));
pars.minsize = 1;
pars.split_index = 'fval';
pars.labels = [];

pars = myparser(X,K,varargin,pars);
[N, dim] = size(X);

%%%%%%%%% HIERARCHICAL ALGORITHM EFFECTIVELY STARTS HERE

% pass{i}: contains binary allocation of observations allocated to i-th tree node based on separating hyperplane 
pass = {};
% location of each node in t
node_index = [1];
% split_index[i]: contains the value of the split index for the tree node with number node_index[i]
split_index = [];
[opthp, pass{1}, split_index(1)] = gpp(X, @depddppp, pars);
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
t = ctree(opthp, struct('data',X,'method','dePDDP','params',pars));

if isempty(K),
	K = N;
end
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
		list{j} = list{id}(pass{id} == 2*i-3);

		% Partition data subset
		[opthp, pass{j}, split_index(j)] = gpp(X(list{j},:), @depddppp, pars);
		% if PP method fails to identify valid separating hyperplane
		if isempty(opthp),
			% Use projection matrix from parent to enable visualisation
			opthp = gsep(t.Node{node_index(id)}.v, -inf, [1:length(list{j})]', @depddppp, pars);
		end
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
