% Class implementing cluster hierarchy in tree data structure
%
% This class implements a tree data structure used to store cluster
% hierarchies. Nodes are mainly accessed through their index. The index of a
% node is returned when it is added to the ctree, and actually corresponds to
% the order of addition. 
% 
%  Based on MATLAB tree-class: Jean-Yves Tinevez <tinevez@pasteur.fr> March 2012
%  Modified by: Nicos Pavlidis 2018

classdef ctree
properties (Access = public)
	% Stores cluster separator at each node
	Node = { [] };

	% Index of the parent node. The root of the ctree has a parent index equal to 0.
	Parent = [ 0 ];

	% Name of clustering algorithm
	method = '';

	% Vector of cluster assignments
	cluster =  [];

	% N-by-D Data matrix used for clustering after preprocessing
	data = [];

	% Parameters of clustering algorithm
	params = struct();
end
    
methods (Access = public)
	function obj = ctree(content, in)
	%Constructor
	%[OBJ, ROOT_ID] = CTREE(CONTENT, IN)
	% t = CTREE(ctree) is the copy constructor for this class
	% It returns a new ctree where the node order and content
	% is duplicated from the ctree argument.
	% 
	% t = CTREE(root_content, in) where 'root_content' is not a ctree,
	% initialises a new ctree with only the root node, and set its
	% content to be 'root_content'. A structure containing parameters 
	% 'in' is also added
	   
		if nargin < 1
			return;
		end
	    
		if isa(content, 'ctree')
			% Copy constructor
			obj.Node = content.Node;
			obj.Parent = content.Parent;

			obj.cluster = content.cluster;
			obj.data = content.data;
			obj.method = content.method;
			obj.params = content.params;
		else
			% New object with only root content
			obj.Node = { content };
			root_ID = 1;
		
			if nargin>1,
				obj.data = in.data;
				n = size(obj.data,1);

				if isfield(in, 'cluster') & length(in.cluster) == n,
					obj.cluster = in.cluster;
				else
					obj.cluster = ones(n,1);
				end

				obj.method = in.method;
				obj.params = in.params;
			end
		end
		% Add default parameters
		if ~isfield(obj.params, 'minsize'),
			obj.params.minsize = 1;
		end
		if ~isfield(obj.params, 'maxit'),
			obj.params.maxit = 50;
		end
		if ~isfield(obj.params, 'ftol'),
			obj.params.maxiter = 1.e-5;
		end
		if ~isfield(obj.params, 'verb'),
			obj.params.maxiter = 0;
		end
	end % CONSTRUCTOR
	
	
	% METHODS
	function [obj, ID] = addnode(obj, parent, data)
	%ADDNODE attach a new node to a parent node
	%[OBJ, ID] = ADDNODE(OBJ, PARENT, DATA)
	% 
	% Creates new node with content 'data', and attach it as a child of the
	% node with index 'parent_index'. Return the modified ctree.
	%
	% [ctree, ID] = ctree.ADDNODE(...) returns the modified ctree and
	% the index of the newly created node.
	    
	    if parent < 0 || parent > numel(obj.Parent)
		error('MATLAB:ctree:addnode', ...
		    'Cannot add to unknown parent with index %d.\n', parent)
	    end
	    
	    if parent == 0,
		    % Replace the whole ctree by overiding the root.
		    obj.Node = { data };
		    obj.Parent = 0;
		    ID = 1;
		    return;
	    end
	    % Expand the cell by
	    obj.Node{end+1} = data;
	    obj.Parent = [obj.Parent; parent];
	    ID = numel(obj.Node);
	end % ADDNODE
	
	function flag = isleaf(obj, ID)
	%ISLEAF True if given Node with index ID is leaf node.
	%FLAG = ISLEAF(OBJ, ID)
		if ID < 1 | ID > numel(obj.Parent)
			error('MATLAB:ctree:isleaf', 'No node with ID %d.', ID)
		end
		flag = ~any( obj.Parent == ID );
	end % ISLEAF
	
	function IDs = findleaves(obj)
	%FINDLEAVES Return the indices of all leaves 
	%IDS = FINDLEAVES(OBJ)
		parents = obj.Parent;
		IDs = (1 : numel(parents)); % All IDs
		IDs = setdiff(IDs, parents); % Remove those which are marked as parent
	end
	
	function content = get(obj, ID)
	%GET Return a copy Node with index ID.
	%CONTENT = GET(OBJ, ID)
		content = obj.Node{ID};
	end

	function obj = set(obj, ID, content)
	%SET Set the content of node with index ID and return the modifed ctree.
	%OBJ = SET(OBJ, ID, CONTENT)
		obj.Node{ID} = content;
	end

	function obj = setfield(obj, ID, field, value)
	%SETFIELD Set the value FIELD of the Node with index ID to VALUE
	%OBJ = SETFIELD(OBJ, ID, FIELD, VALUE)
		if isfield(obj.Node{ID}, field),
			obj.Node{ID} = setfield(obj.Node{ID},field,value);
		else
			error('MATLAB:ctree:setfield', 'No field with name %s.', field)
		end
	end
	
	function IDs = getchildren(obj, ID)
	%GETCHILDREN Return column vector of node indices (IDs) of children of Node ID.
	%IDS = GETCHILDREN(OBJ, ID)
		parent = obj.Parent;
		IDs = find( parent == ID );
		IDs = IDs';
	end
	
	function ID = getparent(obj, ID)
	%GETPARENT  Return the ID of the parent of the given node.
	%ID = GETPARENT(OBJ, ID)
		if ID < 1 | ID > numel(obj.Parent)
			error('MATLAB:ctree:getparent', 'No node with ID %d.', ID)
		end
		ID = obj.Parent(ID);
	end
	
	function IDs = getsiblings(obj, ID)
	%GETSIBLINGS  Return column vector of indices (IDs) of siblings of node with index (ID) (including node ID itself)
	%IDS = GETSIBLINGS(OBJ, ID)

		if ID < 1 | ID > numel(obj.Parent)
			error('MATLAB:ctree:getsiblings', 'No node with ID %d.', ID)
		end
		if ID == 1 % Special case: the root
			IDs = 1;
			return;
		end
		IDs = obj.getchildren(obj.getparent(ID));
	end
	
	function obj = chop(obj, node)
	%CHOP Remove node with index (NODE) and all subnodes from the given ctree
	%OBJ = CHOP(OBJ, NODE)
		iterator = obj.depthfirstiterator(node);
		% Build new parent array
		np = obj.Parent;
		% Remove unwanted nodes
		 np(iterator) = [];

		% Shift parent value: if a parent were after some nodes we removed, we
		% need to shift its value by an amount equal to the number of parent we
		% removed, and that were BEFORE the target parent
		for i = 1:numel(np),
			np(i) = np(i) - sum(np(i) > iterator);
		end
		
		obj.Parent = np;
		obj.Node(iterator) = [];
	end

	function obj = prune(obj, ID)
	%PRUNE Remove everything below node with index (ID)
	%OBJ = PRUNE(OBJ, ID)
		if ID < 1 || ID > numel(obj.Parent),
			error('MATLAB:ctree:prune ', 'No node with ID %d.', ID)
		end
		if obj.isleaf(ID),
			error('MATLAB:ctree:prune ', 'Node with ID %d is a leaf.', ID)
		end
		% nodes to remove
		rm = sort(obj.getchildren(ID));
		obj.Node{ID}.tree_params.leaf = 1;
		for i = length(rm):-1:1,
			obj = obj.chop(rm(i));
		end
		iterator = obj.depthfirstiterator;
		for i = iterator,
			obj.Node{i}.tree_params.nodeid = i;
		end
		obj.cluster = tree2clusters(obj);
	end

	function n = nnodes(obj)
	%NNODES Number of nodes in the ctree. 
	%N = NNODES(OBJ)
		n = numel(obj.Parent);
	end

	function disp(obj)
	%DISP Print tree on terminal
	%DISP(OBJ)
		disp(obj.tostring);
	end

	function d = depth(obj)
	%DEPTH Depth of tree
	%D = DEPTH(OBJ)
		dt = obj.depthtree;
		it = dt.depthfirstiterator;
		d = 0;
		for i = it,
			d = max(d, dt.get(i));
		end
	end

	function IDs = conditioniterator(obj, startNode, condition, sorted)
	%CONDITIONITERATOR returns (IDs) of nodes satisfying CONDITION
	%IDS = CONDITIONITERATOR(OBJ, STARTNODE, CONDITION, SORTED)
		if nargin < 2, 
			startNode = 1;
		elseif isempty(startNode),
			startNode = 1;
		end

		if nargin < 3, 
			condition = @(x) true;
		elseif isempty(condition),
			condition = @(x) true;
		end

		if nargin < 4,
			sorted = false;
		elseif isempty(sorted),
			sorted = false;
		end

		IDs = recurse(startNode);
		function val = recurse(node)
			val = node;
			content = obj.Node{val};
			if obj.isleaf(node) | ~condition(content)
				return;
			else
				children = obj.getchildren(node);
				if sorted && numel(children) > 1
					contents = obj.Node(children);
					[~, sorting_array] = sortrows(contents);
					children = children(sorting_array);
				end
				cellval = cell(numel(children), 1);
				for i = 1 : numel(children)
					content = obj.Node{children(i)};
					if ~condition(content)
						continue
					end
					cellval{i} = recurse(children(i));
				end
				val = [ val cellval{:} ] ;
			end
		end
	end

	function IDs = depthfirstiterator(obj, startNode, sorted)
	%DEPTHFIRSTITERATOR  Index sequence traversing the tree, depth first.
	%IDS = DEPTHFIRSTITERATOR(OBJ, STARTNODE, SORTED)
		if nargin < 2, startNode = 1; end
		if nargin < 3, sorted = false; end
		IDs = recurse(startNode);
		function val = recurse(node)
			val = node;
			if obj.isleaf(node),
				return;
			else
				children = obj.getchildren(node);
				if sorted & numel(children) > 1,
					contents = obj.Node(children);
					[~, sorting_array] = sortrows(contents);
					children = children(sorting_array);
				end
				cellval = cell(numel(children), 1);
				for i = 1:numel(children),
					cellval{i} = recurse(children(i));
				end
				val = [val cellval{:}];
			end
		end
	end

	function dt = depthtree(obj)
	%DEPTHTREE Create a coordinated tree where each node holds its depth.
	%DT = DEPTHTREE(OBJ)
		dt = ctree(obj);
		dt = dt.set(1, 0);
		iterator = obj.depthfirstiterator;
		iterator(1) = []; % Remove root

		for i = iterator,
			parent = dt.Parent(i);
			parentDepth = dt.get(parent);
			dt = dt.set(i, parentDepth+1);
		end
	end

	function obj = perf(obj,labels)
	%Evaluate performance at each node of ctree with respect to Purity and Success Ratio
	%OBJ = PERF(OBJ,LABELS)
	%
	% Input:
	%	(labels): True cluster assignment

		assert(size(obj.data,1) == length(labels));
		labels = fixLabels(labels);
		nc = max(labels);

		% this should be a separate function
		if nc > 1,
			for i=1:nnodes(obj),
				l = labels( abs(obj.Node{i}.idx) );
				obj.Node{i}.tree_params.purity = max(mycrosstab(l, ones(length(l),1)),[],1)/length(l);
				obj.Node{i}.tree_params.SR = success_ratio(0.5*sign(obj.Node{i}.idx) + 1.5, l);
				obj.Node{i}.tree_params.leaf = obj.isleaf(i);
			end
		end
	end

	function nplot(obj, node_id, labels, colours)
	%NPLOT Visualisation of binary partition at specific node in ctree
	%NPLOT(OBJ, NODE_ID, LABELS, COLOURS)
	% Input:
	%	(tree): Divisive hierarchical clustering model obtained from PPCI toolbox
	%	(node_id): Node number
	%	(labels): (Optional) True cluster assignment (used to colour data)
	%	(colours): (Optional) Colours for different clusters (row-wise) (active when 'labels' specified)

		if nargin < 2,
			error('Node number needs to be specified');
		elseif numel(node_id) > 1 | node_id<1,
			error('Node number needs to be positive integer');
		elseif node_id > obj.nnodes,
			error('Node number exceeds tree size');
		end

		%% Fix labels
		defLabels = true;
		if nargin<3,
			labels = [];
			defLabels = false;
		elseif isempty(labels),
			labels = [];
			defLabels = false;
		else
			assert(length(labels)==size(obj.data,1),'Labels vector incompatible with data matrix');
			labels = fixLabels(labels);
		end

		if nargin < 4, colours = []; end

		hFig = figure(1);
		clf;
		set(hFig, 'Position', [0 0 1024 1024]);

		if exist ('OCTAVE_VERSION', 'builtin') > 0,
			plot(obj.Node{node_id}, obj.data, labels, colours);
			return;
		end

		% Plot 2-dimensional data projections and associated densities
		a = subplot(3,4, [1 2 3 5 6 7 9 10 11]);
		a = plot(obj.Node{node_id}, obj.data, labels, colours, a);
	 
		% Plot position of node in cluster cluster hierarchy (tree) 
		% and report summary information
		subplot(3,4, [4 8]);
		hold on;
		width = 1./ length(obj.findleaves);
		height = 1/(obj.depth+0);
		add_subtree(obj,1,0,1,1,height,width,node_id)
		hold off;
		axis([0,1,0,1]);
		set(gca,'Xtick',[],'Ytick',[],'XColor','w','YColor','w')

		subplot(3,4,12);
		axis([0,1,0,1]);

		if defLabels,
			index = abs(obj.Node{node_id}.idx);
			perf = sprintf('purity: %1.3f\nsuccess ratio: %1.3f', ...
				purity(sign(obj.Node{node_id}.idx), labels(index)), ...
				success_ratio(sign(obj.Node{node_id}.idx), labels(index)));
		else 
			perf = '';
		end

		descr = {sprintf('node id: %i',node_id), sprintf('N: %i',length(obj.Node{node_id}.idx))};

		% relevant only for Minimum Density Hyperplanes
		if isfield(obj.Node{node_id}, 'rel_dep'),
			descr = {descr{:}, sprintf('Rel. Depth: %1.3f', obj.Node{node_id}.rel_dep)};
		end
		descr = {descr{:}, perf};

		if obj.isleaf(node_id) & isinf(obj.Node{node_id}.fval),
			descr = {descr{:}, 'no valid separator identified'};
		end

		text(0.025, 1.0, descr);
		axis off;

		function add_subtree(Tr,cur_node,L,U,y,height,width, node)
			
			if ~Tr.isleaf(cur_node),
				% the width of subtrees = # of leaves
				lr = Tr.getchildren(cur_node);
				w1 = length( findleaves( Tr.subtree(lr(1)) ) );
				w2 = length( findleaves( Tr.subtree(lr(2)) ) );
				M = L+ (U-L)*w1 /(w1+w2);

				% plot the line segments first from the bottom of the cluster and then left and right
				%text(0.5*(L+U),y - 0.1*height, sprintf('%i',cur_node));
				if cur_node==1,
					plot( [0.5*(L+U), 0.5*(L+U)], [y,y-0.20*height],'k-','LineWidth',1);
					dec = 0.2;
				else
					dec = 0;
				end

				plot( [0.5*(L+M), 0.5*(M+U)], [y-dec*height, y-dec*height],'k-','LineWidth',1);
				plot( [0.5*(L+M), 0.5*(M+L)], [y-dec*height,y-1.00*height],'k-','LineWidth',1);
				plot( [0.5*(U+M), 0.5*(M+U)], [y-dec*height,y-1.00*height],'k-','LineWidth',1);

				add_subtree(Tr,lr(1),L,M,y-1.00*height,height,(M-L)/w1, node);
				add_subtree(Tr,lr(2),M,U,y-1.00*height,height,(U-M)/w2,node);
			end

			if node == cur_node, 
				plot( 0.5*(L+U), y,'o','MarkerSize',8,'MarkerEdgeColor','g', 'MarkerFaceColor','g'); 
			else
				plot( 0.5*(L+U), y,'o','MarkerSize',8,'MarkerEdgeColor','w', 'MarkerFaceColor','w'); 
			end
			text(0.5*(L+U), y, num2str(cur_node),'HorizontalAlignment','center', 'FontSize',8);
		end
	end

	function plot(obj, labels, colours)
	%PLOT Plot divisive hierarchical clustering model
	%PLOT(OBJ, LABELS, COLOURS)
	%
	% Inputs:
	%	(labels): True cluster assignment (optional)
	%	(colours): Colours for different clusters (row-wise) (optional; active when 'labels' specified)

		if nargin < 3, 
			colours = []; 
		end

		if nargin<2, 
			labels = ones(size(obj.data,1),1);
		elseif isempty(labels),
			labels = ones(size(obj.data,1),1);
		else 
			labels = fixLabels(labels);
		end

		% Set colours used projections
		colours = palette(max(labels),colours);

		hFig = figure(2);
		clf;
		set(hFig, 'Position', [0 0 1024 1024]);
		axis off;
		axis([0,1,0,1]);

		width = 1./ length(obj.findleaves);
		height = 0.8/(obj.depth+1);

		function add_subtree(obj,labels,printID,nodeid,L,U,y,height,width,colours)
			index = abs(obj.Node{nodeid}.idx);
			X = bsxfun(@minus,obj.data(index,:), mean(obj.data(index,:)) );
			[v,w] = obj.Node{nodeid}.tree_plot_data(X);
			X = X*[v,w];
			% normalise in [0,1]
			X = bsxfun(@rdivide, bsxfun(@minus,X,min(X)), max(X)-min(X));

			X(:,1) = width*X(:,1) + 0.5*(L+U-width);
			X(:,2) = height*X(:,2) + y-height;

			if max(labels)>1,
				l = unique(labels(index));
				for i=1:length(l),
					plot(X(labels(index)==l(i),1), X(labels(index)==l(i),2),'o', ...
						'MarkerSize',3.0,'Color',colours(l(i),:));
				end
			else
				if ~obj.isleaf(nodeid),
					left = find(obj.Node{nodeid}.idx < 0);
					if ~isempty(left), 
						plot(X(left,1), X(left,2), 'o', 'MarkerSize',3.0,...
							'Color', colours(1,:));
					end
				
					left = find(obj.Node{nodeid}.idx > 0);
					if ~isempty(left), 
						plot(X(left,1), X(left,2), 'o', 'MarkerSize',3.0,...
							'Color', colours(2,:));
					end
				else
					% if leaf then colour of observations depends on partition at parent node
					cl = obj.getchildren(obj.getparent(nodeid));
					o = ifelse( nodeid==cl(1), colours(1,:), colours(2,:));
					plot(X(:,1), X(:,2), 'o', 'MarkerSize',3.0,'Color', o);
				end
			end

			if printID, 
				text(0.5*(L+U), y-1.0*height, sprintf('%i',nodeid)); 
			end

			if ~obj.isleaf(nodeid),
				% the width of subtrees = # of leaves
				lr = obj.getchildren(nodeid);
				w1 = length( findleaves( obj.subtree(lr(1)) ) );
				w2 = length( findleaves( obj.subtree(lr(2)) ) );
				M = L+ (U-L)*w1 /(w1+w2);

				% plot the line segments first from the bottom of the cluster and then left and right
				plot( [0.5*(L+U), 0.5*(L+U)], [y-height,y-1.125*height],'k-','LineWidth',1);
				plot( [0.5*(L+M), 0.5*(M+U)], [y-1.125*height,y-1.125*height],'k-','LineWidth',1);
				plot( [0.5*(L+M), 0.5*(M+L)], [y-1.125*height,y-1.25*height],'k-','LineWidth',1);
				plot( [0.5*(U+M), 0.5*(M+U)], [y-1.125*height,y-1.25*height],'k-','LineWidth',1);

				add_subtree(obj,labels,printID,lr(1),L,M,y-1.25*height,height,(M-L)/w1,colours);
				add_subtree(obj,labels,printID,lr(2),M,U,y-1.25*height,height,(U-M)/w2,colours);
			end
		end

		hold on;
		add_subtree(obj,labels,1,1,0,1,1,height,width,colours);
		hold off;
	end

	function obj = split(obj, id, varargin)
	%SPLIT Partitions leaf node with number (id)
	%OBJN = SPLIT(OBJ, ID, VARARGIN)
	%	Internal nodes of (obj) cannot be modified through this function. 
	%
	% If no further arguments are specified all the arguments used to construct (OBJ)
	% in the first instance are used to bi-partition the leaf node (ID). The
	% user can modify this behaviour by setting any argument of the clustering
	% algorithm through the same (Name,Value) syntax used in creating the cluster
	% hierarchy originally.
	%
	% Returns ctree object (OBJN)

		if ~obj.isleaf(id),
			error('MATLAB:ctree:split', 'Node with ID %i is not a leaf', id);
			return;
		end

		% Copy parameter settings (should this be at the tree or the node level?)
		pars = obj.params;

		% Number of observations allocated to present node
		N = length(obj.Node{id}.idx);

		% These are used to enable the user to visualise the partition of a specific node
		nc = 1;
		labels = [];
		colours = [];

		% If any of the below conditions hold then either the current node (cluster) cannot
		% be split or the existing split needs to be revised
		if isempty(obj.Node{id}) | isinf(obj.Node{id}.fval) | (N < 2*pars.minsize) | ~isempty(varargin),

			% If node (cluster) was deemed to not admit a valid separating hyperplane
			% attempt to partition it ONLY if optional parameter values have been provided
			if isempty(varargin),
				if N < 2*pars.minsize,
					error('MATLAB:ctree:split', ...
						'Node to be split contains %i observations < 2 * (minsize= %i)',N, pars.minsize);
				else
					error('MATLAB:ctree:split', ...
						'No valid hyperplane for this leaf: Set optional arguments to split node');
				end
			end

			% Validate & Update Arguments
			if mod(length(varargin),2)~=0,
				error('MATLAB:ctree:split', 'Optional arguments specified incorrectly');
			end

			% Update parameter settings (and set labels)
			for i = 1:2:(length(varargin)-1),

				if strcmp(varargin{i},'labels'),
					if length(varargin{i+1}) ~= size(obj.data,1),
						error('MATLAB:ctree:split','Label vector of incorrect length');
					end
					labels = fixLabels(varargin{i+1});

				elseif strcmp(varargin{i},'colours'),
					colours = varargin{i+1};	

				% if field does not exist
				elseif ~isfield(pars, varargin{i}) 
					error('MATLAB:ctree:split', ...
						'Error: "%s" not valid option for "%s"',varargin{i}, obj.method);
				else
					% Update parameter settings
					pars = setfield(pars, varargin{i}, varargin{i+1});
				end
			end

			if N < 2*pars.minsize,
				% Ensure minsize constraint is not violated
				error('MATLAB:ctree:split', ...
					'New minsize %i does not allow cluster of size %i to be split', ...
					pars.minsize, length(obj.Node{id}.idx));
			end

			index = abs(obj.Node{id}.idx);
			% Set labels and colours according to data subset
			if ~isempty(labels),
				la = labels(index);
				colours = palette(max(labels), colours);
				cl = colours(unique(la), :);
				% if labels have been specified plot
				if max(labels)>1,
					pars.verb = 1;
				end
			else
				la = ones(length(index),1);
				cl = palette(1, colours);
			end

			% Split node ID 
			[node, idx, spindex] = obj.Node{id}.split(obj.data(index,:), pars, la, cl);

			% Verify that a valid HP separator has been obtained
			if isinf(node.fval),
				error('MATLAB:ctree:split', 'No valid binary split identified for this cluster');
			end
			n1 = sum(idx==1);
			if min(n1, N-n1) < pars.minsize,
				% This should never occur if the PP algorithm incorporates minsize constraint
				error('MATLAB:ctree:split', 'Cluster split violates minimum size constraint');
			end

			node.tree_params.nodeid = obj.Node{id}.tree_params.nodeid;
			% Set index of observations allocated to node
			node.idx = index .* idx;

			% Overwrite previous node in tree
			obj.Node{id} = node;
		end

		% Leaf node is partitioned through valid separator
		obj.Node{id}.tree_params.leaf = 0;

		% To avoid plotting partitions of leaf nodes
		pars.verb = 0;
		for i = [-1,1],
			index = abs( obj.Node{id}.idx(sign(obj.Node{id}.idx) == i) );
			[node, idx, spindex] = obj.Node{id}.split(obj.data(index,:), pars);

			node.tree_params.nodeid = nnodes(obj) + 1;
			node.tree_params.leaf = 1;
			% Set index of observations allocated to node
			node.idx = index .* idx;

			% extend the tree by splitting the node with number: node_index(id)
			obj = obj.addnode(id, node);
		end
		obj.cluster = tree2clusters(obj);
	end %SPLIT

	function [st, index] = subtree(obj, node, condition)
	%SUBTREE Returns the sub-tree made of all the nodes below the given one
	%[ST, INDEX] = SUBTREE(OBJ, NODE, CONDITION)

		if nargin < 3
			condition = @(s) true;
		end
		% Get indices of the subtree.
		iterator = obj.conditioniterator(node, condition);
		% Copy the content of the tree related to the subtree
		parents = obj.Parent(iterator);
		nodes = obj.Node(iterator);
		% Revamp parent indices
		newParents = NaN(numel(parents), 1);
		newParents(1) = 0; % The new root node
		for i = 2 : numel(parents)
			pr = parents(i);
			newParents(i) = find( iterator == pr, 1, 'first');
		end
		% Create a new tree with the sub-content
		st = ctree;
		st.Node = nodes;
		st.Parent = newParents;
		% Return the link new node index -> old node index
		index = ctree; %iterator;
		index.Parent = newParents;
		index.Node = num2cell(iterator);
	end % SUBTREE


	function str = tostring(obj, sorted)
	%TOSTRING Returns char matrix reprensenting the ctree
	%STR = TOSTRING(OBJ, SORTED)
		if nargin < 2
			sorted = false;
		end
		%% Use a tree to compute spaces;
		% 1. Generate string representation of content as a tree
		strContentTree = treefun(obj, @contentToString);
	    
		% 1.25 Compute space requirements
		spaceTree = treefun(strContentTree, @numel);
	    
		% 1.5 Each children must at least have a size determined by its parent
		iterator = spaceTree.depthfirstiterator;
		for i = iterator
			parent = spaceTree.getparent(i);
			if parent == 0,
				continue
			end

			nSiblings = numel(spaceTree.getsiblings(i));
			parentWidth = spaceTree.get(parent);
			minWidth = ceil(parentWidth / nSiblings);
			thisWidth = spaceTree.get(i);
			if thisWidth < minWidth
				spaceTree = spaceTree.set(i, minWidth);
			end
		end
	    
		% 2. Add 2 for proper spacing
		spaceTree = treefun2(spaceTree, 2, @plus);
		%spaceTree = spaceTree + 2;

		% 3. Build cumulative space tree
		spaceTree = recursivecumfun(spaceTree, @sum);

		% Put at least 1 when there is nothing
		iterator = spaceTree.depthfirstiterator;
		for i = iterator
			if isempty(spaceTree.get(i)) | spaceTree.get(i) == 0,
				   spaceTree = spaceTree.set(i, 1);
			end
		end

		%% Iterate in the tree a first time to put the node content and the vertical bars
		nLevels = obj.depth + 1;
		depth = obj.depthtree;
		iterator = obj.depthfirstiterator(1, sorted);
		str = repmat(' ', 3 * nLevels - 2, spaceTree.get(1));

		% The last column in which something was writen in at given level.
		columns = zeros(nLevels, 1);
		% And a matching tree to store at what position we write. We will need
		% this in the second iteration.
		columnTree = ctree(obj);

		for i = 1 : numel(iterator)

			ID = iterator(i);
			level = depth.get(ID) + 1; % Because the 2 tree have the same structure
			spaceWidth = spaceTree.get(ID);
			ds = ceil( spaceWidth / 2 );

			index = columns(level) + ds;
			row = 3 * (level-1) + 1;

			if level > 1
				% Line 0: the '+'
				if numel(obj.getsiblings(ID)) > 1
					ch = '+';
				else
					ch = '|';
				end
				str(row-2, index) = ch;
				% Line 1: the vertical tick
				str(row-1, index) = '|';
			end
			% Line 2: the content
			contentStr = strContentTree.get(ID);
			contentWidth = numel(contentStr);
			dc = floor(contentWidth / 2);
			str(row, index-dc : index + contentWidth-1-dc) = contentStr;
			columnTree = columnTree.set(ID, columns(level));
			columns(level) = columns(level) + spaceWidth;

			% Is a leaf? Then we move the cursor to the right for the next
			% levels. If we do not do it, the display is going to be messy as
			% soon as the nodes do not have the same number of children.

			if obj.isleaf(ID) & level < nLevels
				for j = level+1 : nLevels
					columns(j) = columns(j) + spaceWidth;
				end
			end

		end

		%% A second iteration to draw horizontal bars
		for i = 1 : numel(iterator)
			ID = iterator(i);
			level = depth.get(ID) + 1; % Because the 2 tree have the same structure
			 if level == nLevels | obj.isleaf(ID)
				continue
			end
			spaceWidth = spaceTree.get(ID);
			ds = floor( spaceWidth / 2 );
			col = columnTree.get(ID);
			row = 3 * (level-1) - 1;

			% Move to the level below, and edit the line with the '+' so that
			% they show siblings
			index = col + ds;
			if numel(obj.getchildren(ID)) > 1
				ch = '+';
			else
				ch = ' ';
			end
			str(row+3, index) = ch;
			childbar = str(row+3, col + 1 : col + spaceWidth);
			ai = find(childbar == '+' | childbar == '|');

			if isempty(ai)
			    continue
			end
			fi = ai(1);
			li = ai(end);

			toReplace = setdiff( fi:li, ai);
			childbar(toReplace) = '-';
			str(row+3, col + 1 : col + spaceWidth) = childbar;
		end
	end %TOSTRING
end % METHODS

methods (Access = private)
	function cumtree = recursivecumfun(obj, fun)
	%RECURSIVECUMFUN  Create a tree where content is calculated recursively from children
	%CUMTREE = RECURSIVECUMFUN(OBJ, FUN)

		% Prepare a blank tree for holding values
		cumtree = ctree(obj);
		descend(1);
		function val = descend(n) % current node index
			if cumtree.isleaf(n)
				val = cumtree.get(n);
			else
				children = cumtree.getchildren(n);
				nChildren = numel(children);
				stock = NaN(nChildren, 1);

				for i = 1:nChildren,
					child = children(i);
					stock(i) = descend(child);
				end
				val = fun(stock);
				% Store value in new tree
				cumtree.Node(n) = {val};
			end
		end
	end % RECURSIVECUMFUN

	
	function newTree = treefun(obj1, fun)
	%TREEFUN  Create a new tree by applying function handle to each node of the source tree.
	%NEWTREE = TREEFUN(OBJ, FUN)
	
	    if ~isa(fun, 'function_handle')
		error('MATLAB:tree:treefun', ...
		    'Second argument, fun, must be a function handle. Got a %s.', ...
		    class(fun));
	    end
	    % First we copy
	    newTree = ctree(obj1, 'clear'); 
	    % Then we override
	    newTree.Node = cellfun(fun, obj1.Node, 'UniformOutput', false);
	end % TREEFUN

	function obj1 = treefun2(obj1, val, fun)
	%TREEFUN2  Two-arguments function on trees, with scalar expansion
	%OBJ = TREEFUN2(OBJ, VAL, FUN)
		if ~isa(obj1, 'ctree')
			tmp = obj1;
			obj = val;
			val = tmp;
		end

		if isa(val, 'ctree')
			content = cellfun(fun, obj1.Node, val.Node, 'UniformOutput', false);
		else
			content = cellfun(@(x) fun(x, val), obj1.Node,'UniformOutput', false);
		end
		obj1.Node = content;
	end % TREEFUN2
end % PRIVATE METHODS METHODS

end

