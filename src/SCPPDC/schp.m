%Class implementing a linear projection subspace of arbitrary dimensions estimated through SCPP
classdef schp

properties (SetAccess = public)
	% matrix of unit-length vector(s) defining subspace
	v;
	% projection index
	fval;
	% Parameters for separating hyperplane
	params;
	% last two attributes are relevant for HPs in cluster hierarchies

	% row-numbers of observations assigned to specific cluster
	% (sign indicates allocation to left (-) or right (+) child)
	idx;

	% Parameters employed by hierarchical divisive clustering algorithms
	tree_params;
end

methods
	function obj = schp(v1,f1,p1,id1,t1)
	% CONSTRUCTOR

		if nargin < 1,
			obj.v = [];
			obj.fval = [];
			obj.params = schp.def_opar();
			obj.idx = [];
			obj.tree_params = [];
			return;
		end

		if nargin==1 & isa(v1, 'schp'),
			% Copy constructor
			obj.v = bsxfun(@rdivide, v1.v, sqrt(sum(v1.v.^2)));
			obj.fval = v1.fval;
			obj.params = schp.def_opar(v1.params);
			obj.idx = v1.idx;
			obj.tree_params = v1.tree_params;
			return;
		end

		if nargin < 5, t1 = []; end;

		if nargin < 4, 
			id1 = []; 
		elseif ~isempty(id1),
			assert(numel(id1)>1 & size(id1,1)>1 & size(id1,2)==1 & ~sum(id1~=floor(id1)),'idx: must be a column vector of integers');
		end;

		if nargin < 3,
			p1 = [];
		end;

		if nargin < 2, 
			f1 = [];
		elseif ~isempty(f1),
			assert(numel(f1)==1 & isreal(f1),'f: must be a real number');
		end

		assert(size(v1,1)>1 & isempty(find(v1==nan)) & isempty(find(v1==inf)), 'v: must be a column vector');
		obj.v = bsxfun(@rdivide, v1, sqrt(sum(v1.^2)));

		obj.fval = f1;
		obj.params = p1;

		obj.idx = id1;

		% Common parameters for all Projection Pursuit algorithms
		obj.params = schp.def_opar(p1);
		
		% Parameters employed by tree algorithms
		obj.tree_params = t1;
		if ~isfield(obj.tree_params, 'nodeid'),
			obj.tree_params.nodeid = [];
		end
		if ~isfield(obj.tree_params, 'leaf'),
			obj.tree_params.leaf = 0;
		end
		if ~isfield(obj.tree_params, 'SR'),
			obj.tree_params.SR = [];
		end
		if ~isfield(obj.tree_params, 'purity'),
			obj.tree_params.purity = [];
		end
	end


	function out = isempty(obj)
	%ISEMPTY check for initialisation
		out = isempty(obj.v);
	end

	function idx = cluster(obj,data)
	%CLUSTER assign observations to element of binary partition
		if isempty(obj.v),
			error('MATLAB:schp:cluster','Cannot cluster with non-initialised SCHP');
		end
		if size(data,2) ~= size(obj.v,1),	
			error('MATLAB:schp:cluster','Data matrix - HP Incompatible dimensions');
		end
		% Generic Spectral Clustering algorithm
		idx = scppNJW(2,obj.v,data,obj.params.params.sigma,[],obj.params.beta, obj.params.delta);
	end

	function [v,w] = tree_plot_data(obj,X)
	%TREE_PLOT_DATA Two-dimensional space for visualisation
		if isempty(obj.v),
			error('Cannot plot SCHP with uninitialised (v)');
		end
		v = bsxfun(@rdivide, obj.v, sqrt(sum(obj.v.^2)));

		[d,pdim] = size(v);
		if pdim ==1,
			w = pcacomp(X - (X*v)*v', 1);
			if isempty(w),
				v = [1;zeros(d-1,1)];
				w = [0;1;zeros(d-2,1)];
			end
		else
			[~,id] = sort(std(X*v), 'descend');
			w = v(:,id(2));
			v = v(:,id(1));
		end

	end

	function hFig = plot(obj,data,labels,colours,hFig,clusters)
	%PLOT Visualisation of binary partition
		if isempty(obj.v),
			error('Cannot plot SCHP with non-initialised projection matrix (v)');
		end

		if nargin<2,
			error('Data matrix undefined');
		end

		if size(data,2) ~= size(obj.v,1),
			error('Data matrix is incompatible with projection vector');

		end

		if nargin < 5, hFig = []; end;

		% if labels are undefined define all observations to belong to the same cluster
		if nargin<3,
			labels = ones(size(data,1),1);
			nc = 1;
			defLabels = false;
		elseif isempty(labels),
			labels = ones(size(data,1),1);
			nc = 1;
			defLabels = false;
		else
			assert(length(labels)==size(data,1), ...
				'Labels vector is incompatible with data matrix');
			% Ensure labels are consecutive integers {1,...,K}
			labels = fixLabels(labels);
			% number of actual clusters
			nc = max(labels);
			defLabels = true;
		end
		% Set colours
		if nargin < 4,
			colours = [];
		end;
		colours = palette(nc,colours);

		% Identify subset of data to be plotted
		if isempty(obj.idx),
			idx = [1:size(data,1)]';
			if nargin < 6,
				clusters = [];
			end
		else
			clusters = 1+(obj.idx>0);
			idx = abs(obj.idx);
		end
		data = data(idx,:);
		labels = labels(idx);

		% estimate scale parameter
		if isa(obj.params.sigma, 'function_handle'),
			obj.params.sigma = obj.params.sigma(data, obj.params);
		end

		% normalise all columns of v to have unit length
		v = bsxfun(@rdivide, obj.v, sqrt(sum(obj.v.^2)));

		% assign cluster labels
		if isempty(clusters),
			clusters = scppNJW(2,v,data,obj.params.sigma,[],obj.params.beta, obj.params.delta);
		end

		[d,pdim] = size(v);
		% if more than two projection directions then use those with maximum variance
		if pdim ==1,
			v2 = pcacomp(data - (data*v)*v', 1);
			if isempty(v2),
				warning('Insufficient observations to compute principal components: Using first two coordinates');
				X = data(:,1:2);
			else
				X = data*[v,v2];
			end
		else
			X = data*v;
			id = [1,2];
			if size(v,2)>2, 
				[~,id] = sort(std(X), 'descend');
			end
			X = X(:, id(1:2));
		end

		% Estimate projection index
		if isempty(obj.fval),
			obj.fval = f_sc(obj.v, data, obj.params);
		end

		subfig = true;
		% This will always return true in Octave
		if isempty(hFig) | ~isa(hFig,'matlab.graphics.axis.Axes'),
			cFig = figure(1);
			clf;
			set(cFig, 'Position', [0 0 1024 1024]);
			hFig = subplot(1,1,1);
			subfig = false;
		end
		hold on;
		if defLabels,
			% ensure reproducibility of colours
			l = unique(labels);
			for i=1:length(l),
				s = find(labels==l(i) & clusters==1);
				plot(hFig, X(s,1), X(s,2),'bo','MarkerSize',3, ...
					'Color', colours(l(i),:));

				s = find(labels==l(i) & clusters==2);
				plot(hFig, X(s,1), X(s,2),'bo','MarkerSize',3, ...
					'Color', colours(l(i),:), 'MarkerFaceColor', colours(l(i),:));
			end
		else
			left = (clusters==1);
			if sum(left)>0,
				plot(hFig, X( left,1), X( left,2), 'bo','MarkerSize',3, 'Color',colours(1,:));
			end
			if sum(~left)>0,
				plot(hFig, X(~left,1), X(~left,2), 'bo','MarkerSize',3, 'Color',colours(2,:));
			end
		end
		%axis([min(X(:,1)), max(X(:,1)), min(X(:,2)), max(X(:,2))]);
		hold off;

		% Figure title
		if ~subfig,
			str = sprintf('(iter: %i) - beta: %1.2f - Eigenvalue: %1.5f ', ...
				obj.params.iter, obj.params.beta, obj.fval);
			if defLabels,
				str = sprintf('%s - SR: %1.3f', str, success_ratio(clusters, labels));
			end
			title(str);
			drawnow;
		end
	end
end

methods(Static)
	function out = def_opar(pars,X)
	%DEF_OPAR: Default optimisation parameters for SCPP
		if nargin < 1,
			pars = struct();
		elseif isempty(pars) | ~isstruct(pars),
			pars = struct();
		end
		out = pars;

		if ~isfield(out, 'iter'), out.iter = 0; end
		if ~isfield(out, 'maxit'), out.maxit = 50; end
		if ~isfield(out, 'ftol'), out.ftol = 1.e-7; end
		if ~isfield(out, 'verb'), out.verb = 0; end
		if ~isfield(out, 'minsize'), out.minsize = 1; end

		if ~isfield(out,'v0'),
			out.v0 = @(x,p)(pcacomp(x,[1,2]));
		end
		if ~isfield(out,'omega'), out.omega = 1.0; end

		if ~isfield(out,'NumMicroClust'), out.NumMicroClust = 200; end
		if ~isfield(out,'weights'), out.weights = []; end

		% Kernel scale parameter (sigma)
		if ~isfield(out,'sigma'), 
			if nargin < 2, 
				out.sigma = @(x,p)(scpp_def_sigma(x));
			else
				out.sigma = scpp_def_sigma(X);
			end
		end 

		% Similarity transformation 
		if ~isfield(out, 'beta'), out.beta = 3.0; end
		if ~isfield(out, 'betaMin'), out.betaMin = 0.5; end
		if ~isfield(out, 'betaStep'), out.betaStep = 0.2; end
		if ~isfield(out, 'delta'), out.delta = 0.01; end
	end

	function [node,idx,spindex] = split(data, pars, labels, colours)
	%SPLIT Bi-partition data through SCPP
		if nargin < 2,
			error('MATLAB:schp:split','split function requires at least 2 inputs');
		end
		if nargin < 4,
			colours = [];
		end
		if nargin < 3,
			labels = [];
		end
		if nargin < 2,
			pars = schp.def_opar();
		end
		[node, idx, spindex] = scpp(data, pars, labels, colours);
	end
end
end
