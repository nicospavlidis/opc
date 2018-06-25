%Class implementing generic hyperplane interface

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

classdef hp
properties (SetAccess = public)
	% unit-length vector orthogonal to hyperplane
	v;
	% displacement from the origin
	b;
	% projection index
	fval;
	% Parameters for separating hyperplane
	params;

	% row-numbers of observations assigned to specific cluster 
	% (sign determines whether observation is allocated to the left (-) or right (+) child)
	idx;

	% Parameters employed by hierarchical divisive clustering algorithms
	tree_params;
end

methods
	function obj = hp(v1,b1,f1,p1,id1,t1)
	% CONSTRUCTOR

		if nargin < 1,
			obj.v = [];
			obj.b = [];
			obj.fval = [];
			obj.params = [];
			obj.idx = [];
			obj.tree_params = [];
			return;
		end

		if isa(v1, 'hp'),
			% Copy constructor
			obj.v = v1.v;
			obj.b = v1.b;
			obj.fval = v1.fval;
			obj.params = v1.params;
			obj.idx = v1.idx;
			obj.tree_params = v1.tree_params;
			return;
		end

		if nargin < 6,
			t1 = []; 
		end;

		if nargin < 5, 
			id1 = []; 
		elseif ~isempty(id1),
			assert(numel(id1)>1 & size(id1,1)>1 & size(id1,2)==1 & ~sum(id1~=floor(id1)),'idx: must be a column vector of integers');
		end;

		if nargin < 4,
			p1 = [];
		end;

		if nargin < 3, 
			f1 = [];
		elseif ~isempty(f1),
			assert(numel(f1)==1 & isreal(f1),'f: must be a real number');
		end

		if nargin < 2,
			b1 = [];
		elseif ~isempty(b1),
			assert(numel(b1)==1 & isreal(b1),'b: must be a real number');
		end;

		assert(size(v1,1)>1 & size(v1,2)==1 & isreal(v1) & ~sum(isnan(v1)) & ~sum(isinf(v1)), 'v: must be a column vector');
		obj.v = v1;

		obj.b = b1;
		obj.fval = f1;
		obj.params = p1;

		obj.idx = id1;

		% Common parameters for all univariate projection pursuit algorithms
		obj.params = hp.def_opar(p1);
		
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
	%CLUSTER assign observations to half spaces defined by HP
		if isempty(obj.v) | isempty(obj.b),
			error('MATLAB:hp:cluster','Cannot cluster with non-initialised HP');
		end
		if size(data,2) ~= size(obj.v,1),	
			error('MATLAB:hp:cluster','Data matrix - HP Incompatible dimensions');
		end
		idx = -ones(size(data,1),1);
		idx(data*obj.v > obj.b) = 1;
	end

	function hFig = plot(obj,data,labels,colours,Fig)
	%PLOT Interface for plotting function for different types of hyperplanes
	end

	function [v,w] = tree_plot_data(obj,X)
	% TREE_PLOT_DATA Returns two dimensional subspace used for visualisation
		if isempty(obj.v), error('Cannot plot uninitialised MCHP object'); end
		v = obj.v;
		v = v./norm(v,2);
		w = pcacomp(X - (X*v)*v', 1);
		if isempty(w),
			warning('Insufficient observations to compute principal components: Using first two coordinates');
			w = [0; 1; zeros(size(X,2)-2,1)];
			v = [1; zeros(size(X,2)-1,1)];
		end

	end

	function [D,L,nc,idx,X,v,w,colours] = preprocess(obj,data,labels,colours)
	%PREPROCESS Preprocess data-matrix, labels, and colours (used in visualisation)
		if isempty(obj.v),
			error('Cannot plot non-initialised HP');
		end

		if nargin<2,
			error('Data matrix undefined');
		end

		if size(data,2) ~= length(obj.v),
			error('Data matrix is incompatible with projection vector');

		end
		% if labels are undefined define all observations to belong to the same cluster
		if nargin<3,
			labels = ones(size(data,1),1);
			nc = 1;
		elseif isempty(labels),
			labels = ones(size(data,1),1);
			nc = 1;
		else
			assert(length(labels)==size(data,1), ...
				'Labels vector is incompatible with data matrix');

			% Ensure labels are consecutive integers {1,...,K}
			labels = fixLabels(labels);
			% number of clusters
			nc = max(labels);
		end

		% Set colours
		if nargin < 4,
			colours = [];
		end;
		colours = palette(nc,colours);

		if isempty(obj.idx),
			idx = [1:size(data,1)]';
		else
			idx = abs(obj.idx);
		end
		D = data(idx,:);
		L = labels(idx);

		% Direction of maximum variance orthogonal to optimal projection vector, v
		[v,w] = obj.tree_plot_data(D);
		% Data projected onto 2-dimensional space (v,w)
		X = D*[v,w];
	end
end

methods(Static)
	function out = def_opar(pars)
	%DEF_OPAR Default Optimisation parameters
		if nargin < 1,
			pars = struct();
		end
		if isempty(pars) | ~isstruct(pars),
			pars = struct();
		end

		out = pars;
		if ~isfield(out, 'iter'),
			out.iter = 0;
		end
		if ~isfield(out, 'minsize'),
			out.minsize = 1;
		end
		if ~isfield(out, 'maxit'),
			out.maxit = 50;
		end
		if ~isfield(out, 'ftol'),
			out.ftol = 1.e-5;
		end
		if ~isfield(out, 'verb'),
			out.verb = 0;
		end
	end

	function [node,idx,spindex] = split(data, pars, labels, colours)
	%SPLIT Performs binary split through minimum density hyperplane
	end
end
end
