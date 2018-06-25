%Generic binary cluster separator class

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

classdef gsep
properties (SetAccess = public)
	% Projection matrix: columns form orthonormal basis
	v;
	% projection index
	fval;

	% row-numbers of observations assigned to specific cluster
	% (sign determines whether observation is allocated to the left (-) or right (+) child)
	idx;

	% function handle to projection pursuit function
	ppfun;

	% Parameters of projection pursuit algorithm for this separator
	params;

	% Parameters employed by hierarchical divisive clustering algorithms
	tree_params;
end

methods
	function obj = gsep(v1,f1,id1,ppf,p1,t1)
	% CONSTRUCTOR

		if nargin < 1,
			obj.v = [];
			obj.fval = [];
			obj.idx = [];
			obj.ppfun = [];
			obj.params = [];
			obj.tree_params = [];
			return;
		end

		if isa(v1,'gsep'),
			% Copy constructor
			obj.v = v1.v;
			obj.fval = v1.fval;
			obj.idx = v1.idx;
			obj.ppfun = v1.ppfun;
			obj.params = v1.params;
			obj.tree_params = v1.tree_params;
			return;
		end

		if nargin < 6, t1 = struct(); end;

		if nargin < 5, p1 = []; end;

		if nargin < 4, ppf = []; end

		if nargin < 3, 
			id1 = [];
		elseif ~isempty(id1),
			assert(size(id1,2)==1 & ~sum(id1~=floor(id1)), ...
				'idx: must be a column vector of integers');
		end;


		if nargin < 2, 
			f1 = [];
		elseif ~isempty(f1),
			assert(numel(f1)==1 & isreal(f1),'f: must be a real number');
		end

		assert(size(v1,1)>1 & ~sum(isnan(v1)) & ~sum(isinf(v1)), 'v: must be a matrix of reals');

		% ensure all column vectors are unit length
		obj.v = bsxfun(@rdivide, v1, sqrt(sum(v1.^2)));
		obj.fval = f1;
		obj.idx = id1;
		obj.ppfun  = ppf;
		obj.params = gsep.def_opar(p1);

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

	function [v,w] = tree_plot_data(obj,X)
	% TREE_PLOT_DATA Returns two dimensional subspace used for visualisation
		if isempty(obj.v), error('Cannot plot uninitialised GSEP object'); end

		v = bsxfun(@rdivide, obj.v, sqrt(sum(obj.v.^2)));
		[d,pdim] = size(v);
		if pdim ==1,
			w = pcacomp(X - (X*v)*v', 1);
			if isempty(w),
				warning('Insufficient observations to compute principal components: Using first two coordinates');
				w = [0; 1; zeros(d-2,1)];
				v = [1; zeros(d-1,1)];
			end
		else
			[~,id] = sort(std(X*v), 'descend');
			w = v(:,id(2));
			v = v(:,id(1));
		end
	end

	function hFig = plot(obj,data,labels,colours,hFig)
	%PLOT Visualisation of binary partition
		if isempty(obj.v), error('Cannot plot uninitialised GSEP object'); end
		if isempty(obj.idx), error('Cannot plot GSEP with no cluster information'); end

		if nargin<2, error('Data matrix undefined'); end

		if size(data,2) ~= size(obj.v,1), error('Data matrix is incompatible with projection vector'); end

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

		index = abs(obj.idx);
		data = data(index,:);
		labels = labels(index);

		% normalise all column vectors to unit-length
		v = bsxfun(@rdivide, obj.v, sqrt(sum(obj.v.^2)));
		[d,pdim] = size(v);
		% if more than two projection directions then use those with maximum variance
		if pdim == 1,
			v2 = pcacomp(data - (data*v)*v', 1);
			if isempty(v2),
				warning('Insufficient observations to compute principal components: Using first two coordinates');
				X = data(:,1:2);
			else
				v1 = v;
				X = data*[v1,v2];
			end
		else
			X = data*v;
			id = [1,2];
			if size(v,2)>2, 
				[~,id] = sort(std(X), 'descend');
			end
			X = X(:, id(1:2));
			v1 = v(:,id(1));
			v2 = v(:,id(2));
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
				s = find(labels==l(i) & obj.idx<0);
				plot(hFig, X(s,1), X(s,2),'bo','MarkerSize',3, 'Color', colours(l(i),:));

				s = find(labels==l(i) & obj.idx>0);
				if pdim > 1,
					plot(hFig, X(s,1), X(s,2),'bo','MarkerSize',3, ...
						'Color', colours(l(i),:), 'MarkerFaceColor', colours(l(i),:));
				else
					plot(hFig, X(s,1), X(s,2),'bo','MarkerSize',3, 'Color', colours(l(i),:));
				end
			end
		else
			left = (obj.idx<0);
			if sum(left)>0,
				plot(hFig, X( left,1), X( left,2), 'bo','MarkerSize',3, 'Color',colours(1,:));
			end
			if sum(~left)>0,
				plot(hFig, X(~left,1), X(~left,2), 'bo','MarkerSize',3, 'Color',colours(2,:));
			end
		end
		%axis([min(X(:,1)), max(X(:,1)), min(X(:,2)), max(X(:,2))]);

		% For one-dimensional projections visualise split point and marginal densities
		if pdim==1 & size(X,1)>1,
			% split point along projection vector
			if length(unique(sign(obj.idx))) > 1,
				if max(X(obj.idx < 0, 1)) < min(X(obj.idx > 0,1)), 
					bmin = 0.5*(max(X(obj.idx<0,1)) + min(X(obj.idx>0,1)));
				else 
					bmin = 0.5*(min(X(obj.idx<0,1)) + max(X(obj.idx>0,1)));
				end
			else
				bmin = [];
			end

			XLim = get(gca,"XLim");
			y = linspace(XLim(1), XLim(2), 200)';

			% Bandwidth selection
			if ~isfield(obj.params,'bandwidth'),
				h = 0.9* std(X(:,1)) * size(X,1)^(-0.2);
			elseif isa(obj.params.bandwidth,'function_handle'),
				h = obj.params.bandwidth(data,obj.params);
			elseif numel(obj.obj.bandwidth)==1 & obj.params.bandwidth>0,
				h = obj.params.bandwidth;
			end

			f = fgt_kde(X(:,1),y,h);

			ax1_pos = get(gca,"position");
			ax2 = axes('Position',ax1_pos,...
				'XAxisLocation','top',...
				'YAxisLocation','right',... 
				'XLim',XLim,...
				'XTick',[],...
				'YLimMode','manual',...
				'YLim', [0, 2*max(f)],...
				'Color','none');

			% Plot KDE for projected density
			line(y, f, 'Parent',ax2,'Color','k','LineWidth',2);
			%plot split point
			if ~isempty(bmin), 
				line([bmin bmin], [0, 2*max(f)], 'Parent', ax2, 'Color','r'); 
			end
			if defLabels,
				binL = binAssign(labels, sign(obj.idx));
				% if not all observations are assigned to the same cluster
				if length(unique(binL)) > 1,
					fi = (size(X(binL==1,1),1)/size(X,1)) * fgt_kde(X(binL==1,1),y,h);
					line(y, fi,  'Parent',ax2, 'Color','r','LineWidth',2, 'LineStyle','--'); 
					line(y, f-fi,'Parent',ax2, 'Color','b','LineWidth',2, 'LineStyle','--');
				end
			end
		end
		hold off;

		% Figure title
		if ~subfig,
			str = '';
			if ~isempty(obj.fval),
				str = sprintf('Projection Index: %1.5f ', obj.fval);
			end
			if defLabels,
				str = sprintf('%s - SR: %1.3f', str, success_ratio(sign(obj.idx), labels));
			end
			if ~isempty(str),
				title(str);
			end
			drawnow;
		end
	end % END PLOT

	function [node,idx,spindex] = split(obj, data, pars, labels, colours)
	%SPLIT Interface for splitting clusters for generic Projection Pursuit methods
		if nargin < 2,
			error('MATLAB:gsep:split','split function requires at least 2 inputs');
		end
		if nargin < 3,
			pars = obj.params;
		elseif isempty(pars),
			pars = obj.params;
		end

		% Apply generic projection pursuit algorithm
		[node, idx, spindex] = gpp(data, obj.ppfun, pars);
	end % END 
end

methods(Static)
	function out = def_opar(pars)
	%DEF_OPAR Default Optimisation parameters
		if nargin < 1 
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
end
end
