%MDHP implements Minimum Density Hyperplane (inherits from HP class)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

classdef mdhp < hp

properties (SetAccess = public)
	% Relative depth
	rel_dep = 0;
end

methods 
	function obj = mdhp(v1,b1,f1,rd1,p1,id1,t1)
	% CONSTRUCTOR
		if nargin < 1,
			args = {};
			rd1 = [];
		elseif isempty(v1),
			args = {};
			rd1 = [];
		% Copy constructor
		elseif nargin==1 & isa(v1, 'mdhp'),
			args = {v1.v, v1.b, v1.fval, v1.params, v1.idx, v1.tree_params};
			rd1 = v1.rel_dep;
		else
			if nargin<7, t1 = []; end;
			if nargin<6, id1 = []; end;
			if nargin<5, p1 = []; end;
			if nargin<4, rd1=[]; end;
			if nargin<3, f1=[]; end;
			if nargin<2, b1 = []; end; 
			args = {v1, b1, f1, p1, id1, t1};
		end

		% Call superclass constructor
		obj = obj@hp(args{:});
		obj.rel_dep = rd1;

		obj.params = mdhp.def_opar(obj.params);
	end

	function hFig = plot(obj,data,labels,colours,hFig)
	%PLOT Plots minimum density hyperplane separator
		if isempty(obj.v),
			error('mdhp:plot Cannot plot MDHP with undefined (v)');
		end

		if nargin < 5, hFig = []; end
		% Define colour palette and labels
		defLabels = true;
		if nargin<4,
			colours = [];
			if nargin <3,
				labels = [];
			end
		end
		if isempty(labels),
			defLabels = false;
		end
		% Data preprocessing
		[data,labels,nc,idx,X,v,w,colours] = preprocess(obj,data,labels,colours);

		if isa(obj.params.bandwidth,'function_handle'),
			h = obj.params.bandwidth(data);
		elseif numel(obj.params.bandwidth)==1 & obj.params.bandwidth>0,
			h = obj.params.bandwidth;
		else
			error('mdhp:plot: Invalid bandwidth parameter');
		end

		alpha = obj.params.alpha;
		eta = obj.params.eta;
		epsilon = obj.params.epsilon;

		% Estimate split point
		if isempty(obj.fval) | isempty(obj.b),
			[obj.fval, bmin] = f_md(v, data, obj.params);
			obj.b = bmin(1);
		end
		bmin = obj.b;
		fmin = obj.fval;

		% estimate relative depth
		if isempty(obj.rel_dep),
			obj.rel_dep = reldepth(v, X(:,1), obj.params);
		end
		dp = obj.rel_dep;

		if ~isempty(alpha),
			s = std(X(:,1));
			% Determine limits for plot
			ylow = min(-alpha*s,  min(X(:,1)));
			yup  = max( alpha*s,  max(X(:,1)));
			y = linspace(ylow, yup, 200)';

			f = fgt_kde(X(:,1),y,h);
			mu = mean(X(:,1));
			L = exp(-0.5)/(sqrt(2*pi)*(h*h)* eta^epsilon);
			modF = f + L*max([zeros(length(y),1), y-mu-alpha*s, mu-alpha*s-y],[],2).^(1+epsilon);
		else
			y = linspace(min(X(:,1)), max(X(:,1)), 200)';
			f = fgt_kde(X(:,1), y, h);
			modF = [];
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
			M = sparse(1:size(X,1), labels, 1);
			scatter(X(:,1), X(:,2),14,M*colours);

			%l = unique(labels);
			%for i=1:length(l),
			%	plot(hFig, X(labels==l(i),1), X(labels==l(i),2),'bo','MarkerSize',3,...
			%		'Color', colours(l(i),:));
			%end
		else
			% If a valid separating hyperplane has been located
			if ~isinf(fmin),
				right = 1 + (X(:,1) > obj.b);
				M = sparse(1:size(X,1), right, 1);
				scatter(X(:,1), X(:,2),14,M*colours(1:2,:));
				%left = (X(:,1) <= bmin);
				%plot(hFig, X( left,1), X( left,2), 'bo','MarkerSize',3, 'Color',colours(1,:));
				%plot(hFig, X(~left,1), X(~left,2), 'bo','MarkerSize',3, 'Color',colours(2,:));

			else
				% Otherwise all observations are assigned to the same cluster
				plot(hFig, X(:,1), X(:,2), 'bo','MarkerSize',3, 'Color',colours(1,:));
			end
		end
		axis([y(1), y(end), min(X(:,2)), max(X(:,2))]);

		ax1_pos = get(gca,"position");
		ax2 = axes('Position',ax1_pos,...
			'XAxisLocation','top',...
			'YAxisLocation','right',... 
			'XLim',[min(X(:,1)),max(X(:,1))],...
			'YLimMode','manual',...
			'YLim', [0, 2*max(f)],...
			'Color','none');

		% Plot KDE for projected density
		line(y, f, 'Parent',ax2,'Color','k','LineWidth',2);
		if ~isinf(fmin),
			%plot separating hyperplane
			line([bmin bmin], [0, 2*max(f)], 'Parent', ax2, 'Color','r'); 
			if ~isempty(modF),
				% plot penalised density
				line(y, modF, 'Parent',ax2, 'Color','b','LineStyle','--');
			end
		end
		% plot of clusters assigned the left and right half-space
		if defLabels,
			idx = ones(size(X,1),1);
			idx(X(:,1)> bmin)=2;
			binL = binAssign(labels,idx);

			% if not all observations are assigned to the same cluster
			if length(unique(binL)) > 1,
				fi = (size(X(binL==1,1),1)/size(X,1)) * fgt_kde(X(binL==1,1),y,h);
				line(y, fi,  'Parent',ax2, 'Color','r','LineWidth',2, 'LineStyle','--'); 
				line(y, f-fi,'Parent',ax2, 'Color','b','LineWidth',2, 'LineStyle','--');
			end
		end
		hold off;

		if ~subfig,
			% Figure title
			str = sprintf('(iter: %i)', obj.params.iter);
			str = sprintf('%s a: %1.2f - Dens: %1.3d - RelDepth: %1.3f', str, alpha, fmin, dp);
			if defLabels,
				str = sprintf('%s - SR: %1.3f', str, success_ratio(idx,labels));
			end
			title(str);
			drawnow;
		end
	end
end

methods(Static)
	function out = def_opar(pars)
	%Default Optimisation parameters for Minimum Density Projection Pursuit
		if nargin < 1,
			pars = struct();
		elseif isempty(pars) | ~isstruct(pars),
			pars = struct();
		end
		out = hp.def_opar(pars);

		if ~isfield(out, 'v0'),
			out.v0 = @(x)(pcacomp(x,1));
		end
		if ~isfield(out, 'bandwidth'),
			if isa(out.v0,'function_handle')
				out.bandwidth = @(x)(0.9*size(x,1)^(-0.2)*std(x*out.v0(x)));
			else
				out.bandwidth = @(x)(0.9*size(x,1)^(-0.2)*std(x*out.v0));
			end
		end

		if ~isfield(out, 'alpha'),
			out.alpha = 1;
		end
		if ~isfield(out, 'eta'),
			out.eta = 0.01;
		end
		if ~isfield(out, 'epsilon'),
			out.epsilon = 1-1.0e-6;
		end
		if ~isfield(out, 'alphamin'),
			out.alphamin = 0;
		end
		if ~isfield(out, 'alphamax'),
			out.alphamax = 1;
		end
	end

	function [node,idx,spindex] = split(data, pars, labels, colours)
	%SPLIT Performs binary split through minimum density hyperplane
		if nargin < 2,
			error('mdhp:split','split function requires at least 2 inputs');
		end
		if nargin < 4,
			colours = [];
		end
		if nargin < 3,
			labels = [];
		end
		if nargin < 2,
			pars = mdhp.def_opar();
		end
		[node, idx, spindex] = mdpp(data, pars, labels, colours);
	end
end
end
