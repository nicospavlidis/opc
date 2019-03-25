%Minimum Normalised Cut Hyperplane (inherits from HP class)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

classdef nchp < hp

methods 
	function obj = nchp(v1,b1,f1,p1,id1,t1)
	% CONSTRUCTOR
		if nargin < 1,
			args = {};
		elseif isempty(v1),
			args = {};
		% Copy constructor
		elseif nargin == 1 & isa(v1, 'nchp'),
			args = {v1.v, v1.b, v1.fval, v1.params, v1.idx, v1.tree_params};
		else
			if nargin<6, t1 = []; end;
			if nargin<5, id1 = []; end;
			if nargin<4, p1 = []; end;
			if nargin<3, f1=[]; end;
			if nargin<2, b1 = []; end; 
			args = {v1, b1, f1, p1, id1, t1};
		end

		% Call superclass constructor
		obj = obj@hp(args{:});

		obj.params = nchp.def_opar(obj.params);
	end

	function hFig = plot(obj,data,labels,colours, hFig)
	%PLOT Plots minimum normalised cut hyperplane
		if obj.params.sigma <= 0,
			error('Cannot plot NCHP with undefined scaling parameter (s)');
		end
		if isempty(obj.v),
			error('Cannot plot NCHP with undefined (v)');
		end

		if nargin < 5, hFig = []; end;
		defLabels = true;
		if nargin<4,
			colours = [];
			if nargin <3,
				labels = [];
			end
		end
		if isempty(labels) || (max(labels)==1 && size(colours,1)==2),
			defLabels = false;
		end
		% Data preprocessing
		[data,labels,nc,idx,X,v,w,colours] = preprocess(obj,data,labels,colours);

		% Estimate optimal split point along v
		if isempty(obj.fval) | isempty(obj.b),
			[obj.fval, obj.b] = f_ncut(obj.v, data, obj.params);
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
			M = sparse(1:size(X,1), labels, 1, size(X,1), size(colours,1));
			scatter(X(:,1), X(:,2),14,M*colours);
			%l = unique(labels);
			%for i=1:length(l),
			%	plot(hFig, X(labels==l(i),1), X(labels==l(i),2),'bo','MarkerSize',3, ...
			%		'Color', colours(l(i),:));
			%end
		else
			right = 1 + (X(:,1) > obj.b);
			M = sparse(1:size(X,1), right, 1);
			scatter(X(:,1), X(:,2),14,M*colours(1:2,:));
			%left = (X(:,1) <= obj.b);
			%plot(hFig, X( left,1), X( left,2), 'bo','MarkerSize',3, 'Color',colours(1,:));
			%plot(hFig, X(~left,1), X(~left,2), 'bo','MarkerSize',3, 'Color',colours(2,:));
		end
		y = linspace(min(X(:,1)), max(X(:,1)), 100)';
		axis([ y(1), y(end), min(X(:,2)), max(X(:,2))]);

		% Silverman's bandwidth recommendation for multimodal density
		h = 0.9* std(X(:,1)) * size(X,1)^(-0.2);
		f = fgt_kde(X(:,1), y, h);

		ax1_pos = get(gca,'position');
		ax2 = axes('Position',ax1_pos,...
			'XAxisLocation','top',...
			'YAxisLocation','right',... 
			'XLim',[min(X(:,1)),max(X(:,1))],...
			'YLimMode','manual',...
			'YLim', [0, 2*max(f)],...
			'Color','none');

		% plot marginal density on v_opt
		line(y, f, 'Parent',ax2,'Color','k','LineWidth',2);

		% plot class-conditional densities
		if defLabels,
			idx = ones(size(X,1),1);
			idx(X(:,1)> obj.b) = 2;
			binL = binAssign(labels,idx);

			if length(unique(binL)) > 1,
				fi = (size(X(binL==1,1),1)/size(X,1)) * fgt_kde(X(binL==1,1), y, h);
				line(y, fi,   'Parent',ax2, 'Color','r','LineStyle','--','LineWidth',2); 
				line(y, f-fi, 'Parent',ax2, 'Color','b','LineStyle','--','LineWidth',2); 
			end
		end
		%plot separating hyperplane
		line([obj.b obj.b], [0, 2*max(f)], 'Parent', ax2, 'Color','r'); 
		hold off;

		if ~subfig,
			% Figure title
			str = sprintf('(iter: %i) - Normalised Cut: %1.5f ', obj.params.iter, obj.fval);
			if nc>1,
				str = sprintf('%s - SR: %1.3f', str, success_ratio(idx, labels));
			end
			title(str);
			drawnow;
		end
	end
end

methods(Static)
	function out = def_opar(pars)
	%DEF_OPAR Default optimisation parameters for normalised cut hyperplane projection pursuit
		if nargin < 1,
			pars = struct();
		end
		% Default parameters for all subclasses of HP
		out = hp.def_opar(pars);

		% Default initialisation NCHP is 1st principal component
		if ~isfield(out, 'v0'),
			out.v0 = @(x)(pcacomp(x,1));
		end
		if ~isfield(out, 'sigma'),
			out.sigma = @(x)(ncut_sigma(x,100));
		end
		out.ftol = 1.0e-7;
	end

	function [node,idx,spindex] = split(data, pars, labels, colours)
	%SPLIT Performs binary split through minimum normalised cut hyperPlane (NCHP)
		if nargin < 2,
			error('MATLAB:nchp:split','split function requires at least 2 inputs');
		end
		if nargin < 4,
			colours = [];
		end
		if nargin < 3,
			labels = [];
		end
		if nargin < 2,
			pars = nchp.def_opar();
		end
		[node, idx, spindex] = ncutpp(data, pars, labels, colours);
	end
end
end
