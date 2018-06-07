function [optS, idx, spindex] = scpp(Data, pars, labels, colours)
%function [optS, idx, spindex] = scpp(X, pars, labels, colours)
%
% Implements Spectral Clustering Projection Pursuit (SCPP) 
% (THE DIVISIVE CLUSTERING ALGORITHM IS IN scppdc.m)
%
% Returns:
%	(optS) Linear subspace that minimises second eigenvalue of normalised Laplacian constructed
%		from projected data.(
%	(idx) Binary cluster assignment {-1,1}
%	(spindex) Value of splitting index criterion 
%
% Inputs:
%	(X) Data matrix
%	(pars) Structure containing all parameters of scpp() algorithm
%	(labels) True clusters; only used for visualisation (optional)
%	(colours) Colormap matrix: only used for visualisation (optional)

[N,dim] = size(Data);
invalid = false;
if dim == 1,
	fprintf('Effective dimensionality of 1: No point in projection pursuit\n');
	invalid = true;
end

if nargin < 2,
	pars = schp.def_opar(Data);
elseif isempty(pars),
	pars = schp.def_opar(Data);
elseif ~isfield(pars,'minsize'),
	pars.minsize = 1;
end

if 2*pars.minsize > N,
	%fprintf('Too few observations to split\n');
	invalid = true;
end

% Process parameters
if nargin < 4,
	colours = [];
end

if nargin>2 & ~isempty(labels),
	labels = fixLabels(labels);
else
	labels = [];
end

% Set initial projection vector(s)
if ~isfield(pars,'v0') 
	pars.v0 = pcacomp(Data,[1,2]);
	v0 = pars.v0(Data);
elseif isempty(pars.v0),
	pars.v0 = pcacomp(Data,[1,2]);
	v0 = pars.v0(Data);
elseif isa(pars.v0,'function_handle'), 
	v0 = pars.v0(Data,pars);
	if size(v0,1) ~= dim,
		error('scpp','Incorrect function for v0');
	end
else
	v0 = pars.v0;
end
v0 = bsxfun(@rdivide, v0, sqrt(sum(v0.^2,1)));


% Kernel scaling parameter
if ~isfield(pars,'sigma'), 
	% Default scaling rule
	pars.sigma = scpp_def_sigma(Data);

elseif isempty(pars.sigma),
	% Default scaling rule
	pars.sigma = scpp_def_sigma(Data);

elseif isa(pars.sigma, 'function_handle'),

	pars.sigma = pars.sigma(Data,pars);
end
if numel(pars.sigma)~=1 | pars.sigma<=0, 
	error('Inappropriate choice of kernel scale parameter (sigma)');
end


% Micro-clustering 
if ~isfield(pars,'NumMicroClust'),
	pars.NumMicroClust = 200;
elseif isempty(pars.NumMicroClust),
	pars.NumMicroClust = 200;
end

if pars.NumMicroClust < N,
	[X, pars.weights, data2cores] = microcluster(Data, pars.NumMicroClust);
else
	X = Data;
	pars.weights = ones(N,1);
	data2cores = [1:N]';
end

% Penalty function
if ~isfield(pars,'omega') 
	pars.omega = 1;
elseif isempty(pars.omega),
	pars.omega = 1;
else
	if numel(pars.omega)~=1 | pars.omega<=0,
		error('Inappropriate choice of penalty term (omega)');
	end
end

% similarity function parameters
if ~isfield(pars,'beta') 
	pars.beta = 3;
elseif isempty(pars.beta),
	pars.beta = 3;
end

if ~isfield(pars,'betaStep'),
	pars.betaStep = 0.1;
elseif isempty(pars.betaStep),
	pars.betaStep = 0.1;
end

if ~isfield(pars,'betaMin'),
	pars.betaMin = 0.5;
elseif isempty(pars.betaMin),
	pars.betaMin = 0.5;
end

if ~isfield(pars,'delta'),
	pars.delta = 0.1;
elseif isempty(pars.delta),
	pars.delta = 0.1;
end

pars.beta = pars.beta + pars.betaStep;
exitflag = 0;
pars.iter = 0;


% If conditions to split cluster are violated
if invalid,
	% Set default output parameters
	idx = -ones(N,1);
	optS = schp(v0,inf,pars);
	spindex = -inf;
	return;
end

while ~exitflag & (pars.beta > pars.betaMin),

	pars.beta = max(pars.beta - pars.betaStep, pars.betaMin);

	% Evaluate fumction value at initial projection direction
	[f, eigenGap] = f_sc(v0, X, pars);
	% If there are outliers and beta is too large it is possible for
	% the second eigenvalue to be repeated. These are pathological cases
	% that would be immediately discarded by the minimum cluster size constraint
	if f > sqrt(eps) & eigenGap > sqrt(eps),
		% identify optimal linear subspace
		[v_opt, f_opt, it, flag] = BFGS(v0, X, pars, data2cores, Data, labels, colours);

		pars.iter = pars.iter + it;
		% bi-partition micro-clusters centres using spectral clustering
		idx = scppNJW(2, v_opt, X, pars.sigma, pars.weights, pars.beta, pars.delta);
		% if micro-clustering has been performed:
		if size(X,1) ~= N,
			% assign observations to clusters
			idx = reverse_assign(idx, data2cores);
		end

		% Get cluster sizes
		cluster_sizes = hist(idx,[1,2]);
		min_clust = min(cluster_sizes);
		
		% if minimum cluster constraint is satisfied exit
		if min_clust > pars.minsize,
			exitflag = 1;
		end
	end
end

% Projection subspace found that respects minimum cluster size constraint
if exitflag,
	% normalise columns to have unit length
	v_opt = bsxfun(@rdivide, v_opt, sqrt(sum(v_opt.^2)));
	% construct linear subspace object
	optS = schp(v_opt, f_opt, pars);

	% set idx in {-1,1}
	idx = 2*idx - 3;

	% compute cluster splitting index
	if isfield(pars,'split_index') & isa(pars.split_index, 'function_handle'),
		spindex = pars.split_index(v_opt, Data, pars);
		if numel(spindex)~=1 | isnan(spindex) | isinf(spindex),
			error('MATLAB:scpp:', 'split_index function returns invalid output');
		end

	elseif isfield(pars,'split_index') & strcmp(pars.split_index,'fval'),
		spindex = -f_opt;
	else 
		spindex = N;
	end
else
	% Failure to find linear subspace that respects minsize constraint
	idx = -ones(N,1);
	optS = schp(v0,inf,pars);
	spindex = -inf;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spectral Clustering Projection Pursuit
function [v_opt, fval, iter, exitflag] = BFGS(v, X, params, data2cores, Data, labels, colours)
global myiterations;

fcn = @(y)f_df_sc(y,X,params);

outF = @(a,b,c)(outFcn(a,b,c, Data, X, data2cores, params, labels, colours));
if ~isOctave(),
	options = optimoptions(@fminunc,'Algorithm','quasi-newton','GradObj','on','HessUpdate','bfgs', ...
		'MaxIter',params.maxit, 'OutputFcn', outF, 'Display','off', ...
		'ObjectiveLimit',params.ftol,'TolX',params.ftol);
else
	options = optimset('GradObj','on','MaxIter',params.maxit,'TolX',params.ftol,...
		'TolFun',params.ftol, 'OutputFcn', outF, 'Display','Off');
end

[dim,pdim] = size(v);
[v_opt, fval, exitflag] = fminunc(fcn, reshape(v,[dim*pdim,1]), options);
iter = myiterations;

v_opt = reshape(v_opt, [dim, pdim]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stop = outFcn(x, optimValues, state, Data, X, data2cores, params, labels, colours)
global Fprev;

% This is necessary because in Octave when the output function causes the optimisation
% to terminate the (output) structure containing convergence statistics is not defined...
global myiterations;

if isOctave(),
	myiterations = optimValues.iter;
else 
	myiterations = optimValues.iteration;
end

stop = false;

dim = size(Data,2);
pdim = length(x)/dim; % dimensionality of the projection space

% Projection matrix
v = reshape(x, [dim,pdim]);
v = bsxfun(@rdivide, v, sqrt(sum(v.^2)));

% Cluster assignment
idx = scppNJW(2, v, X, params.sigma, params.weights, params.beta, params.delta);
% ensure colours in plots don't change due to ordering of cluster labels
if idx(1) == 2,
	idx = ~(idx-1)+1;
end

% Check minimum cluster size constraint
%	c0: Number of observations assigned to cluster 1
c0 = sum(params.weights(idx==1));
if min(c0, size(Data,1)-c0) < params.minsize,
	if params.verb>1,
		fprintf('stopFcn: beta: %1.2f Minimum cluster %i\n', a, min(c0, N-c0));
	end
	stop = true;
	return;
end

% Function tolerance
if ~strcmp(state,'init'),
	stop = ifelse(abs(optimValues.fval - Fprev)< params.ftol, true, false);
	Fprev = optimValues.fval;
end

if params.verb==0, 
	return; 
end

if size(X,1) ~= size(Data,1),
	% assign observations to clusters 
	idx = reverse_assign(idx, data2cores);
end

params.iter = params.iter + myiterations;
hpi = schp(v, optimValues.fval, params);
switch state,
	case 'init'
		plot(hpi,Data,labels,colours,idx);
	case 'iter'
		plot(hpi,Data,labels,colours,idx);
	case 'done'
end
end
