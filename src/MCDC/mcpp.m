function [optHP, idx, spindex] = mcpp(X, pars, labels, colours)
% Maximum Clusterability Projection Pursuit (MDPP) algorithm
%[OPTHP, IDX, SPINDEX] = MCPP(X, PARS, LABELS, COLOURS)
%
% Returns:
%	(OPTHP) Maximum clusterability hyperplane (if more than one initial projection vectors
%		are used then the one that maximises the splitting criterion pars.split_index())
%	(IDX) Binary cluster assignment {-1,1}
%	(SPINDEX) Value of splitting index criterion 
%
% Inputs:
%	(X) N-by-D Data matrix
%	(PARS) Structure containing all parameters of mddc() algorithm
%	(LABELS) True clusters; used only for visualisation (optional)
%	(COLOURS) Colormap matrix used only for visualisation (optional)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

[N,dim] = size(X);
invalid = false;
if dim == 1,
	fprintf('Effective dimensionality of 1: No point in projection pursuit\n');
	invalid = true;
end
if 2*pars.minsize > N,
	fprintf('Too few observations to split\n');
	invalid = true;
end

% Cluster membership
idx = -ones(N,1);
% Split index
spindex = -inf;
% If split constraints are not met
optHP = mchp(ones(dim,1)/sqrt(dim),0,-inf,pars);
if invalid,
	return;
end

if nargin < 4,
	colours = [];
end

if nargin>2 & ~isempty(labels),
	labels = fixLabels(labels);
else
	labels = [];
end


% In the absence of HP selection criterion use Variance Ratio Clusterability index
if ~isfield(pars,'split_index') | isempty(pars.split_index),
	crit = 2;
elseif isa(pars.split_index, 'function_handle'),
	crit = 3;
elseif strcmp(pars.split_index,'size'),
	crit = 1;
elseif strcmp(pars.split_index,'fval'),
	crit = 2;
else
	% If the splitting criterion is incorrectly specified use default
	crit = 2;
end

% Set initial projection vector(s)
if ~isfield(pars,'v0') | isempty(pars.v0),
	pars.v0 = mc_v0(X);
	V = pars.v0(X);
elseif isa(pars.v0,'function_handle'), 
	V = pars.v0(X);
else
	V = pars.v0;
end


% initialise empty MCHP
total_iter = 0;
for i = 1:size(V,2),
	sol = mcpp_internal(X, V(:,i), pars, labels, colours);
	total_iter = total_iter + sol.params.iter;
	if (crit==1 | crit==2) & sol.fval > spindex, 
		% fval is the negative of the Variance Ratio clusterability index
		spindex = -sol.fval;
		optHP = sol;
	elseif crit==3,
		% User defined projection index
		pri = pars.split_index(sol.v, X, pars);
		if pri > spindex;
			spindex = pri;
			optHP = sol;
		end
	end
end

if ~isempty(optHP) > 0,
	% Set cluster assignment
	idx(X* optHP.v > optHP.b) = 1;

	% Set split index when selection criterion is Size
	if crit == 1,
		spindex = N;
	end

	optHP.params.iter = total_iter;
	if pars.verb > 0,
		plot(optHP,X,labels,colours);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [optHP, exitFlag]  = mcpp_internal(Data, v, pars, labels, colours)
%function [optHP, exitFlag] = mcpp_internal(v, Data, labels, pars)
%
% Projection pursuit algorithm for Maximum Clusterability Divisive Clustering
% Returns:
%	(optHP) Optimal Maximum Clusterability 
%	(exitFlag) True if a valid MCHP has been identified
%
% Inputs:
%	(Data) Data matrix
%	(v) Initial projection vector
%	(pars) Parameters of MCDC clustering algorithm:
%	      maxit: Maximum number of iterations
%	      minsize: Minimum cluster size
%	      ftol: Function value tolerance
%	(labels) True cluster assignment (Optional: only used for visualisation)
%	(colours) Colours for each cluster (Optional: only used for visualisation) 

% We use this scaling to avoid terminating in flat regions of f
objfungrad = @(x)scaled_fdf(x,Data,pars);

outF = @(a,b,c)outFcn(a,b,c,Data,pars,labels,colours);
if ~isOctave(),
	options = optimoptions('fminunc','Algorithm','quasi-newton','GradObj','on', ...
		'HessUpdate','bfgs', 'MaxIter',pars.maxit,'OutputFcn', outF, 'Display','off');
else
	options = optimset('GradObj','on','MaxIter',pars.maxit,'TolX',pars.ftol,...
		'TolFun',pars.ftol, 'OutputFcn', outF, 'Display','Off');
end

global myiterations;
pars.iter = 0;

% Optimise objective function
[v_opt, fval, exitFlag] = fminunc(objfungrad, v, options);
v_opt = v_opt./norm(v_opt,2);
pars.iter = myiterations;

[f, bmin] = f_mc(v_opt,Data,pars.minsize);
optHP = mchp(v_opt, bmin, f, pars);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,df] = scaled_fdf(v,X,pars)
[f,df] = f_df_mc(v,X,pars.minsize);
f = 1000*f;
df = 1000*df;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stop = outFcn(x, optimValues, state, Data, pars, labels, colours)
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
if ~strcmp(state,'init'),
	stop = ifelse(abs(optimValues.fval - Fprev) < pars.ftol, true, false);
	Fprev = optimValues.fval;
end

if pars.verb>0,
	x = x./norm(x,2);
	pars.iter = myiterations;
	[f,b] = f_mc(x,Data,pars.minsize);
	switch state,
	case 'init'
		hpi = mchp(x,b,f,pars);
		plot(hpi,Data,labels,colours);
	case 'iter'
		hpi = mchp(x,b,f,pars);
		plot(hpi,Data,labels,colours);
	end
end
