function [optHP, idx, spindex] = ncutpp(Data, pars, labels, colours)
%Minimum normalised cut projection pursuit
%[OPTHP, IDX, SPINDEX] = NCUTPP(X, PARS, LABELS, COLOURS)
%
% Returns:
%	(OPTHP): Minimum normalised cut hyperplane (if more than one initial projection vectors
%		are used then the one that maximises the splitting criterion pars.split_index())
%	(IDX): Binary cluster assignment {-1,1}
%	(SPINDEX): Value of splitting index criterion 
%
% Inputs:
%	(X): Data matrix
%	(PARS): Structure containing all parameters of ncutdc() algorithm
%	(LABELS): True clusters; only used for visualisation (optional)
%	(COLOURS): Colormap matrix: only used for visualisation (optional)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

[N,dim] = size(Data);
invalid = false;
if dim == 1,
	fprintf('Effective dimensionality of 1: No point in projection pursuit\n');
	invalid = true;
end
if 2*pars.minsize > N,
	fprintf('Too few observations to split\n');
	invalid = true;
end

idx = -ones(size(Data,1),1);
optHP = nchp(ones(dim,1)/sqrt(dim),0,inf,pars);
spindex = -inf;
% If split constraints are not met
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

% In the absence of HP selection criterion use Eigenvalue as split index
if ~isfield(pars,'split_index') 
	crit = 2;
elseif isempty(pars.split_index),
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
if ~isfield(pars,'v0'),
	pars.v0 = pcacomp(Data,1);
	V = pars.v0(Data);
elseif isempty(pars.v0),
	pars.v0 = pcacomp(Data,1);
	V = pars.v0(Data);
elseif isa(pars.v0,'function_handle'), 
	V = pars.v0(Data,pars);
else
	V = pars.v0;
end

% Set scaling parameter
if isa(pars.sigma,'function_handle'),
	pars.sigma = pars.sigma(Data);
end

% Projection pursuit over all initialisations of projection vector
total_iter = 0;
for i = 1:size(V,2),
	sol = ncutpp_internal(Data, V(:,i), pars, labels, colours);
	total_iter = total_iter + sol.params.iter;
	if (crit==1 | crit==2) & -sol.fval > spindex,
		spindex = -sol.fval;
		optHP = sol;
	elseif crit==3,
		% User defined projection index
		pri = pars.split_index(sol.v, Data, pars);
		if pri > spindex;
			spindex = pri;
			optHP = sol;
		end
	end
end

% If valid hyperplane separator has been obtained
if ~isempty(optHP) > 0,
	% Set cluster assignment
	idx(Data* optHP.v > optHP.b) = 1;

	% Set split index when selection criterion is Size
	if crit == 1,
		spindex = N;
	end

	optHP.params.iter = total_iter;
	if pars.verb > 0,
		plot(optHP,Data,labels,colours);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [optHP, exitFlag] = ncutpp_internal(Data, v, pars, labels, colours)
%function [optHP, exitFlag] = ncutpp_internal(Data, v, pars, labels, colours)
%
% Projection pursuit algorithm for Normalised Cut Hyperplane Divisive Clustering
% Returns:
% Returns:
%	(optHP) Optimal Maximum Clusterability 
%	(exitFlag) True if a valid NCHP has been identified
%
% Inputs:
%	(Data) Data matrix
%	(v) Initial projection vector
%	(labels) True cluster assignment (only used for visualisation) (if unknown this set to ones)
%	(pars) Structure array containing parameters of NCUTDC:
%		(sigma) Scaling parameter for Laplace kernel
%		(minsize) Minimum cluster size
%		(maxit) Maximum number of iterations
%		(ftol) Function value tolerance
%		(nc) Number of clusters
%		(verb) Verbosity level

objfungrad = @(x)f_df_ncut(x,Data,pars);

outF = @(a,b,c)outFcn(a,b,c, Data, pars, labels, colours);
if ~isOctave(),
	options = optimoptions(@fminunc,'Algorithm','quasi-newton','GradObj','on', ...
		'HessUpdate','bfgs', 'MaxIter',pars.maxit, 'OutputFcn', outF, 'Display','Off');
else
	options = optimset('GradObj','on','MaxIter',pars.maxit,'TolX',pars.ftol,...
		'TolFun',pars.ftol, 'OutputFcn', outF, 'Display','Off');
end

global myiterations;

[v_opt, fval, exitFlag] = fminunc(objfungrad, v, options);
v_opt = v_opt./norm(v_opt,2);

[f, bmin] = f_ncut(v_opt, Data, pars);
pars.iter = myiterations;

optHP = nchp(v_opt,bmin,f,pars);

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
	stop = ifelse(abs(optimValues.fval - Fprev)< pars.ftol, true, false);
	Fprev = optimValues.fval;
end

if pars.verb > 0,
	x = x./norm(x,2);
	pars.iter = myiterations;
	[f,b] = f_ncut(x,Data,pars);
	switch state,
	case 'init'
		hpi = nchp(x,b,f,pars);
		plot(hpi, Data, labels, colours);
	case 'iter'
		hpi = nchp(x,b,f,pars);
		plot(hpi, Data, labels, colours);
	end
end
