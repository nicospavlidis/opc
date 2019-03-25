function [optHP, idx, spindex] = mdpp(Data, pars, labels, colours)
%Minimum Density Projection Pursuit (MDPP) algorithm
%[OPTHP, IDX, SPINDEX] = MDPP(X, PARS, LABELS, COLOURS)
%
% Inputs:
%	(X): Data matrix
%	(PARS): Structure containing all parameters of mddc() algorithm
%	(LABELS): True clusters; used only for visualisation (optional)
%	(COLOURS): Colormap matrix used only for visualisation (optional)
%
% Output:
%	(OPTHP): Minimum density hyperplane (if more than one initial projection vectors
%		are used then the one that maximises the splitting criterion pars.split_index())
%	(IDX): Binary cluster assignment {-1,1}
%	(SPINDEX): Value of splitting index criterion 
%

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
if N < 2*pars.minsize + 3,
	fprintf('Too few observations to split\n');
	invalid = true;
end

% Cluster membership
idx = -ones(N,1);
% Split index
spindex = -inf;
% If split constraints are not met
optHP = mdhp(ones(dim,1)/sqrt(dim), 0, inf, 0, pars);
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

% centre Data at zero
mu = mean(Data);
Data = bsxfun(@minus,Data, mu);

% In the absence of HP selection criterion use Relative Depth
crit = 3;
if isfield(pars,'split_index'),
	if ~isempty(pars.split_index),
		if isa(pars.split_index, 'function_handle'),
			crit = 4;
		elseif strcmp(pars.split_index,'size'),
			crit = 1;
		elseif strcmp(pars.split_index,'fval'),
			crit = 2;
		end
	end
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
if size(V,1)~= dim, 
	error('Initial projection matrix has %i rows whereas data matris has %i columns',size(V,1),dim);
end

% Compute bandwidth for this cluster
if isa(pars.bandwidth,'function_handle'),
	pars.bandwidth = pars.bandwidth(Data,pars);
end

if pars.alphamax < pars.alphamin, 
	error('alphamin %1.2f >= alphamax %1.2f',pars.alphamin,pars.alphamax);
end


%% Main loop of projection pursuit algorithm
total_iter = 0;
for i = 1:size(V,2),
	[sol, exitflag] = mdpp_internal(Data, V(:,i), pars, labels, colours);
	total_iter = total_iter + sol.params.iter;

	if exitflag,
		% if cluster splitting criterion is size or relative depth
		if (crit==1 | crit==3) & sol.rel_dep > spindex, 
			spindex = sol.rel_dep;
			optHP = sol;

		% if cluster splitting criterion is density on separating hyperplane
		elseif crit==2 & (-sol.fval > spindex),
			spindex = -sol.fval;
			optHP = sol;
		% if another user defined cluster splitting criterion is used
		elseif crit==4,
			% User defined projection index
			pri = pars.split_index(sol.v, Data, sol.params);
			if pri > spindex;
				spindex = pri;
				optHP = sol;
			end
		end
	end
end

if ~isinf(optHP.fval) > 0,
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
	% account for de-meaning the data
	optHP.b = optHP.b + mu*optHP.v;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [optHP, exitFlag] = mdpp_internal(Data, x0, pars, labels, colours)
%function [optHP, exitFlag] = mdpp_internal(Data, x0, pars, labels, colours)
%
% Minimum Density Projection Pursuit for 1 initial projection vector
% Returns:
%	(optHP) Optimal Minimum Density Hyperplane
%	(exitFlag) True if a valid MDHP has been identified
% Inputs:
%	(X) Data matrix
%	(v) Initial projection vector
%	(pars) Structure array containing all the parameters of MDDC
%	(labels) True cluster assignment (Optional: only used for visualisation)
%	(colours) Colours for each cluster (Optional: only used for visualisation) 

global bestA; bestA = 0;
global bestX; bestX = [];
global bestD; bestD = 0;
global Fprev; Fprev = inf;

% normalise initial projection vector to unit-length
x0 = x0./norm(x0,2);
x = x0;
pars.iter = 0;
for a=[pars.alphamin:0.1:pars.alphamax],
	pars.alpha = a;
	[x,fx,it] = BFGS(x, Data, pars, labels, colours);
	pars.iter = pars.iter + it;
end

if ~isempty(bestX),

	pars.alpha = bestA;
	[fmin, bmin] = f_md(bestX, Data, pars);
	optHP = mdhp(bestX,bmin,fmin,bestD,pars);

	exitFlag = 1;

	if pars.verb>0,
		plot(optHP,Data, labels,colours);
	end
else
	pars.alpha = 0;
	optHP = mdhp([],0,-inf,0,pars);
	exitFlag = 0;
end

clear bestX;
clear bestA;
clear bestD;
clear Fprev;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, fval, iter] = BFGS(x,Data,pars,labels,colours)
global myiterations;

objfungrad = @(v)f_df_md(v,Data,pars);
outF = @(a,b,c)outFcn(a,b,c,Data,pars,labels,colours);
if ~isOctave(),
	options = optimoptions(@fminunc,'Algorithm','quasi-newton','GradObj','on','HessUpdate','bfgs', ...
		'MaxIter', pars.maxit, 'ObjectiveLimit', pars.ftol, 'OptimalityTolerance', pars.ftol, ...
		'Display',ifelse(pars.verb>0,'iter','off'), 'OutputFcn', outF);
else
	options = optimset('GradObj','on','MaxIter',pars.maxit,'TolX',pars.ftol,...
		'TolFun',pars.ftol, 'OutputFcn', outF, 'Display',ifelse(pars.verb>0,'iter','off'));
end

[x, fval, exitflag] = fminunc(objfungrad, x, options);
iter = myiterations;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stop = outFcn(x, optimValues, state, Data, pars, labels, colours)
global bestX;
global bestA;
global bestD;
global Fprev;

% This is necessary because in Octave when the output function causes the optimisation
% to terminate the (output) structure containing convergence statistics is not defined
global myiterations;

if isOctave(),
	myiterations = optimValues.iter;
else 
	myiterations = optimValues.iteration;
end

stop = false;

if ~strcmp(state,'init'),
	stop = ifelse(abs(optimValues.fval - Fprev) < pars.ftol, 1, 0);
	Fprev = optimValues.fval;
end

switch state,
case 'init'
	if pars.verb>0 && pars.iter==0,
		hpi = mdhp(x,[],[],[],pars);
		plot(hpi,Data,labels,colours);

	end
case 'iter'
	x = x./norm(x,2);
	proj = Data*x;
	[fmin, bmin, ~, kdeGrad] = f_md(x,Data,pars);
	bmin = bmin(1);
	dp = reldepth(x, proj, pars);
	if (abs(kdeGrad)<1.0e-6) && (dp>0) && (min(sum(proj<=bmin), sum(proj>bmin)) >= pars.minsize),
		bestX = x;
		bestA = pars.alpha;
		bestD = dp;
	end

	if pars.verb>0, 
		pars.iter = pars.iter + myiterations;
		hpi = mdhp(x,bmin,fmin,dp,pars);
		plot(hpi,Data,labels,colours);
	end
case 'done'
otherwise
end
