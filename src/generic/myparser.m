function pars = myparser(X,K,args,defp)
%Function used to parse optional arguments in form of Name,Value pairs for a number of OPC algorithms
%PARS = MYPARSER(X,K,ARGS,DEFP)
%
% Returns:
%	(PARS): Structured array containing optional user-specified parameters 
%
% Inputs:
%	(X): N-by-D data matrix
%	(K): Number of clusters
%	(ARGS): Optional arguments
%	(DEFP): Structured array containing default parameters

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

% Check mandatory inputs: Data matrix
if nargin < 1, 
	error('Valid data matrix needs to be provided');
elseif  size(X,2)<2 | sum(sum((isnan(X)))) | sum(sum((isinf(X)))),
	error('Valid data matrix needs to be provided');
end

% Number of clusters: Data matrix
if nargin < 2, 
	error('Incorrect cluster number specification');
elseif ~isscalar(K) | K<=0 | K~=floor(K),
	error('Incorrect cluster number specification');
end

if (rem(length(args),2)==1)
	error('Optional parameters should always go by pairs');
end
pars = defp;

for i=1:2:(length(args)-1),
	if ~ischar(args{i}),
		error('Unknown type of optional parameter name (parameter names must be strings)');
	end
	pars = setfield(pars,args{i}, args{i+1});	
end

%%%%%%%%%%%%% Parameter checking
[N, dim] = size(X);

% minimum cluster size
if isfield(pars,'minsize'),
	if 2*pars.minsize > N, 
		error('Dataset size (%i) is incompatible with minimum cluster size (%i)',N,pars.minsize);
	end
end

% Cluster labels
nc = 1;
if isfield(pars,'labels'),
	if isempty(pars.labels),
		% warning('MATLAB:myparser: No need to undefine labels vector');
		pars.labels = ones(N,1);
	elseif (size(pars.labels,1) ~= N | size(pars.labels,2)~=1),
		error('Dimensions of label vector inconsistent with data matrix');
	elseif sum(floor(pars.labels)~= pars.labels),
		error('Non integer values in label vector');
	else 
		pars.labels = fixLabels(pars.labels);
		nc = max(pars.labels);
	end
end

% Colours
if isfield(pars,'colours'),
	pars.colours = palette(nc,pars.colours);
end

% maximum number of iterations
if isfield(pars,'maxit'),
	if pars.maxit<0,
		error('Minimum number of iterations has to be non-negative');
	end
end

% objective function tolerance
if isfield(pars,'ftol'),
	if pars.ftol<0,
		error('Tolerance level has to be non-negative');
	end
end

% If structure 'param' is specified transfer all its contents to stucture 'pars'
if isfield(pars, 'param'),
	if isstruct(pars.param),
		fields = fieldnames(pars.param);
		for i=1:numel(fields),
			pars.(fields{i}) = pars.param.(fields{i});
		end
		pars = rmfield(pars,'param');
	end
end

% PARSING OF ARGUMENTS THAT REQUIRE PARAMETER STRUCTURE AS INPUT

% Initial projection matrix
if isfield(pars,'v0'),
	if isa(pars.v0, 'function_handle'),
		v = pars.v0(X,pars);
		if size(v,1) ~= dim | norm(v'*v - eye(size(v,2)), 1) > sqrt(eps), 
			error('Inappropriate initialisation of projection matrix');
		end
	elseif isempty(pars.v0) | size(pars.v0,1) ~= dim | norm(pars.v0'*pars.v0 - eye(size(pars.v0,2)), 1) > sqrt(eps),
		error('Inappropriate initialisation of projection matrix');
	end
end


%% BANDWIDTH SELECTION -- there is no point in using a user-defined value here!
if isfield(pars,'bandwidth'),
	if isa(pars.bandwidth,'function_handle'),
		h = pars.bandwidth(X,pars);
		if ~isscalar(h) | h<=0,
			error('Invalid function for bandwidth selection');
		end
	elseif isscalar(pars.bandwidth) & pars.bandwidth>0,
		%mult = pars.bandwidth;
		%pars.bandwidth = @(x,p)(mult* size(x,1)^(-0.2) * std(x* pcacomp(x,1)));
	else
		error('Incorrectly specified argument');
	end
end

% Scaling parameter used in Spectral clustering to estimate Gaussian kernels
if isfield(pars,'sigma'),
	if isa(pars.sigma,'function_handle'),
		sigma = pars.sigma(X,pars);
		if numel(sigma)~=1 | sigma<=sqrt(eps),
			error('Invalid function for scale parameter selection (sigma)');
		end
	elseif numel(pars.sigma)==1 & pars.sigma>0,
	else
		error('Scale parameter needs to be defined as function handle');
	end
end

% Cluster splitting criterion
if isfield(pars,'split_index'),
	if isempty(pars.split_index),
		error('Split index function must be defined');

	elseif isa(pars.split_index, 'function_handle'),
		%Testing if split index function works correctly
		s = pars.split_index(ones(dim,1)/sqrt(dim), X, pars);
		if ~isscalar(s) | isinf(s) | isnan(s),
			error('Cluster splitting index function is not working correctly');
		end

	elseif ~ischar(pars.split_index), 
		error('Cluster splitting index is not valid');

	elseif strcmp(pars.split_index,'fval'),
		% Default setting nothing to do

	elseif strcmp(pars.split_index,'size'),
		pars.split_index = @(v,x,p)(size(x,1));

	elseif strcmp(pars.split_index,'scatter'),
		pars.split_index = @(v,x,p)(total_scatter(x));

	elseif strcmp(pars.split_index,'rdepth'),
		pars.split_index = @(v,x,p)(md_reldepth(v,x,p));

	else
		error('Cluster splitting index is not valid');
	end
end
