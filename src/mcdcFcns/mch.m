function [idx,sol] = mch(X, varargin)
%Maximum Clusterability Hyperplane
%[IDX,SOL] = MCH(X, VARARGIN)
%
%  [IDX,SOL] = MCH(X) bi-partitions the points in the N-by-D data matrix (X) with
%  the hyperplane that maximises the Variance Ratio clusterability criterion.
%
%  MCH returns a vector IDX containing the binary cluster assignment and a
%  Maximum Clusterability Hyperplane (mchp) object SOL. (If S initial projection
%  vectors are specified S maximum clusterability hyperplanes are returned: see 'v0')
%
%  SOL = MCH(X, 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional parameters
%  in the form of name/value pairs.
%
%  OPTIONAL PARAMETERS:
%  'v0' - D-by-S matrix of S initial projection vectors
%   	(default: Vector connecting centroids of 2-means clustering)
%
%  'minsize' - Minimum cluster size (integer)
%	(default minsize = 1)
%
%  'maxit' - Number of BFGS iterations to perform for each value of alpha (default: 50)
%
%  'ftol' - Stopping criterion for change in objective function value over consecutive iterations
%	(default: 1.e-7)
%
%  'verb' - Verbosity. Values greater than 0 enable visualisation during execution
%	(default: 0)
%
%  'labels' - true cluster assignment. Enables the computation of performance over 
%	successive iterations and a better visualisation of how clusters are split
%
%  'colours' - Matrix containing colour specification for observations in different clusters
%	Number of rows must be equal to the number of true clusters (if 'labels' has been specified) or equal to 2.

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if nargin<1, 
	error('Data matrix needs to be specified\n'); 
end;

% Set default parameters
pars = struct();
v0 = mc_v0(X);
pars.v0 = v0;
pars.minsize = 1;
pars.maxit = 50;
pars.ftol = 1.e-7;
pars.labels = [];
pars.colours = [];
pars.verb = 0;

pars = myparser(X,2,varargin,pars);

assert(~isa(pars.v0,'function_handle'), 'Initial projection vector/matrix cannot be defined as function handle');

% If labels are not defined the parses renders this a vector of ones
if max(pars.labels)==1,
	pars.labels = [];
end


[N, dim] = size(X);

if pars.verb >0,
	fprintf('\nPARAMETER VALUES:\n');
	fprintf('v0: %s\n', ifelse(max(1-abs(v0'*pars.v0))==eps, 'Default: Vector connecting 2-means centroids', 'User-defined') );
	fprintf('maxit: %i\n', pars.maxit);
	fprintf('ftol: %e\n', pars.ftol);
	fprintf('minsize: %i\n', pars.minsize);
end


% Estimate MCHs for every initial projection vector
v0 = pars.v0;
idx = zeros(N, size(v0,2));
for i = 1:size(v0,2),
	pars.v0 = v0(:,i);
	[sol(i), idx(:,i)] = mcpp(X, pars, pars.labels, pars.colours);
end

% if only one mchp hyperplane is required
if size(v0,2) == 1,
	sol = sol(1);
end
