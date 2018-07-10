function [idx,sol] = ncuth(X, varargin)
%Minimum Normalised Cut Hyperplane
%[IDX,SOL] = NCUTH(DATA, VARARGIN)
%
%  [IDX,SOL] = NCUTH(X) bi-partitions the points in the N-by-D data matrix X with
%  the hyperplane that minimises the normalised cut criterion (computed 
%  on one-dimensional projections of the data).
% 
%  NCUTH returns a vector IDX containing the binary cluster assignment and a
%  minimum Normalised Cut Hyperplane (nchp) object SOL. (If S initial projection
%  vectors are specified S nchp hyperplanes are returned: see v0 option)
%
%  [IDX,SOL] = NCUTH(X, 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional parameters
%  in the form of name/value pairs.
%
%  OPTIONAL PARAMETERS:
%  'v0' - D-by-S matrix of S initial projection vectors.
%   	(default: v0 = pca(X,'NumComponents',1) : First principal component of X)
%
%  'sigma' - positive numeric scaling parameter (sigma)
%	(default: sigma = 100*\sqrt(\lambda_1)*N^(-0.2), where \lambda_1 is the 1st eigenvalue of cov(X))
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
%
%Reference:
%D.P. Hofmeyr. Clustering by Minimum Cut Hyperplanes.
%IEEE Transactions on Pattern Analysis and Machine Intelligence, 39(8):1547-1560, 2017.

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if nargin<1,
	error('Data matrix needs to be specified\n');
end;

% Set default parameters
pars = struct();
v0 = pcacomp(X,1);
pars.v0 = v0;
pars.sigma = ncut_sigma(X,100);
pars.minsize = 1;
pars.maxit = 50;
pars.ftol = 1.e-7;
pars.verb = 0;
pars.labels = [];
pars.colours = [];

pars = myparser(X,2,varargin,pars);

assert(~isa(pars.v0,'function_handle'), 'Initial projection vector/matrix cannot be defined as function handle');
assert(~isa(pars.sigma,'function_handle'), 'Scaling parameter (sigma) must be a scalar, not a function handle');

if max(pars.labels)==1,
	pars.labels=[];
end

if pars.verb >0,
	fprintf('\nPARAMETER VALUES:\n');
	fprintf('v0: %s\n', ifelse( max(1 - abs(v0'*pars.v0)) <eps, '1st PC', 'User-defined') );
	fprintf('sigma: %f\n', pars.sigma);
	fprintf('minsize: %i\n', pars.minsize);
	fprintf('maxit: %i\n', pars.maxit);
	fprintf('ftol: %e\n', pars.ftol);
end


% Estimate MDHs for every initial projection vector
v0 = pars.v0;
idx = zeros(size(X,1), size(v0,2));
for i=1:size(v0,2),
	pars.v0 = v0(:,i);
	[sol(i), idx(:,i)] = ncutpp(X,pars,pars.labels,pars.colours);
end
% Cluster labels are {1,2}
idx = 0.5*idx + 1.5;

if size(v0,2) == 1,
	sol = sol(1);
end
