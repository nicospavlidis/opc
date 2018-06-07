function [idx,sol] = mdh(X, varargin)
%function [IDX,SOL] = mdh(X, varargin)
%
%MDH: Minimum Density Hyperplane
%  [IDX,SOL] = MDH(X) bi-partitions the points in the N-by-D data matrix X with
%  the hyperplane that minimises the density on a hyperplane criterion
%  (computed from one-dimensional projections of the data). 
%
%  MDH returns a vector IDX containing the binary cluster assignment and a
%  Minimum Density Hyperplane (mdhp) object SOL. 
%
%  [IDX,SOL] = MDH(X, 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional parameters
%  in the form of name/value pairs. 
%
%  OPTIONAL PARAMETERS:
%  'v0' - D-by-S matrix of S initial projection vectors
%   	(default: v0 = pca(X,'NumComponents',1) : First principal component of X)
%
%  'bandwidth' - Bandwidth parameter (default: H = MULT * N^(-0.2) STD(X*PC1), where PC1 is
%		the 1st principal component)
%
%  'alphamin' - The minimum ALPHA over which MDHs are sought: 
%	[mean(X*V) - ALPHA*std(X*V), mean(X*V) + ALPHA*std(X*V)]. 
%	ALPHA starts from (alphamin) and increases by 0.1 every (maxit) iterations until (alphamax) is reached.
%	(default: 0)
%
%  'alphamax' - The maximum ALPHA over which MDHs are sought is 
%	[mean(X*V) - ALPHA*std(X*V), mean(X*V) + ALPHA*std(X*V)]. 
%	ALPHA starts from (alphamin) and increases by 0.1 every (maxit) iterations until (alphamax) is reached.
%	(default: 1)
%
%  'minsize' - Minimum cluster size (integer)
%	(default minsize = 1)
%
%  'maxit' - Number of BFGS iterations to perform for each value of alpha (default: 50)
%
%  'ftol' - Stopping criterion for change in objective function value over consecutive iterations
%	(default: 1.e-5)
%
%  'verb' - Verbosity. Values greater than 0 enable visualisation during execution
%   	 (default: 0)
%
%  'labels' - true cluster assignment. Enables the computation of performance over 
%	successive iterations and a better visualisation of how clusters are split
%
%  'colours' - Matrix containing colour specification for observations in different clusters
%	Number of rows must be equal to the number of true clusters (if 'labels' has been specified) or equal to 2.
%
%Reference:
%N.G. Pavlidis, D.P. Hofmeyr and S.K. Tasoulis. Minimum density hyperplanes.
%Journal of Machine Learning Research, 17(156):1–33, 2016.
%http://jmlr.org/papers/v17/15-307.html.


if nargin<1, 
	error('Data matrix needs to be specified\n'); 
end;

% Set default parameters
pars = struct();
v0 = pcacomp(X,1);
pars.v0 = v0;
h0 = 0.9*size(X,1)^(-0.2)*std(X*pars.v0);
pars.bandwidth = h0;
pars.split_index = 'size';
pars.alphamax = 1.0;
pars.alphamin = 0.0;
pars.minsize = 1;
pars.maxit = 50;
pars.ftol = 1.e-5;
pars.labels = [];
pars.colours = [];
pars.verb = 0;

% We strongly recommend not to modify these
pars.eta = 0.01;
pars.epsilon = 1-1.0e-6;

pars = myparser(X,2,varargin,pars);

assert(pars.alphamin>=0 & pars.alphamax>=pars.alphamin);
assert(~isa(pars.v0,'function_handle'), 'Initial projection vector/matrix cannot be defined as function handle');
assert(~isa(pars.bandwidth,'function_handle'), 'Bandwidth parameter must be a scalar, not a function handle')

[N, dim] = size(X);

if pars.verb >0,
	fprintf('\nPARAMETER VALUES:\n');
	fprintf('v0: %s\n', ifelse( max(1-abs(v0'*pars.v0))<eps, '1st PC', 'User-defined') );
	fprintf('bandwidth: %s - %1.3f\n', ...
		ifelse(abs(h0-pars.bandwidth)<eps, '(default) Silverman rule', 'User-defined'), ...
		pars.bandwidth);
	fprintf('maxit: %i\n', pars.maxit);
	fprintf('ftol: %e\n', pars.ftol);
	fprintf('minsize: %i\n', pars.minsize);
	fprintf('alpha range: (%1.3f, %1.3f)\n', pars.alphamin, pars.alphamax);
end


% Estimate MDHs for every initial projection vector
v0 = pars.v0;
idx = zeros(N, size(v0,2));
for i=1:size(v0,2),
	pars.v0 = v0(:,i);
	[sol(i), idx(:,i)] = mdpp(X,pars,pars.labels,pars.colours);
end
% Cluster labels are {1,2}
idx = 0.5*idx + 1.5;
