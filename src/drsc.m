function [idx,W,fval,sumD,iter] = drsc(X, K, sigma, varargin)
%Dimensionality Reduction for Spectral Clustering
%[IDX,W,FVAL,SUMD,ITER] = DRSC(X, K, SIGMA, VARARGIN)
%
% [IDX, W, FVAL, SUMD, ITER] = DRSC(X, K, SIGMA) produces a clustering of the
% N-by-D data matrix (X) into (K) clusters, by identifying the optimal
% (K-1)-dimensional linear subspace to project the data. (SIGMA) is the
% bandwidth parameter for the Gaussian kernel used to estimate the kernel
% matrix.
%
% [IDX, W, FVAL, SUMD, ITER] = DRSC(X, K, SIGMA) returns the cluster assignment,
% (IDX); the projection matrix (W); a vector of values of the projection index
% at each iteration (FVAL); the sum of squared distances to the cluster centres
% in the optimal linear subspace, (SUMD); and finally the iteration at which
% the algorithm terminated, (ITER). If DRSC fails to converge ITER=0
%
% [IDX, W, FVAL ,SUMD, ITER] = DRSC(X, K, S, 'PARAM1',val1, 'PARAM2',val2, ...)
% specifies optional parameters in the form of Name,Value pairs. 
%
% 'v0' - D-by-Q matrix of initial projection vectors. Q determines dimensionality of
%	projection subspace
%	(default: (K-1) first principal components)
%
% 'maxit' - Number of DRSC iterations (default: 50)
%
% 'maxitdim' - Number of gradient descent iterations for each dimension (default: 50)
%
% 'ftol' - Stopping criterion for change in objective function value over consecutive iterations
%	(default: 1.e-4, suggested in Niu et al. AISTATS 2011)
%
% 'verb' - Verbosity. Values greater than 0 enables progress monitoring during execution
%	(default: 0)
%
%Reference:
%D. Niu, J.G. Dy and M.I. Jordan. Dimensionality Reduction for Spectral Clustering.
%Proceedings of the 14th International Conference on Artificial Intelligence and Statistics,
%volume 15 of JMLR W&CP, pages 552-560, 2011. 

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if nargin < 3, 
	error('drsc:','Scaling parameter is either, or incorrectly specified');
elseif isempty(sigma) | ~isscalar(sigma) | sigma < sqrt(eps),
	error('drsc:','Scaling parameter is either, or incorrectly specified');
end

% Set default parameters
pars = struct();
pars.v0 = pcacomp(X,[1:(K-1)]);
pars.maxit = 50;
pars.maxitdim = 50;
pars.ftol = 1.e-4;
pars.verb = 0;
pars.labels = [];

pars = myparser(X,K,varargin,pars);
[N, dim] = size(X);

% maxitdim
if pars.maxitdim < 0,
	error('MATLAB:drsc:','Incorrect specification of maxitdim');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = pars.v0;
if isa(W,'function_handle'),
	W = W(X);
end

% dimensionality of linear subspace
rdim = size(W,2);

fval = inf*ones(pars.maxit,1);
counter = 1;
exitflag = 0;
for iter = 1:pars.maxit,
	% Step 1: Spectral clustering of projected data to obtain eigenvectors and degrees
	if iter==1,
		% Initial W is equal to eye(dim), recommended in Niu et al. AISTATS (2011)
		[~,~,U,degs] = spclust(X,K,'s',sigma);
	else
		[~,~,U,degs] = spclust(X*W,K,'s',sigma);
	end

	% Step 2: Optimise W through Dimension Growth
	W0 = W;
	W = [];
	% Optimise over each dimension of the projection space incrementally
	for i = 1:rdim,
		% initial projection vector for i-th dimension is initialised with
		% outcome of previous iteration. Causes W at iter=1 to equal PCs
		w = gram_schmidt(W0(:,i), W);
		W = [W, w];

		f0 = inf;
		% Eq. (11) in Niu et al. AISTATS (2011)
		[f,grad] = f_df_drsc(X,W,sigma,U,degs);

		if pars.verb>0,
			fprintf('Iteration %i, Dimension %i, Fval %e:\n',iter,i,f); 
		end

		initer = 0;
		while (abs(f0-f)>pars.ftol) & (norm(grad,2)>pars.ftol) & (initer < pars.maxitdim),

			% line-search to determine stepsize
			gamma = drsc_linesearch(X, W, 0.1, f, grad, sigma, U, degs);
			% No ascend stepsize
			if gamma == 0,
				break;
			end

			% Update projection vector: Eq. (7) in Niu et al. AISTATS (2011)
			W(:,i) = sqrt(1-gamma^2)*W(:,i) + gamma*gram_schmidt(grad,W);

			f0 = f;
			% Avoid gradient computation at last iteration
			if initer < pars.maxitdim-1,
				f = f_df_drsc(X,W,sigma,U,degs);
			else
				[f,grad] = f_df_drsc(X,W,sigma,U,degs);
			end
			initer = initer + 1;

			if pars.verb>0, 
				fprintf('\tInner iter %i Fval %e, Improvement %e\n', initer, f, f-f0); 
			end
		end
	end

	% Check for convergence of projection matrix
	if norm(W-W0, inf) < pars.ftol,
		fprintf('Convergenced in iter %i\n', iter);
		exitflag = iter;
		break;
	end

	fval(counter) = f;
	counter = counter + 1;
end
fval = fval(~isinf(fval));
iter = exitflag;

% Perform spectral clustering on data projected onto W
[idx,sumD] = spclust(X*W,K,'s',sigma);
