function [L, S] = pcp(X, lambda, rho, tol, maxit)
%Principal Component Pursuit by Alternating Direction Method of Multipliers algorithm
%
% Inputs:
%	(X): N-by-D data matrix
%	(lambda): Regularisation parameter ( default: lambda = 1/sqrt(max(N,D)) )
%	(rho): Penalty parameter for Augmented Lagrangian ( default: rho = 0.25*N*D/norm(X(:),1) )
%	(tol): Tolerance ( default: tol = 1.0e-7 )
%	(maxit): Maximum number of iterations ( default: maxit=1000 )
%
% Outputs:
%	(L): Low-rank component of X
%	(S): Noise/ Perturbation matrix
%
%Reference:
%E.J. Candes, X. Li, Y. Ma and J. Wright. Robust Principal Component Analysis?
%Journal of the ACM, 58(3): Article 11, 2011.


%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

[n, d] = size(X);

% default arguments
if nargin < 5 || isempty(maxit),
    maxit = 1000;
end
if nargin < 4 || isempty(tol),
    tol = 1e-7;
end
if nargin < 3 || isempty(rho),
    rho = 0.25*n*d/norm(X(:),1);
end
if nargin < 2 || isempty(lambda),
    lambda = 1/sqrt(max(n,d));
end

% constants
ir = 1./rho;
lr = lambda*ir;
tol = tol*norm(X,'fro');

% initialisation
L = zeros(n,d);
S = zeros(n,d);
Y = zeros(n,d);
for i = 1:maxit,
	% Update L
	[U,S,V] = svd(X - S + ir*Y, 'econ');
	L = U * shrinkage(S,ir) * V';

	% Update S 
	S = shrinkage(X - L + ir*Y, lr);

	% Primal Residual
	R = X - L - S;

	% Update Lagrange multipliers
	Y = Y + rho*R;

	% check termination criterion
	if norm(R, 'fro') <= tol,
		break; 
	end
end
end

function z = shrinkage(x,kappa)
z = max(x - kappa, 0) - max(-x - kappa, 0);
end

function r = Do(X,kappa)
    % shrinkage operator for singular values
    [U, S, V] = svd(X, 'econ');
    Sr = shrinkage(S,kappa);
    r = U*shrinkage(S,kappa)*V';
end
