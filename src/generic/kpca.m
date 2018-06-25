function [pcv, eigv, rotated] = kpca(K, th)
%Kernel Principal Components Analysis
%[PCV, EIGV, ROTATED] = KPCA(K, TH)
% 
% Inputs: 
%	(K): Kernel Matrix
%	(TH): optional argument with default value 1.e-4. 
%		if (th) is a positive integer >= 1 then th PCs are returned,
%		elseif (th) \in (0,1) Principal components with eigenvalue lower than (th) are ignored.
% Outputs:
%	(PCV): Matrix containing the principal component vectors (stored in columns)
%	(EIGV): Vector of eigenvalues corresponding to each eigenvector
%	(ROTATED): Projections (rotations) on principal components
%
%Reference:
%B. Scholkopf, A. Smola, K.-R. Mueller. Nonlinear component analysis 
%as a kernel eigenvalue problem. Neural Computation 10:1299-1319, 1998.

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if nargin<2, 
	th = 1.0e-4;
elseif isempty(th) || th==0,
	th = 1.0e-4;
end

n = size(K,1);

% Centering 
D = sum(K)/n;
E = sum(D)/n;
J = ones(n,1)*D;
Kc = K - J - J' + E*ones(n,n);

[V,D] = eig(Kc/n);

[eigv, order] = sort(diag(D),'descend');

% Select how many PCs to retain
if th < 1,
	k = sum(eigv>th);
	if k==0, 
		error('threshold used not satisfied by any Kernel PC');
	end
else
	k = floor(th);
end


% eigv = eig in kernlab
eigv = eigv(1:k);
V = V(:, order(1:k));

% the principal component vectors (pcv) are divided by sqrt(lambda_i)
pcv = bsxfun(@rdivide, V, sqrt(sum(V.^2)).*sqrt(eigv)');

if nargout == 3,
	rotated = Kc*pcv; 
end

