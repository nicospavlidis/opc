function [L, v1] = LaplacianN(W,weights)
%Normalised Laplacian and vector containing square root of degree of each vertex
%[L, V1] = LAPLACIANN(W,WEIGHTS)
%
% Inputs:
%	(W): Similarity matrix
%	(WEIGHTS): microcluster weights (if empty no microclustering is assumed)
%
% Output:
%	(L): Normalised Laplacian matrix
%	(V1): unit-length vector containing square root of degrees := 1st eigenvector of (L)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if ~isempty(weights),
	W = (weights*weights') .* W;
end
%u1 = sqrt(sum(W,2));
%%sq_rootD = 1./sqrt(sum(W,2));
%L = eye(size(W,1)) - diag(1./u1)*W*diag(1./u1);

u1 = sqrt(sum(W,2));
v1 = u1./sqrt(sum(u1.^2));
u1 = 1./u1;
L = eye(size(W,1)) - (u1*u1').*W;
end

