function [L, v1] = LaplacianN(W,weights)
%function [L, v1] = LaplacianN(W,weights)
% Returns:
%	(L): Normalised Laplacian matrix
%	(v1): unit-length vector containing square root of degrees := 1st eigenvector of (L)
% Inputs:
%	(W): Similarity matrix
%	(weights): microcluster weights (if empty no microclustering is assumed)

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

