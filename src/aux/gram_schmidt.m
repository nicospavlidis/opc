function v = gram_schmidt(v,W)
%Gram-Schmidt orthonormalisation of vector (v) with respect to columns of (W)
%V = GRAM_SCHMIDT(V,W)
%
% Inputs:
%	V: D-by-1 vector 
%	W: D-by-N Matrix whose columns are orthonormal vectors


if isempty(W),
	v = v./norm(v,2);
else
	% Normalise columns of W to unit-length
	W = bsxfun(@rdivide, W, sqrt(sum(W.^2)));
	for i = 1:size(W,2),
		v = v - (v'*W(:,i))* W(:,i);
	end
	v = v./norm(v,2);
end
