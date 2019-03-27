function [coeff,index] = rpcacomp(X,index)
%Returns the robust principal components of (X) specified in vector (index) 
%using the Principal Component Projection Puruit algorithm by Alternating Directions
%
%[COEFF] = RPCACOMP(X,INDEX)
%
% Inputs:
%	(X): Data matrix
%	(INDEX): Vector containing indices of principal component vectors required
%
% Output:
%	(COEFF): Robust Principal Components ordered in the sequence specified by (INDEX)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------


if isempty(index),
	error('rpcacomp: Unspecified number of Robust PCs');
elseif min(index<=0),
	error('rpcacomp: all elements of (index) have to be positive integers');
end

% perform Principal Component Pursuit 
[L, S] = pcp(X);

r = rank(L);
if max(index) < r,
	index = index(index<=r);
end

coeff = pcacomp(L,index);
