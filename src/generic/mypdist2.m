function out = mypdist2(x,y)
%Computes Euclidean pairwise distances between observations (row vectors) in X and Y
%OUT = MYPDIST2(X,Y)
%
% Included for compatibility with GNU Octave
% In MATLAB: mypdist2(x,y) = pdist2(x,y,'euclidean')
%
% Vectors (X) and (Y) must be of the same length

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if ~isOctave(),
	out = pdist2(x,y);
	return;
end

% compute pairwise distances between row vectors in 2 matrices
n1 = size(x,1);
n2 = size(y,1);
%out = sqrt(repmat(sum(x.^2,2),1,n2) + repmat(sum(y.^2,2)',n1,1) -2*x*y');
out = sqrt( bsxfun(@plus, bsxfun(@plus, -2*x*y', sum(y.^2,2)'), sum(x.^2,2)) );
