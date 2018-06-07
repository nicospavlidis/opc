function out = mypdist2(x,y,distance)
%Computes Euclidean pairwise distances between observations (row vectors) in X and Y
%OUT = MYPDIST2(X,Y,DISTANCE)
%
% Included for compatibility with Octave
% In MATLAB: mypdist2(x,y) = pdist2(x,y,'euclidean')
%
% Vectors (X) and (Y) must be of the same length
if nargin<3 || isempty(distance),
	distance = 'euclidean';
end

if ~isOctave(),
	out = pdist2(x,y);
	return;
end

% compute all pairwise distances
out = squareform(pdist([x;y], distance));
n = size(out,2);
nx = size(x,1);
out = out(1:nx, (nx+1):n);
