function out = palette(nc,colours)
%Determines colours used for visualisation
%OUT = PALETTE(NC,COLOURS)
%
% Inputs:
%	(NC): Number of clusters (if NC=1 it assumed that true clusters are unknown)
%	(COLOURS): Matrix of colours (optional argument)
%
% Output:
%	(OUT): Matrix of RGB colours

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if isempty(colours),
	if nc > 1,
		colours = hsv(nc);
	else
		% Default is to use red and blue for binary partitions
		colours = [1 0 0; 0 0 1];
	end
else
	if (nc > 1 & size(colours,1) ~= nc) | (nc==1 & size(colours,1)~=2 & nc==1 & size(colours,1)~=1),
		error('Colours matrix dimensions do not match labels');
	end

	if ischar(colours) & size(colours,2) ~=1,
		error('Character array must have only one column');
	elseif size(colours,2) ~= 3,
		error('RGB colour specification requires three numbers for each colour');
	end
end
out = colours;
