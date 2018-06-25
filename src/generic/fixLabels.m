function L = fixLabels(labels)
%Enforces cluster labels to be in the range 1:K
%L = FIXLABELS(LABELS)
%
% Returns:
%	(L): New label vector
%
% Input:
%	(LABELS): Original labels

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

l = unique(labels);	
L = zeros(length(labels),1);
for i=1:length(l),
	L( labels==l(i) ) = i;
end
