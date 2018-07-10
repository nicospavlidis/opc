function L = fixLabels(labels)
%Enforces cluster labels to be in the range 1:K
%L = FIXLABELS(LABELS)
%
% Input:
%	(LABELS): Cluster assignment
%
% Output:
%	(L): New label vector \in {1,...,K}^N

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

l = unique(labels);	
L = zeros(length(labels),1);
for i=1:length(l),
	L( labels==l(i) ) = i;
end
