function idx = tree2clusters(T)
%Assigns cluster labels from a cluster hierarchy (ctree object) 
%IDX = TREE2CLUSTERS(T)
%
% Output:
%	(IDX): Cluster assignment vector
%
% Input:
%	(T): ctree object

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

% first node is always the root
idx = zeros(size(T.data,1),1);

index = T.findleaves;
for i=1:length(index),
	idx( abs(T.Node{index(i)}.idx) ) = i;
end
assert(min(idx)>0,'tree2clusters: Error -- Observations unassigned to clusters');
