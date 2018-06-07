function idx = tree2clusters(T)
%Assigns cluster labels from a cluster hierarchy (ctree object) 
%IDX = TREE2CLUSTERS(T)
%
% Input:
%	T: ctree object
%
% Output:
%	IDX: Cluster assignment vector

% first node is always the root
idx = zeros(size(T.data,1),1);

index = T.findleaves;
for i=1:length(index),
	idx( abs(T.Node{index(i)}.idx) ) = i;
end
assert(min(idx)>0,'tree2clusters: Error -- Observations unassigned to clusters');
