function v = mc_v0(X)
%Default projection vector for maximum clusterability projection pursuit
%V = MC_V0(X)
% 
% Input:
%	(X) Data matrix
%
% Output:
%	(v) Vector connecting the centroids of 2-means applied on (X)

[~,C] = kmeans(X,2,'EmptyAction','singleton','Replicates',1);
v = (C(1,:) - C(2,:))';
v = v./norm(v,2);
