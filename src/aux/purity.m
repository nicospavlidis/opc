function out = purity(clusters, trueLabels)
%Computes purity (misclassification rate) of cluster assignment
%OUT = PURITY(CLUSTERS, TRUELABELS)
%
% Inputs:
%	clusters: Estimated cluster assignment
%	trueLabels: True cluster assignment

T = mycrosstab(trueLabels, clusters);

[maxCluster, clus2class] = max(T,[],1);
%clusterNum = sum(T);

out = sum(maxCluster)/length(trueLabels);

