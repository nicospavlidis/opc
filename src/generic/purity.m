function out = purity(clusters, trueLabels)
%Computes purity (misclassification rate) of cluster assignment
%OUT = PURITY(CLUSTERS, TRUELABELS)
%
% Inputs:
%	(CLUSTERS): Estimated cluster assignment
%	(TRUELABELS): True cluster assignment

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

T = mycrosstab(trueLabels, clusters);

[maxCluster, clus2class] = max(T,[],1);
%clusterNum = sum(T);

out = sum(maxCluster)/length(trueLabels);

