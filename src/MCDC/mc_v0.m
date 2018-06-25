function v = mc_v0(X)
%Default projection vector for maximum clusterability projection pursuit
%V = MC_V0(X)
% 
% Input:
%	(X) Data matrix
%
% Output:
%	(V) Vector connecting the centroids of 2-means applied on (X)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

[~,C] = kmeans(X,2,'EmptyAction','singleton','Replicates',1);
v = (C(1,:) - C(2,:))';
v = v./norm(v,2);
