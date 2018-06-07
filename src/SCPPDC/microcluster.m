function [centers, weights, data2coresets] = microcluster(X,K)
%
% estimates micro-clusters

if size(X,1) < K,
	centers = X;
	data2coresets = [1:size(X,1)]';
	weights = ones(size(X,1),1);
	return;
end

[data2coresets, centers] = kmeans(X, K,'EmptyAction','singleton','Replicates',1);

% defining thus the number of clusters is safer as potentially the maximum
% number of clusters may not be obtained
num_clust = size(centers,1);
weights = zeros(num_clust,1);
for i=1:num_clust,
	weights(i) = length(find(data2coresets==i));
end

if sum(weights)~=size(X,1) | sum(weights==0),
	error('microcluster','Weight vector incorrectly computed');
end

end
