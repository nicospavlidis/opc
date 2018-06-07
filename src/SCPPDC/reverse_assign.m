function L = reverse_assign(coresets2clusters, data2coresets)

% in the case below there has been no micro-clustering 
if length(coresets2clusters)==length(data2coresets),
	L = coresets2clusters;
	return;
end

L = zeros(length(data2coresets),1);

% for each cluster of coresets, do
for i=1:max(coresets2clusters),

	% find which coresets are assigned to this cluster
	index = find(coresets2clusters==i);

	% find observations assigned to each of the coresets in index
	a = [];
	for j=1:length(index),
		a = [a; find(data2coresets==index(j))];
	end
	% assign all such observations to cluster i
	L(a) = i;
end
