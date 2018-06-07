function out = binAssign(labels, clusters)
%Assigns true cluster labels (labels) of each observation to the binary partition (clusters = {1,2}^n) that contains the majority of its observations
%OUT = BINASSIGN(LABELS, CLUSTERS)
%
% Used in process of assessing quality of binary partition of datasets that
% contain multiple clusters
%
% Inputs: 
%	labels: Vector of true cluster labels
%	clusters: Binary assignment of observations
%
% Output:
%	out: Assignment of each observation to cluster containing 
%		majority of observations of its true label

l = unique(labels);
if length(l)==2,
	out = ones(length(labels),1);
	out(labels==l(2))=2;
	return;
end

% if all obs. are assigned to a single cluster prediction = majority class
if length(unique(clusters))==1,
	out = ones(length(labels),1);
	return;
end

% Creates |labels| - by - |clusters| table
% T(i,j) : Number of observations with label j assigned to cluster i
T = mycrosstab(clusters,labels);
[~,assign] = max(T);
out = labels;
for i=1:length(l),
	out( labels==l(i) ) = assign(i);
end
