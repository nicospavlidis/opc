function out = binAssign(labels, clusters)
%Assigns true cluster label (labels) of each observation to the element of the binary partition (clusters = {1,2}^n) that contains the majority of its observations
%OUT = BINASSIGN(LABELS, CLUSTERS)
%
% Inputs: 
%	(LABELS): Vector of true cluster labels
%	(CLUSTERS): Binary assignment of observations
%
% Output:
%	(OUT): Assignment of each observation to cluster containing 
%		majority of observations of its true label

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

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
