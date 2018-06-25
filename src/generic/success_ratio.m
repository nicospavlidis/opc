function out = success_ratio(clusters, trueLabels)
%Quantifies accuracy with which a binary partition splits data containing two or more clusters
%OUT = SUCCESS_RATIO(CLUSTERS, TRUELABELS)
%
% Inputs:
%	(CLUSTERS): Estimated binary partition
%	(TRUELABELS): True cluster assignment
%
%Reference:
%Pavlidis N.G., Hofmeyr D.P., Tasoulis S.K. (2016) Minimum Density Hyperplanes. 
%Journal of Machine Learning Research, 17(156), 1–33.

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

T = mycrosstab(clusters,trueLabels);
[~,assign] = max(T);

if length(unique(assign))==1,
	out = 0;
	return;
end

suc = min( sum(T(1,assign==1)), sum(T(2,assign==2)) );
er  = sum( min(T,[],1) );

out = suc/(suc+er);

