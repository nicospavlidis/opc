function [v,fval,idx] = bisKmeansPP(X)
%Projection pursuit function template to implement Bisecting K-Means
%[V,FVAL,IDX] = BISKMEANSPP(X)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

	if size(X,1)<2,
		idx = ones(1,1);
		v = ones(size(X,2),1)./sqrt(size(X,2));
		fval = -inf;
		return;
	end

	[idx,C,sumd] = kmeans(X,2,'EmptyAction','singleton','Replicates',1);

	v = (C(1,:) - C(2,:))';
	if sum(abs(v)) < sqrt(eps),
		v = pcacomp(X,1);
		fval = -inf;
		return;
	end
	v = v./norm(v,2);

	% Recommended in Steinbach, Karypis and Kumar. A Comparison of Document Clustering Techniques. 2000
	fval = size(X,1);
end

