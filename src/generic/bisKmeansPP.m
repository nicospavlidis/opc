function [v,fval,idx] = bisKmeansPP(X)
	[idx,C,sumd] = kmeans(X,2,'EmptyAction','singleton','Replicates',1);
	v = (C(1,:) - C(2,:))';
	if sum(abs(v))<sqrt(eps),
		v = pcacomp(X,1);
		fval = -inf;
		return;
	end
	v = v./norm(v,2);

	% Recommended in Steinbach, Karypis and Kumar. A Comparison of Document Clustering Techniques. 2000
	fval = size(X,1);
end

