function [v,fval,idx] = bisKmeansPP(X)
	[idx,C,sumd] = kmeans(X,2,'EmptyAction','singleton','Replicates',1);
	v = (C(1,:) - C(2,:))';
	v = v./norm(v,2);;

	% Recommended in Steinbach, Karypis and Kumar. A Comparison of Document Clustering Techniques. 2000
	fval = size(X,1);
end

