function clEval = cluster_performance(clusters, trueLabels)
%Returns structure containing: Purity, Adjusted Rand Index, Normalised Mutual Information and V-measure
%CLEVAL = CLUSTER_PERFORMANCE(CLUSTERS, TRUELABELS)
%
% Inputs:
%	clusters: Estimated cluster assignment
%	trueLabels: True cluster assignment

N = length(trueLabels);
% Confusion matrix Actual class is in rows
T = mycrosstab(trueLabels, clusters);

maxCluster = max(T,[],1);

clEval.Purity = sum(maxCluster)/N;


ProbClust = sum(T)/N;
ProbClass = sum(T,2)/N;
I=0;
for i = 1:size(T,1),
	for j=1:size(T,2),
		if T(i,j)>0,
			I = I + (T(i,j)/N) * log((T(i,j)/N)/(ProbClass(i)*ProbClust(j)));
		end
	end
end
HClass = -sum(ProbClass.*log(ProbClass));
HClust = -sum(ProbClust.*log(ProbClust));
clEval.NMI = I/sqrt(HClass*HClust);


nis=sum(sum(T,2).^2);
njs=sum(sum(T,1).^2);
A = nchoosek(N,2) + sum(sum(T.^2)) - 0.5*(nis+njs);
nc = (N*(N^2+1)-(N+1)*nis-(N+1)*njs + 2*(nis*njs)/N)/(2*(N-1));
clEval.AdjRand = (A-nc)/(nchoosek(N,2) - nc);



h=0;
rowSum = sum(T,1);
for j=1:size(T,2),
	for i=1:size(T,1),
		if T(i,j) > 0,
			h = h + (T(i,j)/N)*log(T(i,j)/rowSum(j));
		end
	end
end
H_C_K=-1.0*h;


if abs(h)<eps,
	hom = 1;
else
	hom = 1.0 - (H_C_K/HClass);
end

h=0;
colSum = sum(T,2);
for i = 1:size(T,1),
	for j = 1:size(T,2),
		if T(i,j) > 0,
			h=h+ (T(i,j)/N) * log(T(i,j)/colSum(i));
		end
	end
end
H_K_C=-1.0*h;

if abs(h)<eps,
	comp = 1;
else
	comp = 1.0-(H_K_C/HClust);
end

clEval.Vmeasure = (2.0*comp*hom)/(comp+hom);
