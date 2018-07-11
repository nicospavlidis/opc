setpath;

% Load optiditigs dataset: X := data matrix, labels := true cluster assignment                                     
load('datasets/optidigits.mat');
                                                                                                                   
% rng() function is implemented in MATLAB only
rng(987)
% Standard clustering 
km = kmeans(X,10);
sc = spclust(X,10);

fprintf('Performance of k-means\n');
cluster_performance(km,labels)
fprintf('Performance of normalised spectral clustering, Ng et al. (2001)\n');
cluster_performance(sc,labels)

% Common approach of performing dimensionality reduction
% as a pre-processing step prior to clustering
V = pcacomp(X,[1:9]);

kmPCA = kmeans(X*V,10);
fprintf('k-means on data projected onto first 9 PCs\n');
cluster_performance(kmPCA,labels)

scPCA = spclust(X*V,10);
fprintf('normalised spectral clustering on data projected onto first 9 PCs\n');
cluster_performance(scPCA,labels)


% Project onto the first 33 PCs to retain (slightly
% more than) 90% of total data variability
V = pcacomp(X,[1:33]);

kmPCA = kmeans(X*V,10);
fprintf('k-means on data projected onto first 33 PCs\n');
cluster_performance(kmPCA,labels)

scPCA = spclust(X*V,10);
fprintf('normalised spectral clustering on data projected onto first 33 PCs\n');
cluster_performance(scPCA,labels)


%% LDA-k-means

% Visualisation of clustering result at each iteration
% Output illustrates that the projection index is not 
% improving monotonically over the execution of the algorithm
[ldakm, U, fval] = ldakmeans(X,10,'verb',1);
fprintf('Performance of LDA-k-means\n');
cluster_performance(ldakm,labels)

%%%% Visualisation using the actual labels
[ldakm, U, fval] = ldakmeans(X,10,'labels',labels,'verb',1);
fprintf('Performance of LDA-k-means in second execution\n');
cluster_performance(ldakm,labels)


%% PDDP 
[idx1,t1] = pddp(X,10);

% Evaluate cluster performance
fprintf('Performance of PDDP\n');
cluster_performance(idx1, labels)

% Visualise cluster hierarchy
%plot(t1);
%print('-f2', 'documentation/figures/pddp1', '-dpng','-r0');
%
%% Visualise cluster hierarchy with label information
%plot(t1,labels);
%print('-f2', 'documentation/figures/pddp2', '-dpng','-r0');
%
%nplot(t1,1);
%print('-f1', 'documentation/figures/pddp3', '-dpng','-r0');
%
%nplot(t1,1,labels);
%print('-f1', 'documentation/figures/pddp4', '-dpng','-r0');

% Get root node of cluster hierarchy
node1 = t1.get(1); 

% Assess quality of binary partition through success ratio
fprintf('Success ratio at root node of PDDP\n');
success_ratio(sign(node1.idx), labels)


%%% DEPDDP
[idx2, t2] = depddp(X,10); 
% Plot cluster hierarchy
plot(t2);
%print('-f2', 'documentation/figures/depddp1', '-dpng','-r0');

% Plot split at root node
nplot(t2,1);
%print('-f1', 'documentation/figures/depddp2', '-dpng','-r0');

% Set minimum size
[idx2, t2] = depddp(X,10,'minsize',10); 

% Set bandwidth parameter
fh = @(x,p)(0.45*size(x,1)^(-0.2)*std(x* pcacomp(x,1)));
[idx2, t2] = depddp(X,10,'minsize',10,'bandwidth',fh); 


%%% MCDC
[idx3, t3] = mcdc(X,10);
fprintf('Performance of MCDC algorithm\n');
cluster_performance(idx3,labels)

%plot(t3)
%print('-f2', 'documentation/figures/mcdc1', '-dpng','-r0');
%nplot(t3,1);
%print('-f1', 'documentation/figures/mcdc2', '-dpng','-r0');

[idx4, t4] = mcdc(X,10,'v0', @(x,p)(pcacomp(x,1)) );
fprintf('Performance of MCDC algorithm using 1st PC as initialisation\n');
cluster_performance(idx4,labels)

plot(t4)
%print('-f2', 'documentation/figures/mcdc3', '-dpng','-r0');



%% NCUTDC ALGORITHM 
[idx5, t5] = ncutdc(X,10);
fprintf('Performance of NCUTDC\n');
cluster_performance(idx5,labels)

% MDDC
% More detailed illustration of this algorithm in next section
[idx6, t6] = mddc(X,10);
fprintf('Performance of MDDC\n');
cluster_performance(idx6,labels)
nplot(t6,2);
%print('-f1', 'documentation/figures/mddc1', '-dpng','-r0');


% SPECTRAL CLUSTERING: DRSC
% Step 1: micro-clustering
rng(56789);
[d2c,centroids] = kmeans(X,200);

% Step 2: Set scale parameter to default value used by SCPP
s = scpp_def_sigma(X); 

% Apply DRSC
[idx7,W,fval] = drsc(centroids,10,s,'verb',1);

% Assign original observations to clusters
idx7 = reverse_assign(idx7, d2c);

fprintf('Performance of DRSC\n');
cluster_performance(idx7,labels);

%SPECTRAL CLUSTERING: SCPPDC
rng(1098765);
[idx8,t8] = scppdc(X,10);
plot(t8)
%print('-f2', 'documentation/figures/scpp1', '-dpng','-r0');

plot(t8,labels)
%print('-f2', 'documentation/figures/scpp2', '-dpng','-r0');

fprintf('Performance of SCPP\n');
cluster_performance(idx8,labels);


%%%% MODEL VALIDATION AND MODIFICATION
rng(201800630);
[idx,t] = mddc(X,5)

plot(t)
%print('-f2', 'documentation/figures/val1', '-dpng','-r0');

nplot(t,2)
%print('-f1', 'documentation/figures/val2', '-dpng','-r0');

nplot(t,3)
%print('-f1', 'documentation/figures/val3', '-dpng','-r0');

t1 = split(t,2)
t1 = split(t1,10)

plot(t1)
%print('-f2', 'documentation/figures/val4', '-dpng','-r0');

t1 = prune(t1,3);
plot(t1)
%print('-f2', 'documentation/figures/val5', '-dpng','-r0');


% Consider alternative parameter settings for PP algorithm:
% initialise at 2nd PC 
t2 = split(t1,3,'v0',@(x,p)(pcacomp(x,2)), 'verb',1)

nplot(t2,3)
%print('-f1', 'documentation/figures/val6', '-dpng','-r0');

% initialise at 3rd PC, and increase range of admissible MDHs
split(t1,3,'v0',@(x,p)(pcacomp(x,3)), 'alphamax',1.2, 'verb',1)

% Use the MDH obtained after initialisation at the 3rd PC
t1 = split(t1,3,'v0',@(x,p)(pcacomp(x,3)), 'alphamax',1.2)
nplot(t1,3)
%print('-f1', 'documentation/figures/val7', '-dpng','-r0');

% Visualise resulting clusters and split those that
% appear to contain more than one clusters
nplot(t1,9);
t1 = split(t1,9);

nplot(t1,11);
t1 = split(t1,11);

nplot(t1,12);
t1 = split(t1,12);

nplot(t1,15);
t1 = split(t1,15);

nplot(t1,10);
t1 = split(t1,10);

idx = tree2clusters(t1);
fprintf('Performance of final clustering model\n');
cluster_performance(idx, labels)

% Final clustering model
plot(t1,labels);
%print('-f2', 'documentation/figures/val8', '-dpng','-r0');



% EXTENSIONS
% MAXIMUM MARGIN
load('datasets/optidigitsTest.mat');
index = find(labels==3 | labels==9);
X = normalise(X(index,:),1);
labels = labels(index);

%%% Estimate MDH with default bandwidth
[idx,hp0] = mdh(X);
%er = 1 - purity(idx, labels)

plot(hp0,X);
%print('-f1', 'documentation/figures/mm1', '-dpng','-r0');

hp = hp0;
a = hp.params.alpha;
v0 = 0*hp.v;
while abs(hp.v'*v0) < 1-1.e-10,
	v0 = hp.v;
	h = 0.5*hp.params.bandwidth;
	[idx,hp1]= mdh(X,'v0',v0,'alphamin',a,'alphamax',a,'bandwidth',h);
	if isinf(hp1.fval),
		break;
	else
		hp = hp1;
	end
	%er = 1 - purity(idx, labels)
	%plot(hp,X,labels);
end
plot(hp,X);
%print('-f1', 'documentation/figures/mm2', '-dpng','-r0');

fprintf('Misclassification error of large margin hyperplane obtained through MDH\n');
er = 1 - purity(idx, labels)

%%% Maximum margin example using NCUTH
load('datasets/optidigitsTest.mat');
index = find(labels==3 | labels==9);
X = normalise(X(index,:),1);
labels = labels(index);

%%%% Estimate NCutH with default scaling parameter
[idx,hp0] = ncuth(X);

plot(hp0,X);
%print('-f1', 'documentation/figures/mm3', '-dpng','-r0');

hp = hp0;
v0 = 0*hp.v;
while abs(hp.v'*v0) < 1-1.e-10,
	v0 = hp.v;
	s = 0.5*hp.params.sigma;
	[id,hp1]= ncuth(X,'v0',v0,'sigma',s);

	% Numerical problems can occur for very low scaling parameter:
	% isinf(fval) signals projection pursuit failed 
	if isinf(hp1.fval),
	        break;
	else
	        hp = hp1;
	end
	%er = 1 - purity(id, labels)
end
plot(hp,X);
%print('-f1', 'documentation/figures/mm4', '-dpng','-r0');
fprintf('Misclassification error of large margin hyperplane obtained through MDH\n');
er = 1 - purity(id, labels)


%% Kernel PCA and Nonlinear Clustering

load('datasets/halfmoons.mat');
% Select random subset of 200 observations to estimate KPCA
% (only used to illustrate kpca_predict function)
s = randperm(size(X,1),200);

sigma = 6;
K = exp(-sigma*squareform(pdist(X(s,:)).^2));
pcv = kpca(K);
Knew = exp(-(mypdist2(X, X(s,:)).^2)*sigma);
X2 = kpca_predict(K,Knew,pcv);

% Estimate optimal bi-partition based on ncuth() on the original space
idx1 = ncuth(X);
% Estimate optimal bi-partition based on ncuth() on the kernel defined feature space
idx2 = ncuth(X2);

n = length(idx1);
% Visualise the two partitions
cFig = figure(1);
set(cFig, 'Position', [0 0 1024 1024]);

M = sparse(1:n, idx1, 1);
scatter(X(:,1), X(:,2),[],M*[1 0 0; 0 0 1]);
%print('-f1', 'documentation/figures/kpca1', '-dpng','-r0');

clf;
M = sparse(1:n, idx2, 1);
scatter(X(:,1), X(:,2),[],M*[1 0 0; 0 0 1]);
%print('-f1', 'documentation/figures/kpca2', '-dpng','-r0');


% EXTENDING OPC
load('datasets/optidigits.mat');
[idx,t] = gppdc(X,10, @(x,p)(bisKmeansPP(x)), 'split_index', @(v,x,p)(total_scatter(x)));
