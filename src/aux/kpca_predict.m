function predict = kpca_predict(K,Knew,pcv)
%Returns projections of feature vectors giving rise to kernel matrix (Knew) onto principal components vectors (pcv) of the kernel matrix (K)
%PREDICT = KPCA_PREDICT(K,KNEW,PCV)
%
% Inputs:
%	K: Kernel matrix on which principal component vectors (pcv) were estimated through KPCA
%	Knew: Kernel matrix of inner products (in feature space) between observations that 
%		are to predicted and observations that gave rise to (K)
%	pcv: Principal components vectors of (K) as estimated through KPCA
% Output:
%	predict: Projection of feature vectors in (Knew) onto (pcv)
%
% Reference:
% B. Scholkopf, A. Smola, K.-R. Mueller. Nonlinear component analysis 
% as a kernel eigenvalue problem. Neural Computation 10:1299-1319, 1998.

n = size(Knew,1);
l = size(K,1);
Knew = Knew - ones(n,l)*K/l - Knew*ones(l,l)/l + ones(n,l)*K*ones(l,l)/(l^2);
predict = Knew * pcv;
