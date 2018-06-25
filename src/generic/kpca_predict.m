function predict = kpca_predict(K,Knew,pcv)
%Returns projections of feature vectors giving rise to kernel matrix (Knew) onto principal components vectors (pcv) of the kernel matrix (K)
%PREDICT = KPCA_PREDICT(K,KNEW,PCV)
%
% Inputs:
%	(K): Kernel matrix on which principal component vectors (PCV) were estimated through KPCA function
%	(KNEW): Kernel matrix of inner products (in feature space) between observations that 
%		are to be predicted and observations that gave rise to (K)
%	(PCV): Principal components vectors of (K) as estimated through KPCA
% Output:
%	(PREDICT): Projection of feature vectors in (KNEW) onto (PCV)
%
%Reference:
%B. Scholkopf, A. Smola, K.-R. Mueller. Nonlinear component analysis 
%as a kernel eigenvalue problem. Neural Computation 10:1299-1319, 1998.


%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

n = size(Knew,1);
l = size(K,1);
Knew = Knew - ones(n,l)*K/l - Knew*ones(l,l)/l + ones(n,l)*K*ones(l,l)/(l^2);
predict = Knew * pcv;
