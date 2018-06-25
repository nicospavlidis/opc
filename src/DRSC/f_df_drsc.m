function [f, df] = f_df_drsc(Data, W, sigma, U, degs)
%Function value and derivative with respect to last projection vector (column) of W
%[F, DF] = F_DF_DRSC(DATA, W, SIGMA, U, DEGS)
%
% Inputs:
%	(Data): N-by-D data matrix
%	(W): Projection matrix
%	(sigma): scale parameter for Gaussian kernel used to 
%		estimate similarity matrix (K)
%	(U): First (K) eigenvectors of D^{-1/2} *K*D^{-1/2}
%	(degs): Vector of degrees of each vertex: D = diag(degs)
%
% Outputs:
%	(f) Function value
%	(df) Derivative

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

% to avoid division by zero
degs = 1./sqrt( max(degs,sqrt(eps)) );

% Function value: Eq.(11) in Niu,Dy and Jordan AISTATS (2011)
%	there is a typo in that equation: d_i, and d_j should be \sqrt(d_i), \sqrt(d_j)
K = exp(-(pdist(Data*W).^2)./ (2*sigma^2));
p = pdist(U, @(Xi,Xj)(Xi * Xj')) .* K;
f = sum(p .* pdist(degs, @(xi,xj)(xi*xj)));

if nargout > 1,
	% Eq.(12) in Niu,Dy and Jordan AISTATS (2011)
	[N, dim] = size(Data);
	Diff = zeros(dim, nchoosek(N,2));
	for i=1:dim,
		Diff(i,:) = pdist(Data, @(Xi,Xj)bsxfun(@minus,Xi(:,i),Xj(:,i)));
	end;

	w = W(:,end);
	p2 = -pdist(degs, @(xi,xj)(xi*xj))./(sigma^2); 
	df = sum( bsxfun(@times, Diff, p .* p2 .* pdist(Data, @(Xi,Xj)(bsxfun(@minus,Xi,Xj)*w))), 2);
end

% The above implementation is a MATLAB optimised version of the following for loops
% which implement Eqs. (11) and (12) in Niu,Dy and Jordan AISTATS (2011)
%N = size(Data,1);
%value = 0;
%gradw = zeros(size(W,1),1);
%for i=1:N-1,
%	for j=(i+1):N,
%
%		D = Data(i,:) - Data(j,:);
%		upd = (U(i,:)*U(j,:)' * degs(i)*degs(j)) * exp(-norm(D*W,2)^2/(2*sigma^2));
%
%		value = value + upd;
%		gradw = gradw - (upd/sigma^2) *(D*W(:,end))*D';
%	end
%end
