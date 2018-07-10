function [f, Df] = f_df_sc(v,X, pars)
%Function value and derivative of second smallest eigenvalue of normalised Laplacian
%[F,DF] = F_DF_SC(V, X, PARS)
%
% Inputs:
% 	(V): Projection vector
%	(X): N-by-D Data matrix
%	(PARS): parameter struct containing 
%		(sigma) Scaling parameter for Gaussian kernel
%		(minsize): minimum cluster size
%		(beta) (delta): parameters for similarity transformation
%
% Output: 
%	(F): Second smallest eigenvalue of Normalised Laplacian
%       (DF): derivative of (f) w.r.t. projection matrix (v)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

a = pars.beta;
d = pars.delta;
weights = pars.weights;
scale = pars.sigma;

dim = size(X,2);  % dimensionality of original space
pdim = length(v)/dim; % dimensionality of the projection space

% columns correspond to projection vectors
V0 = reshape(v,[dim,pdim]);
% normalise to unit-length
V = bsxfun(@rdivide, V0, sqrt(sum(V0.^2,1)));

% projections
proj = X*V;
% transformed projections Eq.(19)
%%%%%%%%%%%%%%%%% HERE
%p = proj; 
p = sim_transform(proj, a, d, weights);
%prDist = pdist(p);

% Similarity matrix (weighted by microclustering)
W = exp(-(squareform(pdist(p)).^2)./(2*scale^2));
if ~isempty(weights)
	W = (pars.weights*pars.weights') .* W;
end
% D^{-1/2}
u1 = sqrt( sum(W,2) );

% Normalised Laplacian
%L = eye(size(W,1)) - diag(1./u1)*W*diag(1./u1);
%v1 = u1./sqrt(sum(u1.^2));
% Compute second smallest eigenvalue of normalised Laplacian matrix 
% by deflating smallest eigenvalue (0)
% Clear warning messages
%warning('');
%[U,D] = eigs(L + v1*v1',2, 'sm');
%[warnMsg, warnId] = lastwarn;
% Catch warnings
%if ~isempty(warnMsg), keyboard end
%[lambda, index] = sort(diag(D));

% it is numerically more stable to work with A defined below
% Unit length eigenvector of normalised Laplacian (L) and (A)
v1 = u1./sqrt(sum(u1.^2));
%A =  ((1./u1)*(1./u1)').*W;
[U,D] = eigs(((1./u1)*(1./u1)').*W - v1*v1',2, 'lm');

% eigenvales of L = 1 - eigenvalues of (A)
% therefore largest eigenvalue of (A) is smallest eigenvalue of (L)
[lambda,index] = sort(diag(D), 'descend');
lambda = 1 - lambda;
U = U(:,index(1));


% function value: second smallest eigenvalue of L
f = lambda(1);
% u := eigenvector corresponding to second smallest eigenvalue
u = U(:,index(1));
% Differentiability check
if lambda(2) - lambda(1) < sqrt(eps), 
	warning(sprintf('Repeated eigenvalue. eigengap: %d; f(x): %e\n',lambda(2) - lambda(1), f)); 
end

% penalty term
if pdim > 1,
	f = f + pars.omega * sum(sum(triu(V'*V,1).^2));
end

if nargout == 1,
	return;
end

% coeffs matrix is used in Eq. (13) and (14)
u = u./sqrt( sum(W,2) ); % D = sum(W,2) => u = D^(-1/2)*u
distU = squareform(pdist(u));
sumU = bsxfun(@plus, repmat(u.^2,1,length(u))', u.^2);
coeffs = (distU.^2 - f*sumU).* W/(scale^2);

% Derivative w.r.t. each column of V
Df = zeros(pdim,dim);
for i=1:pdim,

	%for j=1:size(p,1), dlp(j) = sum( coeffs(j,:).*(p(:,i)'-p(j,i)) ); end
	%dlp == sum( coeffs .* bsxfun(@minus,p(:,i),p(:,i)'), 1)
	dLambda_dp = sum(coeffs .* bsxfun(@minus,p(:,i),p(:,i)'));

	% Derivative of v=w/norm(w,2) w.r.t. w
	DwV = diag(ones(dim,1)./norm(V0(:,i),2)) - V0(:,i)*(V0(:,i)'./(norm(V0(:,i),2)^3));
	% Chain rule
	%%%%%%%%%%%%%%%%% HERE
	%Df(i,:) = dLambda_dp * X* DwV;
	Df(i,:) = dLambda_dp * dP_dv(proj(:,i), X, weights, a, d) * DwV;

	if pdim > 1,
		noti = [1:i-1,i+1:pdim];
		% (1) V(:,noti)'*V(:,i) : <v_i,v_j> for all j ~= i
		% (2) 2*diag(V(:,noti)'*V(:,i))* V(:,noti)' : multiply jth proj. vector with 2*<v_i, v_j> for all j~=i
		% (3) Rows of 2*diag(V(:,noti)'*V(:,i)) * V(:,noti)' * DwV correspond to
		% d/dw{ <v_i,v_j>^2 } = 2*<v_i,v_j>* v_j' * d/dw{v_i}
		% (4) Sum along columns to obtain derivative of \sum_{j \neq i} <v_i,v_j>^2 w.r.t. w_i
		Df(i,:)=Df(i,:) + pars.omega*2*sum(diag(V(:,noti)'*V(:,i)) * V(:,noti)' *DwV, 1);
	end
end
Df = reshape(Df',[dim*pdim,1]);

% Finite Differences: to remove similarity transformation uncomment (HERE) 
%pert=1.0e-4;
%finDiff = zeros(length(v),1);
%for i=1:length(v),
%	vn = v;
%	vn(i) = v(i) + pert;
%	f1 = f_sc(vn,X,pars);
%
%	vn(i) = v(i) - pert;
%	f2 = f_sc(vn,X,pars);
%
%	finDiff(i) = (f1-f2)/(2*pert);
%end
%max(abs(Df - finDiff))
%%keyboard

end

