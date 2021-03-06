function out = df_md(v,X,pars,bmin)
%Derivative of penalised density on hyperplane criterion
%OUT = DF_MD(V,X,PARS,BMIN)
%
% Inputs:
%	(V): projection vector
%	(X): N-by-D Data matrix
%	(PARS): Structure array containing parameter settings for MDH (alpha, eta, epsilon, h)
%	(BMIN): (optional argument) location of optimal split on (v)
%
% Output:
%	(OUT): derivative of penalised density on hyperplane criterion

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if nargin<4,
	[~,bmin] = f_md(v,X, pars);
elseif isempty(bmin),
	[~,bmin] = f_md(v,X, pars);
end

alpha = pars.alpha;
eta = pars.eta;
epsilon = pars.epsilon;

h = pars.bandwidth;
if isa(h,'function_handle'),
	h = pars.bandwidth(X,pars);
	assert(numel(h)==1 & h>0,'f_md: Error incorrect bandwidth');
end

normV = norm(v,2);

proj = X*(v./normV);

N = size(proj,1);
%out = df_dp(x, N, p, h)';
% derivative w.r.t. to centers of 1d kde
out = normpdf(proj,bmin,h)'.*(bmin-proj)'/(h*h*N);

% adding penalty term
s= std(proj);
if bmin > alpha*s,
	C = (1+epsilon)*exp(-0.5)/(sqrt(2*pi)*(h*h)*(eta^epsilon));
	out = out + C*(bmin - alpha*s)^epsilon * (-alpha)/((N-1)*s) * proj';
elseif bmin < -alpha*s,
	C = (1+epsilon)*exp(-0.5)/(sqrt(2*pi)*(h*h)*(eta^epsilon));
	out = out + C*(-bmin - alpha*s)^epsilon * (-alpha)/((N-1)*s) * proj';
end

normV3 = normV^3;
%DwV = diag(ones(length(v),1)./normV) - v*(v'./normV3);
%out = out * X * DwV;
out = out * X;
out = out./normV - (out*v)*(v'./normV3);
