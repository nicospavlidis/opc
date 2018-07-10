function [df, dkde] = dpkde(x, proj, h, alpha, eta, epsilon)
%Penalised density on hyperplane: Projection index for MDH
%[DF, DKDE] = DPKDE(X, PROJ, H, ALPHA, ETA, EPSILON)
%
% Inputs:
%	(X): location of split along the projected data
%	(PROJ): Univariate projection of dataset (X*v)
%	(H): bandwidth parameter
%	(ALPHA): range over which minimisers of 1D density are sought
%	(ETA): Term in penalty function controlling maximum distance between minimisers of the
%		kde and the penalised density integral (recommended value: 0.01)
%	(EPSILON): Term in penalty function controlling smoothness (recommended value: 1)
%
% Output: 
%	(DF): derivative of one-dimensional penalised density on hyperplane w.r.t. split point (x)
%	(DKDE): derivative of one-dimensional kernel density estimator w.r.t. split point (x)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

s = std(proj);
L = exp(-0.5)/(sqrt(2*pi)*(h*h)* eta^epsilon);
if length(x) == 1,
	if x < -alpha*s,
		der = -(1+epsilon) * (-alpha*s-x)^epsilon * L;
	elseif x > alpha*s,
		der =  (1+epsilon) * (x-alpha*s)^epsilon * L;
	else
		der = 0;
	end
	%dkde = sum( (proj-x).*normpdf(proj, x, h) )/(h*h*length(proj));
	dkde = dkdeDxC(x,proj,h);
	df = dkde + der;
else
	N = length(proj);
	q = [proj, ones(N,1)];
	out = figtree(proj', sqrt(2)*h, q, x', eps)./(sqrt(2*pi)*N*h^3);
	dkde = out(:,1) - x.*out(:,2);

	% derivative of penalty term
	[~, index] = max([zeros(length(x),1), x-alpha*s, -x-alpha*s], [],2);
	der = zeros(length(x),1);
	der(index==2) =  L *(1+epsilon)*( x(index==2)-alpha*s).^epsilon;
	der(index==3) = -L *(1+epsilon)*(-x(index==3)-alpha*s).^epsilon;

	df = dkde + der;
end
