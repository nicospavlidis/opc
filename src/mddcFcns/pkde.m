function [f, fkde] = pkde(x, proj, h, alpha, eta, epsilon)
%Penalised density on of one-dimensional (projected) dataset at (x)
%[F, FKDE] = PKDE(X, PROJ, H, ALPHA, ETA, EPSILON)
%
% Inputs:
%	(X): location of split along the projected data
%	(PROJ): Univariate data projection
%	(H): bandwidth parameter
%	(ALPHA): range over which minimisers of 1D density are sought
%	(ETA): Term in penalty function controlling maximum distance between minimisers of the
%		kde and the penalised density integral (recommended value: 0.01)
%	(EPSILON): Term in penalty function controlling smoothness (recommended value: 1)
%
% Output: 
%	(F): Penalised density on hyperplane at point (X)
%	(FKDE): Density on hyperplane at (X)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

s = std(proj);
L = exp(-0.5)/(sqrt(2*pi)*(h*h)* eta^epsilon);

if length(x)==1,
	fkde = kdeC(x,proj,h);
else
	fkde = fgt_kde(proj,x,h);
end
f = fkde + L* max([zeros(length(x),1), x-alpha*s, -x-alpha*s], [],2).^(1+epsilon);
