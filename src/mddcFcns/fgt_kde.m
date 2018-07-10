function out = fgt_kde(proj, Y, h)
%Interface to improved Fast Gauss Transform to compute one-dimensional KDE with Gaussian kernels
%OUT = FGT_KDE(X, Y, H)
%
% Input:
%	(X): Column vector representing one dimensional sample
%	(Y): Points at which the value of the KDE will be evaluated
%	(H): Bandwidth parameter
%
% Output:
%	(OUT): Estimated density at (Y)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

N=size(proj,1);
% Gauss transform computes exp{ -||x - m||/h^2 }
out = figtree(proj', sqrt(2)*h, ones(N,1), Y', eps)./(N*h*sqrt(2*pi));
