function out = fgt_kde(proj, Y, h)
%function OUT = FGT_KDE(X, Y, H)
%
% Improved Gauss Transform to compute one-dimensional KDE with Gaussian kernels
% Input:
%	(x) Column vector representing one dimensional sample
%	(y) Points at which the value of the KDE will be evaluated
%	(h) Bandwidth parameter
% Returns:
%	(out): Estimated density at (y)


N=size(proj,1);
% Gauss transform computes exp{ -||x - m||/h^2 }
out = figtree(proj', sqrt(2)*h, ones(N,1), Y', eps)./(N*h*sqrt(2*pi));
