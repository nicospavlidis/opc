function out = fgt_kde(x, y, h)
%Fast Gauss Transform to compute one-dimensional KDE with Gaussian kernels
%OUT = FGT_KDE(X, Y, H)
%
% Input:
%	(x) Column vector representing one dimensional sample
%	(y) Column vector of points at which the value of the KDE will be evaluated
%	(h) Bandwidth parameter


N=size(x,1);
h2 = sqrt(2)*h;

% Gauss transform computes exp{ -||x - m||/h^2 }
%out = figtree(x', sqrt(2)*h, ones(N,1), y', eps)./(N*h*sqrt(2*pi));

[xc , A_k] = fgt_model(x' , ones(1,N) , h2, 30, min(2*sqrt(N),N), 8);
out = fgt_predict(y', xc, A_k, h2, 30)'./(N*h*sqrt(2*pi));

