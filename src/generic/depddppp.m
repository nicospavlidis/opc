function [v, fval, idx] = depddppp(X, pars)
%Projection Pursuit function for dePDDP algorithm
%[V, FVAL, IDX] = DEPDDPPP(X, PARS)
%
% Inputs:
%	(X): N-by-D data matrix
%	(PARS): Parameter structure containing algorithm parameters
%
% Outputs:
%	(V): Projection vector: 1st Principal Component of X
%	(FVAL): Function of KDE at split point
%	(IDX): Binary cluster assignment


[N,dim] = size(X);

% Default return values
v = pcacomp(X,1);
fval = -inf;
idx = ones(N,1);

if dim == 1,
	fprintf('Effective dimensionality of 1: No point in projection pursuit\n');
	return;
end

if isfield(pars,'minsize'),
	if 2*pars.minsize > N,
		return;
	end
	minsize = pars.minsize;
else
	minsize = 1;
end

% Set bandwidth
pars.h = pars.bandwidth;
if isa(pars.h,'function_handle'),
	pars.h = pars.h(X,pars);
end

% no bound on range over which minimisers are sought
pars.alpha = [];

% Project onto 1st principal component
proj = X*v;

sorted = sort(proj);
y = linspace(sorted(1), sorted(end), 500)';
f = fgt_kde(proj,y,pars.h);

df = diff(f);
% location of modes on estimated kde
modes = find(df(1:end-1)>=0 & df(2:end)<0) + 1;
% sanity check
if isempty(modes),
	error('MATLAB:dePDDPpp:No modes on 1D kde: This should not have occured');
end

% if KDE is unimodal cluster cannot be split
if length(modes) == 1,
	return;
end

% Finding lowest minimiser subject to minimum size constraint
%	index of y's that satisfy minsize constraint
index = find(y>sorted(minsize) & y<sorted(N-minsize+1));
%	find minimisers within y-range 
df = diff(f(index));
%	location of antimodes in the original f,y vectors
antimodes = find(df(1:end-1)<0 & df(2:end)>0) + index(1);

if isempty(antimodes),
	%warning('No valid hyperplane separator for minimum cluster size: %i\n',minsize);
	return;
end

%	lowest local minimum
[fval, id] = min(f(antimodes));
%	split point
b = y(antimodes(id));

% cluster assignment in {1,2}^N
idx = (proj > b) + 1;

% Since by default we want to split cluster with lowest minimum value of KDE at split point
fval = -fval; 
