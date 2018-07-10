function out = md_reldepth(v, X, pars)
%Relative Depth
%OUT = RELDEPTH(V, X, PARS)
% 
% Inputs:
%	(V): Projection vector
%	(X): Data matrix (if X is a column vector then it is assumed that 
%		it is the projected dataset: Data*v)
%	(PARS): Structure array that contains (bandwidth) and optionally range (alpha)
%
% Output:
%	(OUT): Relative depth 

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if ~isfield(pars,'bandwidth'),
	error('Bandwidth parameter is undefined'),

elseif isa(pars.bandwidth,'function_handle') & size(X,2)>1,
	h = pars.bandwidth(X,pars);

elseif numel(pars.bandwidth)==1,
	if pars.bandwidth<=0
		error('Bandwidth parameter must be positive scalar'),
	else 
		h = pars.bandwidth;
	end
else 
	error('Incorrect specification of bandwidth parameter');
end

if ~isfield(pars,'alpha'), pars.alpha = []; end
if ~isfield(pars,'minsize'), pars.minsize = 1; end

[N,dim] = size(X);
if dim > 1,
	X = X*(v./norm(v,2));
end

s = std(X);
y = linspace(min(X), max(X), 200)';
f = fgt_kde(X,y,h);

% Determine range over which modes and minimisers will be sought
range = [];
if ~isempty(pars.alpha),
	if pars.alpha > 0,
		range = find(y>=(-pars.alpha*s) & y<=pars.alpha*s);
		[fmin,split] = min(f(range));
		ymin = y(split+range(1)-1);
	else
		fmin = kdeC(0,X,h);
		ymin = 0;
		range = 1;
	end
end

% If an alpha value is not specified then we search for the lowest
% minimiser in the range between the first and last modes of the KDE
if isempty(range), 
	sorted = sort(X);
	df = diff(f);
	modes = find(df(1:end-1)>=0 & df(2:end)<0) + 1;

	if length(modes) == 1,
		out = 0;
		return;
	end
	% Finding lowest minimiser subject to minimum size constraint
	%	index of y's that satisfy minsize constraint
	index = find(y>sorted(pars.minsize) & y<sorted(N-pars.minsize+1));
	%	find minimisers within y-range 
	df = diff(f(index));
	%	location of antimodes in the original f,y vectors
	antimodes = find(df(1:end-1)<0 & df(2:end)>0) + index(1);

	if isempty(antimodes),
		out = 0;
		return;
	end
	%	lowest local minimum
	[fmin, id] = min(f(antimodes));
	ymin = y(antimodes(id));
end
% increase by tiny amount to avoid division by zero below
fmin = fmin + 1.0e-10;

L = max( max(f(y<ymin)), fmin);
R = max( max(f(y>ymin)), fmin);
out = min(L,R)/fmin - 1;
if out < 0,
	warning('md_reldepth has produced negative output')
	keyboard;
end
%assert(out>=0);
