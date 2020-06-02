function [idx,U] = kurtclust(X,varargin)
%KURTOSIS PROJECTION PURSUIT CLUSTERING
%
% [IDX,U] = KURTCLUST(X) produces a flat partitioning of the N-by-D data
% matrix (X). The number of clusters is estimated by the algorithm. The algorithm identifies 
% 2(D) unit-length vectors. The first (D) are directions of maximum kurtosis (directions in which 
% outliers are separated out). The second (D) are directions of minimum kurtosis (directions
% in which the marginal distribution is bimodal and hence directions that separate clusters)
%
% [IDX,U] = kurtclust(X,K) returns the cluster assignment, (IDX); and the
% D-by-2D projection matrix (U) used to perform dimensionality reduction; 
% 
% [IDX,U,FVAL,C,ITER] = ldaKmeans(X, K, 'PARAM1',val1, 'PARAM2',val2, ...)
% specifies optional parameters in the form of Name,Value pairs. These are only
% used for visualisation
%
%  'maxit' - Number of iterations to perform for each one-dimensional projection pursuit step
%	(default: 50)
%
%  'ftol' - Tolerance level for projection pursuit algorithm (default: 1.e-5)
%
%  'verb' - Verbosity. Values greater than 0 enable visualisation during execution
%	Enabling this option slows down the algorithm considerably
%	(default: 0)
%
%  'labels' - true cluster labels. Specifying these enables the computation of performance over 
%	successive iterations and a better visualisation of how clusters are split
%
%  'colours' - Matrix containing colour specification for observations in each of the (K) clusters
%
%Reference:
%D. Pena and F. J. Prieto. Cluster Identification Using Projections.
%Journal of the American Statistical Association, 96(456):1433-1445, 2001.

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2019
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

pars = struct();
pars.maxit = 50;
pars.ftol = 1.e-5;
pars.verb = 0;
pars.v0 = @(x,p)(pcacomp(x,1));
pars.labels = [];
pars.colours = [];

% use an arbitrary number of clusters to ensure that it passes the parser function
pars = myparser(X,2,varargin, pars);
n = size(X,1);
d = rank(X);


% Optimisation options
if ~isOctave(),
	options = optimoptions(@fminunc,'Algorithm','quasi-newton','GradObj','on', 'HessUpdate','bfgs', ...
		'MaxIter',pars.maxit, 'OptimalityTolerance', pars.ftol, ...
		'Display',ifelse(pars.verb > 1,'iter','off'));
else
	options = optimset('GradObj','on','MaxIter',pars.maxit,'TolX',pars.ftol,...
		'Display',ifelse(pars.verb>1,'iter','off'));

end
iter = 0;

% Initialise projection matrix
U = zeros(d, 2*d);
% Remove mean from the data
Y = bsxfun(@minus,X,mean(X));
% Compute d directions that maximise kurtosis: Directions in which outliers manifest
objfungrad = @(v)f_df_kurt(v,X,1);
for i = 1:d,
	
	% Update plotting function options
	if pars.verb > 1,
		pars.title = sprintf('Outlier detection phase dim = %i',i);
		outF = @(a,b,c)outFcn(a,b,c,Y,i,pars);
		if ~isOctave(),
			options = optimoptions(options,'OutputFcn', outF);
		else
			options = optimset(options,'OutputFcn', outF);
		end
	elseif pars.verb,
		fprintf('Outlier detection phase dim = %i\n',i);
	end

	% Identify optimal univariate projection
	[v, f, exitFlag, output] = fminunc(objfungrad, pars.v0(Y), options);

	U(:,i) = v./norm(v,2);
	iter = iter + output.iterations;

	% Project data onto null space of projection vector
	Y = Y - (Y * U(:,i))*U(:,i)';
end

% Compute d directions that minimise kurtosis: Directions that favour bimodality
objfungrad = @(v)f_df_kurt(v,X,0);
% Remove mean from the data
Y = bsxfun(@minus,X,mean(X));
for i = 1:d,
	% Update plotting function options
	if pars.verb > 1,
		pars.title = sprintf('Cluster detection phase dim = %i',i);
		outF = @(a,b,c)outFcn(a,b,c,Y,i,pars);
		if ~isOctave(),
			options = optimoptions(options,'OutputFcn', outF);
		else
			options = optimset(options,'OutputFcn', outF);
		end
	elseif pars.verb,
		fprintf('Cluster detection phase dim = %i\n',i);
	end

	% Identify optimal univariate projection
	[v, f, exitFlag, output] = fminunc(objfungrad, pars.v0(Y), options);

	U(:,d+i) = v./norm(v,2);
	iter = iter + output.iterations;

	% Project data onto null space of projection vector
	Y = Y - (Y * U(:,d+i))*U(:,d+i)';
end

keyboard

% Algorithm to identify significant gaps
% Reset data matrix
Y = bsxfun(@minus,X,mean(X));
% Step 1: Project 
Pu = Y*U;
% Step 2: Standardise
Pu = bsxfun(@minus, Pu, mean(Pu));
Pu = bsxfun(@rdivide, Pu, std(Pu)); 

% Sort in ascending order column-wise
Pu = sort(Pu,1);
% Compute cdf of Normal distribution
Pu = normcdf(Pu);
% Compute gaps
W = Pu(2:n,:) - Pu(1:n-1,:);

idx = ones(n,1);
kappa = 1 - 0.1^(1/n) * d^(-10/(3*n));
for j=1:size(W,2),
	index = find( W(:, j) > kappa );
	for i = 1:length(index),
		k = index(i)+1;
		if idx(k-1) == idx(k),
			idx(k:n) = idx(k:n) + 1;
		end
	end
end
% this should not be necessary but just in case
idx = fixLabels(idx);

%Identify the number of clusters in the data
pos = cell(1,max(idx));
for i=1:max(idx),
	pos{i} = find(idx == i);
end
% order clusters in descending size
[ksize, order] = sort(cellfun(@length, pos), 'descend');
% reorder elements in cell array
pos = pos(order);

cutoff = chi2inv(0.99, d);
i = 1;
while i < length(pos),
	% indices of observations in current cluster
	id = pos{i};
	if rank(Y(id,:)) < d,
		break;
	end
	mu = mean( Y(id, :) );
	invS = inv( cov(Y(id, :)) );

	revise = 0;
	j = i+1;
	while j <= length(pos),
		% indices in cluster j
		jd = pos{j};
		delta = diag(bsxfun(@minus, Y(jd,:), mu) * invS * bsxfun(@minus, Y(jd, :), mu)');
		index = find(delta <= cutoff);

		if ~isempty(index),
			id = [id; jd(index)];

			%jd(find( delta < cutoff)) = [];
			pos{j} = setdiff(jd, jd(index));
			if isempty(pos{j}),
				pos = { pos{[1:j-1,j+1:end] } };
				j = j - 1;
			end
			revise = 1;
		end
		j = j+1;
	end

	if ~revise,
		i = i+1;
		[ksize, order] = sort(cellfun(@length, pos), 'descend');
		pos = pos(order);
	else
		pos{i} = sort(id);
	end
end

idx = ones(n,1);
for i = 2:length(pos),
	idx(pos{i}) = i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, df] = f_df_kurt(v,X,maximise)

if nargin<3,
	maximise = false;
end

% Eq. (11) from the paper assumes that the Covariance matrix of X is
% invertible. This is not be necessarily true (and it certainly won't be if
% n<d). In this case we need to differentiate the actual equation for kurtosis 
%
%n = size(X,1);
%
%nv = norm(v,2);
%invS = inv(chol(cov(X)));
%w = invS*(v./nv);
%
%Xw = X*invS*(v./nv);
%% mean(Xw.^4) == kurtosis(Xw,1)*( ((n-1)/n) * var(Xw) )^2,
%f = mean(Xw.^4);
%df = 4* mean(diag((Xw.^3))*X) * (invS/nv - (invS*v)*(v'./nv^3));
%
%if maximise,
%	f = -f;
%else
%	df = -df;
%end

[n,d] = size(X);
nv = norm(v,2);
w = v./nv;
Xw = X*w;

nvar = sum(Xw.^2);
nkurt = sum(Xw.^4);
f = n*nkurt/(nvar*nvar);
if maximise,
	f = -f;
end

if nargout > 1,
	df_numer = 4*sum(diag(Xw.^3) * X)*(eye(d) - w*w')/nv;
	df_denom = -2*nvar^(-3)*2*sum(diag(Xw)*X)*(eye(d) - w*w')/nv;

	df = n*(df_numer/(nvar*nvar) - df_denom*nkurt);
	if maximise,
		df = -df;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stop = outFcn(x, optimValues, state, X, column, pars, labels, colours)
stop = false;
if isOctave(),
	iter = optimValues.iter;
else
	iter = optimValues.iteration;
end

v = x./norm(x,2);
% direction of maximum variance perpendicular to v
w = pcacomp(X - (X*v)*v', 1);
X = X*[v, w];

hFig = figure(1);
clf;
set(hFig, 'Position', [0 0 1600 800]);
hold on;
if max(pars.labels) > 1,
	l = unique(pars.labels);
	for i=1:length(l),
		s = find(pars.labels==l(i));
		plot(X(s,1), X(s,2),'bo','MarkerSize',3, 'Color', pars.colours(l(i),:))
		% , 'MarkerFaceColor', colours(l(i),:));
	end
else
	plot(X(:,1), X(:,2),'bo','MarkerSize',3);
end
str = sprintf('%s iter %i Kurtosis %1.5f',pars.title, iter, abs(optimValues.fval));
title(str);
axis([min(X(:,1)), max(X(:,1)), min(X(:,2)), max(X(:,2))]);


% Plot KDE for projected density
y = linspace(min(X(:,1)), max(X(:,1)),200)';
f = fgt_kde(X(:,1), y, 0.9*size(X,1)^(-0.2)*std(X(:,1)));
ax1_pos = get(gca,'position');
	ax2 = axes('Position',ax1_pos,...
	'XAxisLocation','bottom',...
	'YAxisLocation','right',... 
	'XLim',[min(X(:,1)),max(X(:,1))],...
	'YLimMode','manual',...
	'YLim', [0, 2*max(f)],...
	'Color','none');
line(y, f, 'Parent',ax2,'Color','k','LineWidth',2);

hold off;
drawnow;
