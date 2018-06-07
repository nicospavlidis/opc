function [fmin, bmin, df, dkde] = f_md(v,X,pars)
%function [fval, split, df, dkde] = f_md(v,X,pars)
%
% Returns: 
%	(fval) penalised density integral for optimal hyperplane orthogonal to (v)
%	(split) optimal split point on (v) after (v) is scaled to unit-length
%	(df) derivative of penalised density integral at (split)
%	(dkde) derivative of 1D kernel density estimator at (split)	
% Inputs:
% 	(v) Projection vector
%	(X) Data matrix
%	(pars) Structure array containing MDH parameters (alpha, eta, epsilon,h)
%
% Used by MATLAB optimisation algorithm

alpha = pars.alpha;
eta = pars.eta;
epsilon = pars.epsilon;

h = pars.bandwidth;
if isa(h,'function_handle'),
	h = pars.bandwidth(X,pars);
	assert(numel(h)==1 & h>0,'f_md: Error incorrect bandwidth');
end

proj = sort( X*(v./norm(v,2)) );

if alpha < 1.0e-3,
	bmin = 0;
	fmin = kdeC(bmin,proj,h);
	if nargout > 2, 
		[df, dkde] = dpkde(bmin, proj, h, alpha, eta, epsilon);
	end
	return
end

%%%%%%%%%%% Create necessary structures for fast computation of derivatives
N = length(proj);
s = std(proj);

h2 = sqrt(2)*h;
[xc1, A_k1] = fgt_model(proj', proj'    , h2, 30, min(N,2*sqrt(N)), 8);
[xc2, A_k2] = fgt_model(proj', ones(1,N), h2, 30, min(N,2*sqrt(N)), 8);

denom = 1/(N*(h^3)*sqrt(2*pi));
L = exp(-0.5)/(sqrt(2*pi)*(h*h)* eta^epsilon);

df_handle = @(y)df_fgt(y, xc1, A_k1, xc2, A_k2, h2, 1/(N*(h^3)*sqrt(2*pi)), s,alpha,epsilon,L);
f_handle = @(y)f_fgt(y, xc2, A_k2, h2, 1/(N*h*sqrt(2*pi)), s, alpha, epsilon, L);

%%%%%%%%%%%

%keyboard;
x = linspace(-alpha*s, alpha*s,200)';
x(end) = alpha*s;

bd = [-alpha*s, alpha*s];
% Newton-Raphson to identify minimiser at boundary of feasible region
for i=1:2,
	j=1;
	xcur = bd(i);
	%grad = dpkde(xcur, proj, h, alpha, eta, epsilon);
	grad = df_handle(xcur);

	%if verbose, fprintf('Bd: %1.5f\tInitial gradient %1.5f\n', bd(i), grad); end
	if (-1)^(i-1)*grad > 0,
		while (abs(grad) > 1.0e-8) & (j < 11),
			d1f = dkdeDxC(xcur, proj, h);
			d2f = d2kdedx(xcur, proj, h);

			if i==2,
				dp = (1+epsilon) * (xcur - alpha*s)^epsilon * L;
			else
				dp = -1*(1+epsilon)*(-xcur - alpha*s)^epsilon * L;
			end
			d2p = (1+epsilon) * L;
			Dx = -(d1f+dp)/(d2f + d2p);
			xcur = xcur + Dx;
			%grad = dpkde(xcur, proj, h, alpha, eta, epsilon);
			grad = df_handle(xcur);

			% in practice not more than 3 iterations are enough
			j = j+1;
		end
	end
	bd(i) = xcur;
end

f1 = f_handle(bd(1)); %pkde(bd(1),proj,h,alpha,eta,epsilon);
f2 = f_handle(bd(2)); %pkde(bd(2),proj,h,alpha,eta,epsilon);
if abs(f1 - f2) < eps,
	fmin = f1;
	bmin = bd;
else
	[fmin,id] = min([f1,f2]);
	bmin = bd(id);
end

% estimate the gradient at grid points
%df = dpkde(x, proj, h, alpha, eta, epsilon);
df = df_handle(x);

for i=2:length(x)-2,
	%fprintf('%d\t%d\n',f(i),df(i));
	fx = -1;
	dfx = -1;
	% gradient is zero
	if abs(df(i)) < 1.0e-8 || abs(df(i+1)) < 1.0e-8,
		% set minimiser
		if abs(df(i)) < abs(df(i+1)),
			xmin = x(i);
			dfx = df(i);
		else
			xmin = x(i+1);
			dfx = df(i+1);
		end
		% estimate minimum 
		fx = f_handle(xmin); %pkde(xmin, proj, h, alpha, eta, epsilon);

	% change of slope from negative to positive-> candidate for minimiser
	elseif df(i)<0.0 && df(i+1)>0.0,

		%fhandle = @(y)dpkde(y, proj, h, alpha, eta, epsilon);

		% fzero function is also available in Octave
		%options = optimset('FunValCheck','on');
		%options = optimset('TolX',eps);
		[xmin, dfx, flag] = fzero(df_handle, [x(i), x(i+1)]);
		if flag~=1,
			fprintf('Convergence problems in fzero. flag = %i\n',flag);
		end
		fx = f_handle(xmin); %pkde(xmin, proj, h, alpha, eta, epsilon);
	end

	if fx ~= -1,
		%if verbose==1, fprintf('%1.5f\tdfdx(x): %1.5f\tf(x): %1.5f\n', xmin, dfx, fx); end
		if abs(dfx)>1.0e-6,
			fprintf('Convergence problems in global_min: Grad at minimum= %f\n', dfx);
		end

		% Unique minimiser/ Differentiability check
		if abs(fx - fmin) < 1.0e-12,
			ident = false;
			for j = 1:length(bmin),
				if abs(xmin-bmin(j)) < 1.0e-6,
					ident = true;
					break
				end
			end
			if ~ident, 
				bmin = [bmin; xmin];
			end

		elseif fx < fmin,
			% bmin = zeros(1,1);
			bmin = [xmin];
			fmin = fx;
		end
	end
end

if fmin < 1.0e-4 & length(bmin)>1,
	[~,max_gap] = max(diff(proj));
	mid = 0.5*(proj(max_gap+1) + proj(max_gap));
	fmid = f_handle(mid); %pkde(mid, proj, h, alpha, eta, epsilon);
	if fmid < fmin,
		bmin = [mid];
		fmin = fmid;
	end
end

if length(bmin)>1,
	warning('More than one global minimisers: ');
	%keyboard;
	for i=1:length(bmin),
		fprintf('%1.8f  ', bmin(i));
	end
	fprintf('\n');
end

if nargout>2,
	%[df, dkde] = dpkde(bmin(1), proj, h, alpha, eta, epsilon);
	[df, dkde] = df_handle(bmin(1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = f_fgt(x, xc,A_k, h2,denom, s, alpha, epsilon, L)
fkde = denom*fgt_predict(x', xc, A_k, h2, 30)';
f = fkde + L* max([zeros(length(x),1), x-alpha*s, -x-alpha*s], [],2).^(1+epsilon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [df, dkde] = df_fgt(x, xc1,A_k1,xc2,A_k2, h2,denom,s,alpha,epsilon,L)
dkde = (fgt_predict(x', xc1,A_k1, h2, 30)' - x.*fgt_predict(x', xc2,A_k2,h2, 30)')*denom;
% derivative of penalty term
if length(x) == 1,
	if x < -alpha*s,
		der = -(1+epsilon) * (-alpha*s-x)^epsilon * L;
	elseif x > alpha*s,
		der =  (1+epsilon) * (x-alpha*s)^epsilon * L;
	else
		der = 0;
	end
else
	[~, index] = max([zeros(length(x),1), x-alpha*s, -x-alpha*s], [],2);
	der = zeros(length(x),1);
	der(index==2) =  L *(1+epsilon)*( x(index==2)-alpha*s).^epsilon;
	der(index==3) = -L *(1+epsilon)*(-x(index==3)-alpha*s).^epsilon;
end
df = dkde + der;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = d2kdedx(x, proj, h)
out = sum( (proj-x-h).*(proj-x+h) .* normpdf(proj, x, h) )/(h^4 * length(proj));
