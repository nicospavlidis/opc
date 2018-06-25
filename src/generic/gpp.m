function [optS, idx, spindex] = gpp(X, ppf, pars)
%Generic Projection Pursuit 
%[OPTS, IDX, SPINDEX] = GPP(X, PPF, CLF, PARS)
%
% Returns:
%	(OPTS): Binary separator object of class GSEP
%	(IDX): Binary assignment vector in {-1,1}
%	(SPINDEX): Split index. By default equal to projection index (FVAL)
%		If PARS.SPLIT_INDEX is empty SPINDEX=FVAL
%		otherwise (SPINDEX) is calculated from the user-defined 
%		function (PARS.SPLIT_INDEX(V,X,PARS))
%
% Inputs:
%	(X): Data matrix
%	(PPF): Function handle to projection pursuit function of the form:
%		[V,FVAL] = F(X,PARS), where (V) is the optimal projection matrix, and 
%		(FVAL) is the projection index
%	(CLF): Binary Cluster assignment function of the form:
%		CLUSTERS = F(V,X,PARS), where (clusters) is in {1,2}

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

[N,dim] = size(X);
invalid = false;
if dim == 1,
	fprintf('Effective dimensionality of 1: No point in projection pursuit\n');
	invalid = true;
end

if isfield(pars,'minsize') & 2*pars.minsize > N,
	invalid = true;
end

if invalid,
	optS = gsep(ones(dim,1)/sqrt(dim), -inf, -[1:size(X,1)]',ppf,pars);
	idx = -ones(N,1);
	spindex = -inf;
	return;
end

% Apply projection pursuit algorithm to obtain optimal 
[v,fval,idx] = ppf(X, pars);
% transform cluster labels from {1,2} to {-1,1}
idx = (2*idx-3);

spindex = fval;
if ~isempty(pars.split_index),
	if isa(pars.split_index,'function_handle'),
		spindex = pars.split_index(v,X,pars);
	elseif strcmp(pars.split_index, 'size'),
		spindex = size(X,1);
	end
end
optS = gsep(v, fval, idx.*[1:size(X,1)]', ppf, pars);
