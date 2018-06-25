function [v, fval, idx] = pddppp(Data, pars)
%Projection Pursuit function for PDDP algorithm
%[V, FVAL, IDX] = PDDPPP(X, PARS)
%
% Returns:
%	(V): Projection vector: 1st Principal Component of (X)
%	(FVAL): Total scatter of (X)
%	(IDX): Binary cluster assignment
%
% Inputs:
%	(X): N-by-D data matrix
%	(PARS): Parameter structure containing algorithm parameters

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

[N,dim] = size(Data);

invalid = false;
if dim == 1,
	fprintf('Effective dimensionality of 1: No point in projection pursuit\n');
	invalid = true;
end

if isfield(pars,'minsize') & 2*pars.minsize > N,
	%fprintf('Too few observations to split\n');
	invalid = true;
end

v = pcacomp(Data,1);
if invalid,
	fval = -inf;
	idx = ones(N,1);
	return;
end

proj = Data*v;
idx = (proj > mean(proj)) + 1;
fval = total_scatter(Data);
