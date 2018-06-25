function tsc = total_scatter(X)
%Total scatter used as split index criterion by PDDP
%TSC = TOTAL_SCATTER(X)
%
% Inputs:
%	(X): Data matrix
%
%Reference:
%D. Boley. Principal Direction Divisive Partitioning. Data Mining and Knowledge Discovery, 2(4):325-344, 1998.

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

tsc = norm(bsxfun(@minus, X, mean(X)),'fro')^2;
