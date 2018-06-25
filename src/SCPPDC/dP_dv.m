function DT = dP_dv(proj, X, weights, a, d)
%Derivative of transformed data projections along projection vector V
%DT = DP_DV(PROJ, X, WEIGHTS, A, D)
%
% Returns:
%	(OUT): Derivative of transformed data projections along (V)
%
% Inputs:
%	(PROJ): Projected data 
%	(X): Data matrix
%	(WEIGHTS): Weight assigned to each observation (in case of micro-clustering)
%	(A): alpha parameter used to determine range
%	(D): parameter used by distance transformation function

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if isempty(weights),
	weights = ones(size(proj,1),1);
end
N = sum(weights);

c1 = (d*(1-d))^(1./d);
DT = X;

sigma = sqrt( sum(weights.*(proj.^2))/N );
dsigma_dv = (1./(sigma*N)) * sum(diag(weights.*proj)* X, 1);

ixlo = find(proj < -a*sigma);
DT(ixlo,:) = bsxfun(@times, d*(1-d)*(-a*sigma - proj(ixlo) + c1).^(-d), bsxfun(@plus, X(ixlo,:), a*dsigma_dv));
DT(ixlo,:) = bsxfun(@plus, DT(ixlo,:), -a*dsigma_dv);

ixhi = find(proj >  a*sigma);
DT(ixhi,:) = bsxfun(@times, (1-d)*d*(-a*sigma + proj(ixhi) + c1).^(-d), bsxfun(@plus, X(ixhi,:), -a*dsigma_dv));
DT(ixhi,:) = bsxfun(@plus, DT(ixhi,:), a*dsigma_dv);

end
