function out = sim_transform(proj, a, d, weights)
%Transformation of points to compute pairwise similarities of projected data Eqs.(17)-(18)
%OUT = SIM_TRANSFORM(PROJ, A, D, WEIGHTS)
%
% Returns:
%	(OUT): Points after transformation
%
% Inputs:
%	(PROJ): Projected data
%	(A): Range parameter (alpha)
%	(D): Delta parameter of transformation function
%	(WEIGHTS): Number of original observations associated with each micro-cluster center
%		(If micro-clustering has been used)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

% no micro-clustering has been performed
if isempty(weights),
	weights = ones(size(proj,1),1);
end

% number of observations = sum of elements of micro-clusters
N = sum(weights);

pdim = size(proj,2);

out = proj;
c1 = (d*(1-d))^(1./d);
c2 = d * (d*(1-d))^((1.-d)/d);
for i=1:pdim,
	% standard deviation along each projection direction
	sigma = sqrt( sum(weights.*(proj(:,i).^2))/N );

	ixlo = find(proj(:,i) < -a*sigma);
	ixhi = find(proj(:,i) >  a*sigma);

	out(ixlo,i) = -a*sigma - d*(-a*sigma - proj(ixlo,i) + c1).^(1-d) + c2;
	out(ixhi,i) =  a*sigma + d*( proj(ixhi,i) - a*sigma + c1).^(1-d) - c2;
end
end
