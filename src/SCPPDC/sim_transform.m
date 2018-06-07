% DEBUGGED!
function out = sim_transform(proj, a, d, weights)
%
% Computes transformation of pairwise similarities of projected data Eqs.(17)-(18)
%


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
