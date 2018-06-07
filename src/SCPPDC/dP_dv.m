function DT = dP_dv(proj, X, weights, a, d)
%
% computes derivative of transformed projections along w.r.t. v, where v is a single projection vector
% 
if isempty(weights),
	weights = ones(size(proj,1),1);
end
N = sum(weights);

c1 = (d*(1-d))^(1./d);
DT = X;

sigma = sqrt( sum(weights.*(proj.^2))/N );
% debugged
dsigma_dv = (1./(sigma*N)) * sum(diag(weights.*proj)* X, 1);

% FINITE DIFFERENCE DEBUGGING OF ABOVE CODE
%pert=1.0e-4;
%finDiff = zeros(size(X,2), size(v,2));
%for j=1:size(v,2),
%	for i=1:size(v,1),
%		vn = v;
%		vn(i,j) = v(i,j) + pert;
%		s1 = sqrt( sum( weights.*(X*vn(:,j)).^2 )/N );
%
%		vn(i,j) = v(i,j) - pert;
%		s2 = sqrt( sum( weights.*(X*vn(:,j)).^2 )/N );
%
%		finDiff(i,j) = 0.5*(s1-s2)/pert;
%	end
%end
%maxD = max(max(abs(dsigma_dv-finDiff)));
%fprintf('Maximum discrepancy in dSigma_dV: %d\n', maxD);

ixlo = find(proj < -a*sigma);

%dim = size(X,2);
%for i=1:length(ixlo),
%	k = ixlo(i);
%	for j=1:dim,
%		DT(k,j) = -a*dsigma_dv(j) - d*(1-d)* (-a*dsigma_dv(j) - X(k,j))*(-a*sigma - proj(k) + c1)^(-d);
%	end
%end
%
%% debugged
DT(ixlo,:) = bsxfun(@times, d*(1-d)*(-a*sigma - proj(ixlo) + c1).^(-d), bsxfun(@plus, X(ixlo,:), a*dsigma_dv));
DT(ixlo,:) = bsxfun(@plus, DT(ixlo,:), -a*dsigma_dv);


ixhi = find(proj >  a*sigma);
%for i=1:length(ixhi),
%	k = ixhi(i);
%	for j=1:dim,
%		DT(k,j) = a*dsigma_dv(j) + d*(1-d)* (-a*dsigma_dv(j) + X(k,j))*(-a*sigma + proj(k) + c1)^(-d);
%	end
%end
% debugged
DT(ixhi,:) = bsxfun(@times, (1-d)*d*(-a*sigma + proj(ixhi) + c1).^(-d), bsxfun(@plus, X(ixhi,:), -a*dsigma_dv));
DT(ixhi,:) = bsxfun(@plus, DT(ixhi,:), a*dsigma_dv);

end

%% Finite Difference debugging (this requires the data and the projection matrix)
%Max = 0;
%MaxJ = 0;
%for j=1:size(v,2),
%	DT = dP_dv(proj(:,j), X, weights, a, d);
%	pert=1.0e-5;
%	finDiff = zeros(size(X));
%	for i=1:length(v),
%		vn = v;
%		vn(i,j) = v(i,j) + pert;
%
%		pr = X*vn(:,j);
%		p1 = sim_transform(pr, a, d, weights);
%
%		vn(i,j) = v(i,j) - pert;
%		pr = X*vn(:,j);
%		p2 = sim_transform(pr, a, d, weights);
%
%		finDiff(:,i) = 0.5*(p1-p2)./pert;
%	end
%	fprintf('column %i, %d\n',j, max(max(abs(finDiff - DT))));
%	if max(max(abs(finDiff - DT))) > Max,
%		Max = max(max(abs(finDiff - DT)));
%		MaxJ = j;
%	end
%end
