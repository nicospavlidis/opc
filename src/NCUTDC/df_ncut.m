function out = df_ncut(v,X,pars)
%Derivative of minimum normalised cut projection index with respect to projection vector
%DF = DF_NCUT(V,X,SIGMA,MINN)
%
% Inputs:
%	v: projection vector
%	X: N-by-D Data matrix
%	pars: parameter struct containing 
%		(minsize): minimum cluster size and
%		(sigma) Scaling parameter for Laplace kernel
%
% Output: 
%	df: derivative of Normalised Cut projection index

if isa(pars.sigma,'function_handle')
	sigma = pars.sigma(X);
else
	sigma = pars.sigma;
end

minN = pars.minsize;

[x, index] = sort( X*( v./(norm(v,2)*sigma) ) );
n = length(x);

CP = exp(x(1) - x);
CSPi = cumsum(1./CP);
c = cumsum( CP(n:-1:1) );
S = c((n-1):-1:1) .* CSPi(1:n-1);
updown = cumsum(CP.*CSPi);
downup = cumsum( [c((n-1):-1:1); 0]./CP );
DL = updown+downup;

r = S.*( 1./DL(1:(n-1)) +  1./(DL(n) - DL(1:(n-1))) );
[~,w] = min( r(minN:(n-minN)) );
w = w + minN - 1;

Ek = exp(x - x(w));
Eki = 1./Ek;
Ck = cumsum(Ek);
Cki = cumsum(Eki);

ds = zeros(1,n);
if w > 2,
	ds(1) = (Cki(n) - Cki(w))*Ek(1)*(1./DL(w) + 1./(DL(n)-DL(w))) ...
		- S(w)*( (2*Ek(1)*(Cki(w) - Cki(1)) + Ek(1)*(Cki(n) - Cki(w)))/(DL(w)^2) ...
		+ Ek(1)*(Cki(n)-Cki(w))/(DL(n)-DL(w))^2 );

	ds(2:w-1) = ( (Cki(n) - Cki(w))*(1./DL(w) + 1./(DL(n)-DL(w))) )*Ek(2:w-1) - S(w)*( ...
		(2*Ek(2:w-1).*(Cki(w)-Cki(2:w-1)) - 2*Eki(2:w-1).*Ck(1:w-2) + Ek(2:w-1)*(Cki(n)-Cki(w)))/(DL(w)^2) ...
		+ Ek(2:w-1)*( (Cki(n)-Cki(w))/((DL(n)-DL(w))^2) ) );
	
	ds(w) = (Cki(n)-Cki(w))*(1./DL(w) + 1./(DL(n)-DL(w))) - S(w)*( (Cki(n)-Cki(w)-2*Ck(w-1))/DL(w)^2 ...
		+(Cki(n) - Cki(w))/(DL(n)-DL(w))^2 );

	ds(w+1:n) = -Ck(w)*(1./DL(w) + 1./(DL(n)-DL(w)))*Eki(w+1:n) - S(w)*( -(Ck(w)/(DL(w)^2))*Eki(w+1:n) ...
		+ (2*Ek(w+1:n).*(Cki(n)-Cki(w+1:n)) - 2*Eki(w+1:n).*(Ck(w:n-1)-Ck(w)) ...
		- Ck(w)*Eki(w+1:n))/(DL(n)-DL(w))^2 );

else
	
	d1 = (Cki(n)-Cki(w))*Ek(1);
	d2 = 2*Ek(1)*(Cki(w)-Cki(1)) + d1;
	ds(1) = d1./DL(w) + d1./(DL(n)-DL(w)) - S(w)*(d2/(DL(w)^2) + d1/((DL(n)-DL(w))^2));

	if w == 2,
		d1 = Cki(n) - Cki(w);
		d2 = Cki(n) - Cki(w) - 2*Ck(w-1);
		ds(2) = d1*(1./DL(w) + 1./(DL(n)-DL(w))) - S(w)*(d2/(DL(w)^2) + d1/((DL(n)-DL(w))^2));
	end

	d1 = -Ck(w)*Eki(w+1:n);
	d3 = 2*Ek(w+1:n).*(Cki(n) - Cki(w+1:n)) - 2*Eki(w+1:n).*(Ck(w:n-1)-Ck(w)) - Ck(w)*Eki(w+1:n);
	ds(w+1:n) = d1.*(1./DL(w) + 1./(DL(n)-DL(w))) - S(w)*( d1./(DL(w)^2) + d3./((DL(n)-DL(w))^2) );
end
dv = (X(index,:) - (X(index,:)*v) * (v'./(v'*v)) )./(sigma*norm(v,2));
out = ds*dv;
