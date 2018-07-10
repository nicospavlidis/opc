function [f,df] = f_df_md(v,X, pars)
%Function value and derivative for penalised density 
%[F,DF] = F_DF_MD(V,X, PARS)
%
% Inputs:
%	(V): Projection vector
%	(X): Data matrix
%	(PARS): Structure array containing all parameters of MDH algorithm
%
% Output: 
%	(F): the value of projection index for MDH 
%	(DF): derivative of projection index w.r.t. projection vector (v)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

[f, bmin] = f_md(v,X,pars);
if nargout > 1,
	df = df_md(v,X,pars,bmin);
end
