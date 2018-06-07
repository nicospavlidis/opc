function [f,df] = f_df_md(v,X, pars)
%function [f,df] = f_df_md(v,X, pars)
%
% Returns: 
%	(f) the value of projection index for MDH 
%	(df) derivative of projection index w.r.t. projection vector (v)
% Inputs:
% 	(v) Projection vector
%	(X) Data matrix
%	(pars) Structure array containing all parameters of MDH algorithm
%
% Used by MATLAB optimisation algorithm

[f, bmin] = f_md(v,X,pars);
if nargout > 1,
	df = df_md(v,X,pars,bmin);
end
