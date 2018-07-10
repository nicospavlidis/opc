function str = contentToString(sp)
%Function used to print ctree Node to terminal
%STR = CONTENTTOSTRING(SP)
%
% Inputs:
%	(SP): Node from ctree object
%
% Output:
%	(STR): String description of Node
%
%  Based on MATLAB tree class function orinally created by Jean-Yves Tinevez <tinevez@pasteur.fr> March 2012
%  Modified by: Nicos Pavlidis 2018

	str = sprintf('(%i) n:%i', sp.tree_params.nodeid, length(sp.idx));
	if isfield(sp.tree_params,'purity') & isfield(sp.tree_params,'leaf') & ... 
		~isempty(sp.tree_params.purity) & sp.tree_params.leaf,

		str = strcat(str, sprintf(' P:%1.2f', sp.tree_params.purity));

	elseif isfield(sp.tree_params,'SR') & ~isempty(sp.tree_params.SR),

		str = strcat(str, sprintf(' SR:%1.2f', sp.tree_params.SR));
	end
end


