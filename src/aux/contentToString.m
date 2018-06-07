function str = contentToString(sp)
	str = sprintf('(%i) n:%i', sp.tree_params.nodeid, length(sp.idx));
	if isfield(sp.tree_params,'purity') & isfield(sp.tree_params,'leaf') & ... 
		~isempty(sp.tree_params.purity) & sp.tree_params.leaf,

		str = strcat(str, sprintf(' P:%1.2f', sp.tree_params.purity));

	elseif isfield(sp.tree_params,'SR') & ~isempty(sp.tree_params.SR),

		str = strcat(str, sprintf(' SR:%1.2f', sp.tree_params.SR));
	end
end

