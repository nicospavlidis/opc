function [coeff,score,latent] = pcacomp(X,index)
%Returns the principal components of (X) specified in vector (index)
%[COEFF,SCORE,LATENT] = PCACOMP(X,INDEX)
%
% Returns:
%	(COEFF):  [pca(X,'NumComponents',index(1)), pca(X,'NumComponents',index(2)), ... ]
%	(SCORE):  Representation of X in the principal component space
%	(LATENT): Eigenvalues of COV(X'*X)
%
% Inputs:
%	(X): Data matrix
%	(INDEX): Vector containing indices of principal component vectors required

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if isempty(index),
	error('pcacomp: Unspecified number of PCA components');
elseif min(index<=0),
	error('pcacomp: index has to be positive');
end

% Determine whether we are using Octave or MATLAB
if exist ('OCTAVE_VERSION', 'builtin') == 0,
	if nargout ==1,
		coeff = pca(X,'NumComponents',max(index));
	else
		[coeff,score,latent] = pca(X,'NumComponents',max(index));
	end

	if ~isempty(coeff),
		coeff = coeff(:,index);
		if nargout > 1,
			score = score(:,index);
			latent = latent(index);
		end
	end
else
	% Octave
	if nargout==1,
		coeff = princomp(X);
	else
		[coeff,score,latent] = princomp(X);
	end

	coeff = coeff(:,index);
	if nargout>1,
		score = score(:,index);
		latent = latent(index);
	end
end
