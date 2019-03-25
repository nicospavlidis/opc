function [coeff,score,latent] = pcacomp(X,index)
%Returns the principal components of (X) specified in vector (index)
%[COEFF,SCORE,LATENT] = PCACOMP(X,INDEX)
%
% Inputs:
%	(X): Data matrix
%	(INDEX): Vector containing indices of principal component vectors required
%
% Output:
%	(COEFF):  [pca(X,'NumComponents',index(1)), pca(X,'NumComponents',index(2)), ... ]
%	(SCORE):  Representation of X in the principal component space
%	(LATENT): Eigenvalues of COV(X'*X)

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
if ~isOctave(),
	if nargout == 1,
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
	% Octave (code below is much faster than 'princomp')
	Xc = bsxfun(@minus,X,mean(X));
	[~,S,coeff] = svds(Xc,max(index));
	if ~isempty(coeff),
		coeff = coeff(:,index);
		score = Xc*coeff;
		latent = diag(S).^2/(size(Xc,1)-1);
		latent = latent(index);
	end

	%if nargout==1,
	%	%coeff = princomp(X);
	%else
	%	%[coeff,score,latent] = princomp(X);
	%end

	%if nargout>1,
	%	score = score(:,index);
	%	latent = latent(index);
	%end
end
