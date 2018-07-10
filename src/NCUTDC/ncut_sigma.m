function s = ncut_sigma(X,smult)
%Default scaling parameter employed by Gaussian kernel in minimum normalised cut projection pursuit
%S = NCUT_SIGMA(X,SMULT)
%
% Input:
% 	(X): Data matrix
% 	(SMULT): Multiplier (recommended value 100)
%
% Output:
%	(S): Default scaling parameter for NCUTPP

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if nargin < 2,
	smult = 100;
end
[~,~,eigVal]  = pcacomp(X,1);
s = smult*sqrt(eigVal(1))*size(X,1)^(-0.2);
