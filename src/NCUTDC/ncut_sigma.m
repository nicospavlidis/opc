function s = ncut_sigma(X,smult)
%Default scaling parameter employed by Gaussian kernel in SCPP
%S = NCUT_SIGMA(X,SMULT)
%
% Default scaling parameter for minimum normalised cut projection pursuit
% Input:
% 	(X) Data matrix
% 	(smult) Multiplier (recommended value 100)

if nargin < 2,
	smult = 100;
end
[~,~,eigVal]  = pcacomp(X,1);
s = smult*sqrt(eigVal(1))*size(X,1)^(-0.2);
