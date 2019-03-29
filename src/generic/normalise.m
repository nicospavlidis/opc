function [Data,S] = normalise(Data, standard);
%Removes observations (rows) with missing values and variables (columns) with no variation
%DATA = NORMALISE(X, STANDARD);
% 
% Inputs:
%	(X): N-by-D data matrix (observations stored row-wise)
%	(STANDARD): If (STANDARD)=1, then variables/ columns are standardised. 
%		If (STANDARD)=='whiten' the whitening transformation is performed.
%		default = 0
%
% Output:
%	(X): N-by-D data matrix
%	(S):	In (STANDARD)=='whiten' (S) is the Cholesky decomposition of (X)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if nargin < 2 || isempty(standard),
	standard=0;
end

% Eliminate columns with no variation
inactive = (max(Data) - min(Data))==0;
if sum(inactive) > 0,
	Data = Data(:, ~inactive);
end

Data = bsxfun(@minus, Data, mean(Data));
if strcmp(standard,'whiten'),
	% Cholesky decomposition of data matrix
	S = chol(cov(Data));
	Data = Data*pinv(S);
elseif standard,
	Data = bsxfun(@rdivide, Data, std(Data)); 
end

%Data = Data - repmat(mean(Data), size(Data,1), 1);
%Data = Data./repmat(std(Data), size(Data,1), 1);

