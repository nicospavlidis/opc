function Data = normalise(Data, standard);
%Removes observations (rows) with missing values and variables (columns) with no variation
%DATA = NORMALISE(DATA, STANDARD);
% 
% Inputs:
%	(DATA): N-by-D data matrix (observations stored row-wise)
%	(STANDARD): If (STANDARD)=1, then variables are standardised. Default = 0

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if nargin < 2,
	standard=0;
end

% Eliminate columns with no variation
inactive = (max(Data) - min(Data))==0;
if sum(inactive) > 0,
	Data = Data(:, ~inactive);
end

Data = bsxfun(@minus, Data, mean(Data));
if standard,
	Data = bsxfun(@rdivide, Data, std(Data)); 
end

%Data = Data - repmat(mean(Data), size(Data,1), 1);
%Data = Data./repmat(std(Data), size(Data,1), 1);

