function Data = normalise(Data, standard);
%Removes observations (rows) with missing values and variables (columns) with no variation
%DATA = NORMALISE(DATA, STANDARD);
% 
% Inputs:
%	Data: N-by-D data matrix (observations stored row-wise)
%	standard: Optional argument with default value 0. If (standard)=1, then variables are standardised

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

