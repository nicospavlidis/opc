function T = mycrosstab(x,y)
%Cross tabulation table for vectors X,Y
%T = MYCROSSTAB(X,Y)
%
% Included for compatibility with Octave
% In MATLAB: mycrosstab(x,y) = crosstab(x,y)
%
% Vectors (X) and (Y) must be of the same length
if ~isOctave(),
	T = crosstab(x,y);
	return;
end

assert(length(x)==length(y),'Incompatible vector lengths');
xl = unique(x);
yl = unique(y);

T = zeros(length(xl), length(yl));
for i=1:length(xl),
	for j=1:length(yl),
		T(i,j) = length( find(x==xl(i) & y==yl(j)) );
	end
end
