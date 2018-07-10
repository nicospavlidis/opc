function out = ifelse(a,b,c)
%Shorthand for ternary operator: if-then-else
%OUT = IFELSE(A,B,C)
%
% Ternary operator for MATLAB like the one in GNU Octave

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if a,
	out = b;
else
	out = c;
end
