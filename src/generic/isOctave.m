function flag = isOctave()
%Determines whether the environment is GNU Octave (returns TRUE) or MATLAB
%FLAG = ISOCTAVE()

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

flag=1;
if exist('OCTAVE_VERSION', 'builtin') == 0,
	flag=0;
end

