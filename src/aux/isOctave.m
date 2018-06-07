function flag = isOctave()
% Determines whether the environment is Octave (returns TRUE) or MATLAB
%FLAG = ISOCTAVE()

flag=1;
if exist('OCTAVE_VERSION', 'builtin') == 0,
	flag=0;
end

