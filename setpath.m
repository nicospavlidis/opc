addpath(genpath(pwd))

if isOctave,
	warning('off', 'Octave:possible-matlab-short-circuit-operator');
	pkg load statistics;
	pkg load optim;
end

