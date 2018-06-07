try
	if exist ('OCTAVE_VERSION', 'builtin') == 0,
		echo on;
		mex -compatibleArrayDims kdeC.cpp
		mex -compatibleArrayDims dkdeDxC.cpp
		mex -compatibleArrayDims dval.c
		mex -compatibleArrayDims fgt_model.c
		mex -compatibleArrayDims fgt_predict.c
		echo off;
	else
		echo on;
		mex kdeC.cpp
		mex dkdeDxC.cpp
		mex dval.c
		mex fgt_model.c
		mex fgt_predict.c
		echo off;
	end
catch ME
	idSegLast = regexp(ME.identifier, '(?<=:)\w+$', 'match');
	if(~isempty(idSegLast))
		disp('Mex compilation error, please be sure to setup your compiler by typing ''mex -setup'' ')
	end
end
