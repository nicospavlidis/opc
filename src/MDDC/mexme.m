if exist ('OCTAVE_VERSION', 'builtin') == 0,
	mex -compatibleArrayDims kdeC.cpp
	mex -compatibleArrayDims dkdeDxC.cpp
else
	mex kdeC.cpp
	mex dkdeDxC.cpp
end
