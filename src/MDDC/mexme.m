% Script to compile C++ functions for 1D kernel density estimation
% with Gaussian kernels

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------

if exist ('OCTAVE_VERSION', 'builtin') == 0,
	mex -compatibleArrayDims kdeC.cpp
	mex -compatibleArrayDims dkdeDxC.cpp
else
	mex kdeC.cpp
	mex dkdeDxC.cpp
end
