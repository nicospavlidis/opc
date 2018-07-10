function step = drsc_linesearch(Data, W, max_step, f, grad, sigma, U, degs)
%Line search subroutine to determine stepsize for gradient ascent performed by DRSC
%STEP = DRSC_LINESEARCH(DATA, W, MAX_STEP, F, GRAD, SIGMA, U, DEGS)
%
%Line-search for gradient ascent within DRSC satisfying only first Wolfe condition 
%(using both takes excessively long)
%
% Inputs:
%	(DATA): N-by-D data matris
%	(W): 	Projection matrix (ascent is performed w.r.t. last column of W)
%	(MAX_STEP): Maximum stepszie	
%	(F):	Function value at current poit
%	(GRAD):	Gradient at current point
%	(SIGMA): Scale parameter for Gaussian kernel (used to estimate similarity matrix)
%	(U):	Top K eigenvectors of D^{-1/2} *K*D^{-1/2}
%	(DEGS): Vector of degrees of each vertex: D = diag(degs)
%
% Output:
%	(STEP): Stepsize satisfying 1st Wolfe condition
%		(Returns 0 if no such stepsize is found)

%-------------------------------------------------------------------------------------
% Copyright @ Nicos Pavlidis, 2018
% OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
%-------------------------------------------------------------------------------------


% Eq. (6) in Niu,Dy and Jordan AISTATS (2011)
projGrad = gram_schmidt(grad, W);

step = alphamax;
fail = 1;
for i=1:20,
	% Eq. (7) in Niu,Dy and Jordan AISTATS (2011)
	Wtry = [W(:,1:end-1), sqrt(1-step^2)*W(:,end) + step*projGrad];
	f1 = f_df_drsc(Data, Wtry, sigma, U, degs);

	if f1 > f + 1.0e-4*step*(grad'*projGrad),
		fail=0;
		break;
	end
	step = 0.5*step;
end
if fail==1,
	step=0;
end

