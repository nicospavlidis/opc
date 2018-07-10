#include <iostream>
#include <math.h>
#include <stdio.h>
#include "mex.h"

//-------------------------------------------------------------------------------------
// Copyright @ Nicos Pavlidis, 2018
// OPC is licensed under the BSD-3-Clause License - see the LICENSE.md file for details
//-------------------------------------------------------------------------------------

using namespace std;

// debugged
double kde(const double &x, const int &N, const double *proj, const double &h)
{
	const double PI = 3.141592654;
	double h2= 2*h*h;
	double y = 0.;
	for (int i=0; i < N; i++) {
		y += exp((proj[i]-x)*(x - proj[i])/h2);
	}
	return y/(N*sqrt(2.*PI)*h);
}


/* the gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	/* ensure input is a column vector */
	if (mxGetN(prhs[0]) != 1 && mxGetM(prhs[0]) !=1) { mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notColumnVector","Input matrix must be Column Vector."); }

//	/* Examine input (right-hand-side) arguments. */
//	mexPrintf("\nThere are %d right-hand-side argument(s).", nrhs);
//	for (int i=0; i<nrhs; i++) { mexPrintf("\n\tInput Arg %i is of type:\t%i ",i,mxGetM(prhs[i])); }
//
//	/* Examine output (left-hand-side) arguments. */
//	mexPrintf("\n\nThere are %d left-hand-side argument(s).\n", nlhs);
//	if (nlhs > nrhs) { mexErrMsgIdAndTxt( "MATLAB:mexfunction:inputOutputMismatch", "Cannot specify more outputs than inputs.\n"); }
//
//	for (int i=0; i<nlhs; i++)  {
//		plhs[i]=mxCreateDoubleMatrix(1,1,mxREAL);
//		*mxGetPr(plhs[i])=(double)mxGetNumberOfElements(prhs[i]);
//	}

	// for (int i=1; i <5; i++) { if (mxGetN(prhs[i]) != 1 || mxGetM(prhs[i]) != 1) { cout << "Input " << i << " must be a real number:\nInputs: (proj, h, alpha,beta,gamma)"<< endl; } }

	/* Get number of rows of input vector (mxGetN is for columns) */ 
//	size_t N = mxGetM(prhs[0]);
//	/* create a pointer to the real data in the input matrix  */
//	double *proj = mxGetPr(prhs[0]);
//	double h = *mxGetPr(prhs[1]);
//	double alpha = *mxGetPr(prhs[2]);
//	double beta = *mxGetPr(prhs[3]);
//	double gamma = *mxGetPr(prhs[4]);
//	double bandMultiplier = *mxGetPr(prhs[5]);
//	double kern = *mxGetPr(prhs[6]);

	double x = *mxGetPr(prhs[0]);
	int N = mxGetM(prhs[1]);
	double *proj  = mxGetPr(prhs[1]);
	double h = *mxGetPr(prhs[2]);

	double k = kde(x,N,proj,h);
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0]) = k;
}
