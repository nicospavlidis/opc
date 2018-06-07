/*

 Evaluate weighted gaussian RBF between vectors x,y
 

 Usage
 -------

 v            = dval(x , y , [w] , [sigma] );

 Inputs
 -------

 x             Source data (d x Nx)
 y             Target data (d x Ny)
 w             Weigths (1 x Nx) ( default w = ones(1 , Nx) ) 
 sigma         Noise Standard deviation of the kernel (default sigma = 1)


 Ouputs
 -------

 v             density (1 x Ny)

	
 To compile
 ----------

 mex -output dval.dll dval.c

 mex  -f mexopts_intel10.bat -output dval.dll dval.c

 If compiled with the "OMP" compilation flag

 mex -DOMP -f mexopts_intel10.bat -output dval.dll dval.c "C:\Program Files\Intel\Compiler\11.1\065\mkl\ia32\lib\mkl_core.lib" "C:\Program Files\Intel\Compiler\11.1\065\mkl\ia32\lib\mkl_intel_c.lib" "C:\Program Files\Intel\Compiler\11.1\065\mkl\ia32\lib\mkl_intel_thread.lib" "C:\Program Files\Intel\Compiler\11.1\065\lib\ia32\libiomp5md.lib"

Example 1 
---------


d       = 10;
Nx      = 100;
Ny      = 10000;
x       = randn(d , Nx);
y       = randn(d , Ny);
v       = dval(x , y);


Example 2 
---------


d       = 10;
Nx      = 100;
Ny      = 10000;
x       = randn(d , Nx);
y       = randn(d , Ny);
u       = rand(1 , Nx);
sigma   = 2;
v       = dval(x , y , u , sigma);

 
Example 3 
---------


d       = 2;
R       = [2 , 0.4 ; 0.4  3];
Nx      = 100;
sigma   = 3;
vect    = (-5:0.3:5);
Ny      = length(vect);
w       = (1/Nx)*ones(1 , Nx);
  
x       = (chol(R)'*randn(d , Nx));

[X , Y] = meshgrid(vect);
y       = [X(:) , Y(:)]';

v       = dval(x , y , w , sigma);

densite = reshape( v , Ny , Ny);

figure
set(gcf , 'renderer' , 'opengl');
surfc(X , Y , densite)
shading interp
lighting phong

light
alpha(0.5);
hold on
plot(x(1 , :) , x(2 , :) , 'r+' , 'markersize' , 10);
hold off
view(2)
colorbar


  
  Author : Sébastien PARIS : sebastien.paris@lsis.org
  -------

*/


#include <math.h>
#include "mex.h"
#ifdef OMP
 #include <omp.h>
#endif

#ifndef MAX_THREADS
#define MAX_THREADS 64
#endif
#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif


/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/

void dval(double * , double * , double * , double  , double * , int  , int  , int );

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/

void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{	
	double *x , *y , *w;
	double sigma = 1.0;
	double *v;
	const int  *dimsx , *dimsy;
	int numdimsx , numdimsy ;
	int i , d , Nx , Ny;

	
	/*--------------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------------*/
	/* -------------------------- Parse INPUT  -------------------------------------- */
	/*--------------------------------------------------------------------------------*/	
	/*--------------------------------------------------------------------------------*/

	if ((nrhs < 2))
	{
			mexErrMsgTxt("Usage : v = dval(x , y , [w] , [sigma]);"); 
	}
	
	/* ----- Input 1 ----- */
	
	x           = mxGetPr(prhs[0]);
	numdimsx    = mxGetNumberOfDimensions(prhs[0]);
	dimsx       = mxGetDimensions(prhs[0]);
	
	d           = dimsx[0];
	Nx          = dimsx[1];
	
	/* ----- Input 2 ----- */
	
	y           = mxGetPr(prhs[1]);
    numdimsy    = mxGetNumberOfDimensions(prhs[1]);
	dimsy       = mxGetDimensions(prhs[1]);
	Ny          = dimsy[1];
	
	/* ----- Input 3 ----- */
	
	if ((nrhs < 3) || mxIsEmpty(prhs[2]))
	{	
		w    = (double *)mxMalloc(Nx*sizeof(double));	
		for (i = 0 ; i < Nx ; i++)
		{
			w[i] = 1.0;	
		}
	}
	else
	{	
		w           = mxGetPr(prhs[2]);	
	}
	
	/* ----- Input 4 ----- */
	
	if (nrhs == 4)
	{			
		sigma       = (double)mxGetScalar(prhs[3]);	
	}
	
	/*--------------------------------------------------------------------------------*/
	/*---------------------------------------,----------------------------------------*/
	/* -------------------------- Parse OUTPUT  ------------------------------------- */
	/*--------------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------------*/
	
	/* ----- output 1 ----- */

	plhs[0]        = mxCreateDoubleMatrix(1 , Ny , mxREAL);	
	v              = mxGetPr(plhs[0]);
		
	/*---------------------------------------------------------------------------------*/
	/*---------------------------------------------------------------------------------*/
	/* ----------------------- MAIN CALL  -------------------------------------------- */
	/*---------------------------------------------------------------------------------*/
	/*---------------------------------------------------------------------------------*/	
	/*---------------------------------------------------------------------------------*/
	
	dval(x , y , w , sigma , v , d , Nx , Ny);
	
	/*-----------------------------------------------*/
	/*-----------------------------------------------*/
	/* ------------ END of Mex File ---------------- */
	/*-----------------------------------------------*/
	/*-----------------------------------------------*/

	if ((nrhs < 3) || mxIsEmpty(prhs[2]))
	{			
		mxFree(w);
	}
}

/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/

void dval(double *x , double *y , double *w , double sigma , double *v , int d , int Nx , int Ny )
{
	int i , j , l , id , jd;	
	register double temp , res ;
	double tempv;
	double cte = -1.0/(sigma*sigma);
	int num_threads;
	
#ifdef OMP 
    num_threads          = min(MAX_THREADS,omp_get_num_procs());
    omp_set_num_threads(num_threads);
#endif

	for (i = 0; i < Ny ; i++) 
	{
		id    = i*d;	
		tempv = 0.0;
#ifdef OMP 
/* #pragma omp parallel for default(none) private(j,l,temp,jd,res) shared(x,y,w,d,id,Nx,cte,tempv) */
#pragma omp parallel for default(none) firstprivate(l,temp,res) lastprivate(j,jd,tempv) shared(x,y,w,d,id,Nx,cte)
#endif
		for (j = 0 ; j < Nx ; j++) 
		{
			jd  = j*d;
			res = 0.0;
			for (l = 0 ; l < d ; l++)
			{
				temp  = (y[l + id] - x[l + jd]);	
				res  +=temp*temp;
			}
			tempv += w[j]*exp(cte*res);	
		}
		v[i] = tempv;
	}
}
