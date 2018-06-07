# OPC: Optimal Projections for Clustering

Author: Nicos G. Pavlidis
E-mail: n(.)pavlidis(at)lancaster(.)ac(.)uk
Date:     2018-06-07

=============================================
Optimal Projections for Clustering
=============================================

OPC is a MATLAB and GNU Octave package that implements a number of
dimensionality reduction for clustering algorithms.

-----------------------------------------------------------------------
INSTALLATION INSTRUCTIONS
-----------------------------------------------------------------------

0) OPC the MATLAB Optimization and Statistics Toolboxes and the MATLAB compiler. 

========================================================================================================
1.1) INSTALLING MINGW-W64 COMPILER IN WINDOWS (prior to compiling FGT)

In Windows you can install the MinGW-w64 compiler from the MATLAB add-ons menu see:

https://uk.mathworks.com/help/releases/R2015b/matlab/matlab_external/install-mingw-support-package.html

If this does not work you can still install this compiler manually and use it in MATLAB. 

+ Download the TDM-GCC MinGW installer:

  https://sourceforge.net/projects/tdm-gcc/

  The most recent version is currently version 5.1.0-2.

  https://sourceforge.net/projects/tdm-gcc/files/TDM-GCC%20Installer/tdm64-gcc-5.1.0-2.exe/download

+ After running the installer Mathworks recommends disabling the option "Check
  for updated files on the TDM-GCC server" that appears on the first screen
  after executing the installer. 

+ In the first screen of the installer select "Create" to create a new
  installation and ensure that you install the compiler in the directory:

  C:\TDM-GCC-64

  and not in "Program Files" or any directory whose path contains spaces, as suggested by Mathworks:

  https://uk.mathworks.com/help/releases/R2015b/matlab/matlab_external/install-mingw-support-package.html

+ After installation MATLAB will still not be able to find the MingGW-w64
  compiler because the path is not set correctly. To fix this use the setenv()
  function to set the environment variable MW_MINGW64_LOC to the directory that
  the compiler is installed:

>> setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')
>> mex -setup
MEX configured to use 'MinGW64 Compiler (C)' for C language compilation.

+ You will need to set this environment variable after every restart of MATLAB
  if you wish to compile MEX functions
========================================================================================================


Once the mex is correctly configured to use a C and C++ compiler:

+ Change directory to 'libs/fgt/' and run the mexme_fgt script to compile FGT. Then change
directory to 'src/mdh' and execute the mdh_compile script.

>> cd libs/fgt/
>> mexme_fgt
>> cd ../../src/mdh/
>> mdh_compile

The last script compiles simple C++ implementations of univariate kernel
density estimators using the Gaussian kernel, kdeC.cpp, and its derivative,
dkdeDxC.cpp. 

2) The code uses (a slightly modified version of) the 'tree' data structure
MATLAB implementation by Jean-Yves Tinevez. This is included in the
directory 'libs/@tree' and does not require compilation. For the original
library and documentation refer to:

http://tinevez.github.io/matlab-tree/index.html

3) The code uses the CBREWER function in MATLAB. A copy of this is included
in the directory 'libs/cbrewer'. For more information refer to:

https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab

-----------------------------------------------------------------------
EXAMPLES
-----------------------------------------------------------------------

Sample code for using the library in Matlab is provided in 'examples'
directory. The examples are documented thoroughly in the accompanying paper. 


-----------------------------------------------------------------------
DATASETS
-----------------------------------------------------------------------

Datasets are located in the 'datasets/' directory. All filenames have
extension .mat and load an N-by-D data matrix (X), with observations stored in
rows, and an N-by1 vector (labels) containing the true cluster assignments.
