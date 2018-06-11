OPC: Optimal Projections for Clustering
=====================


Author: Nicos G. Pavlidis
E-mail: n(.)pavlidis(at)lancaster(.)ac(.)uk
Date:     2018-06-07

OPC is an open source package for MATLAB and GNU Octave that implements
clustering methods that seek the optimal low dimensional subspace to identify
clusters.


#### DEPENDENCIES

To install the package you need a C/C++ compiler. For Linux/ Mac we strongly
recommend using the GCC compiler. For Microsoft Windows the documentation
describes how to set up the MinGW-w64 compiler.

OPC depends on the following two packages, which are included in the package:

1. A cluster tree class, called ctree, which is a modification of the MATLAB class [tree](https://tinevez.github.io/matlab-tree/)
implemented by Jean-Yves Tinevez.

2. The [Fast Gauss Transform mex implementation](https://uk.mathworks.com/matlabcentral/fileexchange/17438-fast-gaussian-transform-mex-implementation?focused=5194134&tab=example)
(FGT) by Sebastien Paris.

To create a function reference in HTML format the 
[M2HTML](https://github.com/pdollar/toolbox/tree/master/external/m2html)
documentation system is required.

#### INSTALLATION

1. Download the latest release [here](https://github.com/nicospavlidis/opc/)

2. Uncompress the `opc-master.zip` file.

``` bash
unzip opc-master.zip
cd opc-master
```

3. From MATLAB/ Octave execute the `mexme_fgt` script to compile Fast Gauss Transform and other C++ code. This needs to be done
only once.

``` matlab
cd('DOWNLOAD-PATH/opc-master/src/libs/fgt/')
mexme_fgt
```

4. Each time that you restart MATLAB/ Octave you need to add to the search path the root OPC directory and all its subdirectories.
This is performed by the `setpath.m` script.

``` matlab
cd('DOWNLOAD-PATH/opc-master/')
setpath
```

5. Test installation: The script `reproduction_script.m` reproduces all the documentation material. If this script
exits without an error OPC is configured correctly.

``` matlab
reproduction_script
```

More detailed instructions for installing OPC and creating the HTML function reference
are provided in the file `documentation/document.pdf`.

### CONTACT

Please file bug reports at [github.com/nicospavlidis/opc](https://github.com/nicospavlidis/opc).
For any other questions, comments, or concerns, please contact [Nicos Pavlidis](http://www.stanford.edu/~hallac/).

#### LICENSE

This project is licensed under the BSD-3-Clause License - see the [LICENSE.md](LICENSE.md) file for details

#### ACKNOWLEDGMENTS

* David Hofmeyr
* Sotiris Tasoulis
* Michael Epitropakis
* Hankui Peng
