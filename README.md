OPC: Optimal Projections for Clustering library in MATLAB/Octave
=====================

Author: Nicos G. Pavlidis

E-mail: n(.)pavlidis(at)lancaster(.)ac(.)uk

Date:     2018-06-07


#### INTRODUCTION

OPC is an open source MATLAB and GNU Octave package that implements clustering
methods that seek the optimal low dimensional subspace to identify clusters.

Whenever the data contains irrelevant features, or correlations among subsets
of features exist (which is typical in high-dimensional data), or when clusters
are defined in different subspaces, the spatial data structure becomes less
informative about the underlying clusters. Under these conditions clustering
algorithms need to simultaneously solve two interrelated problems: (i) identify
the subspace in which clusters can be distinguished, and (ii) associate
observations to clusters. 

OPC focuses on methods which seek low dimensional subspaces that are optimal
with respect to specific clustering criteria. This distinguishes the methods in
OPC from generic dimensionality reduction techniques that optimise objective
functions that are not related to any clustering criterion, and are therefore
not guaranteed to preserve the cluster structure.


#### <a name="alg"> ALGORITHMS IMPLEMENTED IN OPC </a>

- **Partitioning clustering algorithms** 
	- `ldakmeans`: Linear Discriminant Analysis k-means<br/>
		C. Ding and T. Li.
		[Adaptive dimension reduction using discriminant analysis and k-means](http://users.cs.fiu.edu/~taoli/tenure/Ding-Li-ICML2007.pdf).
		*Proceedings of the 24th International Conference on Machine Learning*, pages 521--528, 2007.



#### DEPENDENCIES

* To install the package you need a C/C++ compiler. For Linux/ Mac we strongly
recommend using the GCC compiler. For Microsoft Windows the documentation
describes how to set up the MinGW-w64 compiler.

* OPC depends on the [improved Fast Gauss Transform](http://legacydirs.umiacs.umd.edu/~morariu/figtree/).


* A cluster tree class, called ctree, which is a modification of the MATLAB class [tree](https://tinevez.github.io/matlab-tree/)
implemented by Jean-Yves Tinevez.

* In MATLAB OPC depends on the optimization and statistics toolboxes

* In GNU Octave OPC depends on the optim and statistics packages

* To create a function reference in HTML format the 
[M2HTML](https://github.com/pdollar/toolbox/tree/master/external/m2html)
documentation system is required.


#### INSTALLATION

* Download the latest OPC release from this page.

* Uncompress the `opc-master.zip` file.

``` bash
unzip opc-master.zip
cd opc-master
```

* Compile the C++ FIGTree library by following the instructions on 
the README.md file located in `opc-master/src/libs/figtree-0.9.3/`, or equivalently
on the [FIGTree GitHub repository](https://github.com/vmorariu/figtree). 

* After the C++ code is compiled you need to compile the MATLAB/ Octave interface
   to the FIGTree library as well as two C++ functions for kernel density estimation
   included in OPC. For ease
   of use the script `install.m` in the root OPC directory performs these
   tasks. (Note that the `install.m` script assumes that the FIGTree library is in the
   folder `DOWNLOAD-PATH/opc-master/src/libs/figtree-0.9.3/`. If this is modified you
   need to edit this script to provide the correct path)

``` matlab
>> cd('DOWNLOAD-PATH/opc-master/')
>> install
```
* In GNU Octave you also need to install and load the optim and statistics packages, which can be
found at the extra packages for GNU Octave [repository](https://octave.sourceforge.io/packages.php).

##### Setting the path

After each restart of MATLAB/ Octave it is necessary to add to the search path the root OPC
directory and all its subdirectories.  This is performed by the `setpath.m` script.

``` matlab
>> cd('DOWNLOAD-PATH/opc-master/')
>> setpath
```

#### TESTING

After installation you can execute the script `reproduction_script.m` (located in
the root OPC directory) to reproduce all the examples in the documentation (`documentation/documentation.pdf`). 
If this script exits without an error OPC is configured correctly.

``` matlab
cd('DOWNLOAD-PATH/opc-master/')
reproduction_script
```

More detailed instructions for installing OPC are provided in the file `documentation/document.pdf`.

#### DOCUMENTATION

* The file  `documentation/document.pdf` contains a detailed documentation of OPC
including more detailed installation instructions, and a large number of
examples that illustrate how to use and extend the code. The source TeX file and all the
figures are also included in the `documentation/` directory.

* An HTML function reference for OPC can be constructed using the
[M2HTML](https://github.com/pdollar/toolbox/tree/master/external/m2html). For this
purpose a custom template and a modified `m2html.m` script are provided in the
folder `documentation/m2html_template`.

* The script file `reproduction_script.m` located in the root OPC directory
contains all the examples discussed in the documentation.


#### CONTACT

Please file bug reports at [github.com/nicospavlidis/opc](https://github.com/nicospavlidis/opc/).
For any other questions, comments, or concerns, please contact [Nicos Pavlidis](http://www.lancaster.ac.uk/lums/people/nicos-pavlidis/)

#### LICENSE

This project is licensed under the BSD-3-Clause License - see the [LICENSE.md](LICENSE.md) file for details

#### ACKNOWLEDGMENTS

The following people contributed to OPC (in alphabetical order):

* Michael Epitropakis
* David Hofmeyr
* Dimitris Kostaras
* Hankui Peng
* Sotiris Tasoulis
