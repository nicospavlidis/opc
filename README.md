## OPC: Optimal Projections for Clustering
### An open source MATLAB/Octave library for Dimensionality Reduction for Clustering

Author: [Nicos G. Pavlidis](http://www.lancaster.ac.uk/lums/people/nicos-pavlidis)<br/>
<!-- E-mail: n(.)pavlidis(at)lancaster(.)ac(.)uk<br/> -->
Date:     2018-07-11


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
	- `ldakmeans` Linear Discriminant Analysis k-means
		-  C. Ding and T. Li.
		[Adaptive dimension reduction using discriminant analysis and k-means](https://icml.cc/imls/conferences/2007/proceedings/papers/521.pdf).
		*Proceedings of the 24th International Conference on Machine Learning*, pages 521-528, 2007.

	- `drsc` Dimension Reduction for Spectral Clustering
		- D. Niu, J. G. Dy, and M. I. Jordan.
		[Dimensionality reduction for spectral clustering](http://proceedings.mlr.press/v15/niu11a/niu11a.pdf).
		*Proceedings of the 14th International Conference on Artificial Intelligence and Statistics*,
		pages 552-560, 2011.

		- D. Niu, J. G. Dy, and M. I. Jordan.
		[Iterative discovery of multiple alternative clustering views](https://doi.org/10.1109/TPAMI.2013.180).
		*IEEE Transactions on Pattern Analysis and Machine Intelligence*, 36(7):1340-1353, 2014.

	- `kurtClust` Kurtosis Projection Pursuit Clustering
		- D. Pena and F. J. Prieto.
		[Cluster Identification Using Projections](https://doi.org/10.1198/016214501753382345).
		*Journal of the American Statistical Association*, 96(456):1433-1445, 2001.

- **Divisive hierarchical clustering algorithms**

	- `pddp` Principal direction divisive partitioning
		- D. Boley. [Principal direction divisive partitioning](https://doi.org/10.1023/A:1009740529316).
		*Data Mining and Knowledge Discovery*, 2(4):325-344, 1998.

	- `depddp` Density-enhanced principal direction divisive partitioning
		- S. K. Tasoulis, D. K. Tasoulis, and V. P. Plagianakos.
		[Enhancing principal direction divisive clustering](https://doi.org/10.1016/j.patcog.2010.05.025).
		*Pattern Recognition*, 43(10):3391-3411, 2010.

	- `mddc` Minimum density divisive clustering
		- N. G. Pavlidis, D. P. Hofmeyr, and S. K. Tasoulis.
		[Minimum density hyperplanes](http://jmlr.csail.mit.edu/papers/v17/15-307.html).
		*Journal of Machine Learning Research*, 17(156):1-33, 2016.

	- `ncutdc` Minimum normalised cut divisive clustering
		- D. Hofmeyr. [Clustering by minimum cut hyperplanes](https://doi.org/10.1109/TPAMI.2016.2609929).
		*IEEE Transactions on Pattern Analysis and Machine Intelligence*, 39(8):1547-1560, 2017.
	
	- `mcdc` Maximum clusterability divisive clustering
		- D. P. Hofmeyr and N. G. Pavlidis.
		[Maximum clusterability divisive clustering](https://doi.org/10.1109/SSCI.2015.116).
		*IEEE Symposium Series on Computational Intelligence*, pages 780-786, 2015.

	- `scppdc` Minimum spectral connectivity projection pursuit
		- D. P. Hofmeyr, N. G. Pavlidis and I. A. Eckley.
		[Minimum spectral connectivity projection pursuit: Divisive clustering using optimal projections for spectral clustering](https://doi.org/10.1007/s11222-018-9814-6).
		*Statistics and Computing*, online since 2018.

	- Bisecting k-means, and bisecting LDA-k-means as illustrations of generic projection pursuit divisive
	clustering framework (`gppdc`).
		- M. Steinbach, G. Karypis, and V. Kumar.
		[A comparison of document clustering techniques](http://glaros.dtc.umn.edu/gkhome/fetch/papers/docclusterKDDTMW00.pdf).
		*Workshop on text mining, ACM SIGKDD International
		  Conference on Knowledge Discovery in Databases (KDD)*, pages 525-526, 2000.


#### DEPENDENCIES

* To install the package you need a C/C++ compiler. We recommend using the GCC compiler. 

* OPC depends on the [improved Fast Gauss Transform](http://legacydirs.umiacs.umd.edu/~morariu/figtree/).

* OPC contains a cluster tree class, called ctree. This is a modification of the MATLAB 
class [tree](https://tinevez.github.io/matlab-tree/) implemented by Jean-Yves Tinevez.

* In MATLAB OPC depends on the optimization and statistics toolboxes.

* In GNU Octave OPC depends on the optim and statistics packages.

<!-- * In MATLAB it is possible to create a function reference in HTML format using the
[M2HTML](https://github.com/pdollar/toolbox/tree/master/external/m2html). Detailed instructions
are provided in the PDF documentation file, `documentation/document.pdf`. -->


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

```Matlab
>> cd('DOWNLOAD-PATH/opc-master/')
>> install
```
* In GNU Octave you also need to install and load the optim and statistics packages, which can be
found at the extra packages for GNU Octave [repository](https://octave.sourceforge.io/packages.php).

##### Setting the path

After each restart of MATLAB/ Octave it is necessary to add to the search path the root OPC
directory and all its subdirectories.  This is performed by the `setpath.m` script.

```Matlab
>> cd('DOWNLOAD-PATH/opc-master/')
>> setpath
```

#### TESTING

After installation you can execute the script `reproduction_script.m` (located in
the root OPC directory) to reproduce all the examples in the documentation (`documentation/documentation.pdf`). 
If this script exits without an error OPC is configured correctly.

```Matlab
cd('DOWNLOAD-PATH/opc-master/')
reproduction_script
```

More detailed instructions for installing OPC are provided in the file `user_guide.pdf`.

#### DOCUMENTATION

* The [online documentation guide](http://www.lancaster.ac.uk/staff/pavlidin/software/opc/documentation/document/onlinedoc.html)
  contains detailed installation instructions and a large number of nontrivial
  examples that illustrate how to how to use and extend OPC. 
  <br>
  A PDF version of the documentation is also available in the file `user_guide.pdf`.  

  <!-- The source TeX file and all the figures are also
  included in the `documentation/` directory.-->

* An online function reference for OPC is available
  [here](http://www.lancaster.ac.uk/staff/pavlidin/software/opc/documentation/reference/index.html).

* The script file `reproduction_script.m` located in the root OPC directory
  contains all the examples discussed in the documentation.


#### CONTRIBUTING

The preferred way to contribute to OPC is to fork the [main repository](https://github.com/nicospavlidis/opc)
on GitHub. Detailed instructions are available through the online
[Git documentation](https://git-scm.com/documentation). 

Please use GitHub issues to file bug reports and feature requests. 

#### LICENSE

This project is licensed under the BSD-3-Clause License - see the [LICENSE.md](LICENSE.md) file for details

#### ACKNOWLEDGMENTS

The following people contributed to OPC (in alphabetical order):

* Michael Epitropakis
* David Hofmeyr
* Dimitris Kostaras
* Hankui Peng
* Sotiris Tasoulis
