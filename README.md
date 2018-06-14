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

OPC depends on the following two open source packages, which are included in the package:


1. The [improved Fast Gauss Transform](http://legacydirs.umiacs.umd.edu/~morariu/figtree/).


2. A cluster tree class, called ctree, which is a modification of the MATLAB class [tree](https://tinevez.github.io/matlab-tree/)
implemented by Jean-Yves Tinevez.

To create a function reference in HTML format the 
[M2HTML](https://github.com/pdollar/toolbox/tree/master/external/m2html)
documentation system is required.

#### INSTALLATION

1. Download the latest OPC release from this page.

2. Uncompress the `opc-master.zip` file.

``` bash
unzip opc-master.zip
cd opc-master
```

3. Compile the C++ FIGTree library by following the instructions on 
the README.md file located in `opc-master/src/libs/figtree-0.9.3/`, or equivalently
on the [FIGTree GitHub repository](https://github.com/vmorariu/figtree). 

4. After the C++ code is compiled you need to compile the MATLAB/ Octave interface
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


5. After each restart of MATLAB/ Octave it is necessary to add to the search path the root OPC
   directory and all its subdirectories.  This is performed by the `setpath.m` script.

``` matlab
>> cd('DOWNLOAD-PATH/opc-master/')
>> setpath
```

6. Testing the installation: The script `reproduction_script.m` 
in the root OPC folder reproduces all the examples in the documentation. If this script
exits without an error OPC is configured correctly.

``` matlab
cd('DOWNLOAD-PATH/opc-master/')
reproduction_script
```

More detailed instructions for installing OPC are provided in the file `documentation/document.pdf`.

#### DOCUMENTATION

* The file  `document/documentation.pdf` contains a detailed documentation of OPC
including more detailed installation instructions, and a large number of
examples that illustrate how to use and extend the code. The source TeX file and all the
figures are also included in the `document/` directory.

* An HTML function reference for OPC can be constructed using the
[M2HTML](https://github.com/pdollar/toolbox/tree/master/external/m2html). For this
purpose a custom template and a modified `m2html.m` script are provided in the
folder `document/m2html_template`.

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
