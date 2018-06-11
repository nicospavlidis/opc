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

#### INSTALLATION
1. Download the latest release [here](https://github.com/nicospavlidis/opc/)
2. Uncompress the `opc-*.zip` file.

``` bash
unzip opc-master.zip
cd opc-master
```

3. From MATLAB/ Octave execute the `mexme_fgt` script to compile Fast Gauss Transform and other C++ code. This needs to be done
only once

``` matlab
cd('DOWNLOAD-PATH/opc-master/src/libs/fgt/')
mexme_fgt
```

4. Each time that you restart MATLAB/ Octave you need to add to the search path the root OPC folder and all its subdirectories.
This is performed by the `setpath.m` script

``` matlab
cd(''DOWNLOAD-PATH/opc-master/')
setpath
```

5. Test installation: The script `reproduction_script.m` reproduces all the documentation material. If this script
exits without any errors OPC is configured correctly.

``` matlab
reproduction_script
```


#### Authors

* **Nicos Pavlidis** 

#### License

This project is licensed under the BSD-3-Clause License - see the [LICENSE.md](LICENSE.md) file for details

#### Acknowledgments

* David Hofmeyr
* Sotiris Tasoulis
* Hankui Peng
