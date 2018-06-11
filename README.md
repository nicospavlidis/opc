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

``` sh
unzip opc-master.zip
cd opc-master
```

3. Run `mexme_fgt` to compile Fast Gauss Transform and other C++ code.

``` matlab
cd('opc-master/src/libs/fgt/')
```

4. Test the installation. The included unit tests will test basic functionality to ensure that SnapVX and its dependencies are working correctly.

        cd Tests
        chmod u+x test_basic.sh
        ./test_basic.sh

Note: to run SnapVX locally, without installing it system-wide, just run setup.py with the --user flag.




#### Authors

* **Nicos Pavlidis** 

#### License

This project is licensed under the BSD-3-Clause License - see the [LICENSE.md](LICENSE.md) file for details

#### Acknowledgments

* David Hofmeyr
* Sotiris Tasoulis
* Hankui Peng
