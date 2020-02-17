# dixonfix

[![Build Status](https://travis-ci.org/biomedia-mira/dixonfix.svg?branch=master)](https://travis-ci.org/biomedia-mira/dixonfix)

*dixonfix* is a command line tool to correct fat-water swaps in Dixon MRI. It makes use of graph cut optimization to determine a pixel-wise swap labelling.

```
@inproceedings{glocker2016correction,
  title={Correction of fat-water swaps in Dixon MRI},
  author={Glocker, Ben and Konukoglu, Ender and Lavdas, Ioannis and Iglesias, Juan Eugenio and Aboagye, Eric O and Rockall, Andrea G and Rueckert, Daniel},
  booktitle={International Conference on Medical Image Computing and Computer-Assisted Intervention},
  pages={536--543},
  year={2016},
  organization={Springer}
}
```

If you make use of *dixonfix*, it would be great if you cite this paper in any resulting publications.

## Dependencies

*dixonfix* depends on several third-party libraries:

* [Eigen](eigen.tuxfamily.org)
* [Boost](http://www.boost.org/) (tested up to v1.58)
* [ITK](http://itk.org) (tested up to [v4.13.2](https://sourceforge.net/projects/itk/files/itk/4.13/InsightToolkit-4.13.2.tar.gz))

## Build instructions

Eigen is a header-only library and can be simply installed via:

```
#!bash

$ mkdir 3rdparty
$ cd 3rdparty
$ wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
$ mkdir -p eigen && tar xvf 3.3.7.tar.gz -C eigen --strip-components=1
```

You can download and install ITK in the same `3rdparty` folder via:

```
#!bash

$ wget https://sourceforge.net/projects/itk/files/itk/4.13/InsightToolkit-4.13.2.tar.gz
$ tar xvf InsightToolkit-4.13.2.tar.gz
$ cd InsightToolkit-4.13.2
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=../../itk ..
$ make -j4
$ make install
```

Alternatively, you can check out these [ITK install instructions](https://itk.org/Wiki/ITK/Getting_Started/Build/Linux).

You can install Boost and TBB via `apt-get`:

```
#!bash

$ sudo apt-get install libboost-all-dev
```

Note, you might have to specify a specific version via `apt-get install <package>=<version>`.

*dixonfix* comes with a CMake configuration file. From the top folder where `CMakeLists.txt` is located (same as this README), do the following to build all internal libraries and executables:

```
#!bash

$ mkdir build
$ cd build
$ export THIRD_PARTY_DIR=<folder_containing_eigen_and_itk>
$ cmake ..
$ make -j4

```

## Usage

Run `./dixonfix -h` to see a list of command line arguments.
