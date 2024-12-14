# DFTCXX

[![Build](https://github.com/ifilot/dftcxx/actions/workflows/build.yml/badge.svg)](https://github.com/ifilot/dftcxx/actions/workflows/build.yml)

DFTCXX calculates the electronic structure of simple molecules within the
framework of Density Functional Theory (DFT) at the LDA level of theory. It is
mainly written for educational purposes. The source code has been documented
(i.e. commented) relatively extensively to provide students the opportunity to
read and understand the algorithm.

> [!TIP]  
>  * Want to know more about electronic structure calculations? Have a look at my
>    [free lecture book](https://ifilot.pages.tue.nl/elements-of-electronic-structure-theory/).
>  * Interested in a Python-based Density Functional Theory code? Have a look at [PyDFT](https://github.com/ifilot/pydft)

## Compilation

DFTCXX depends on a couple of libraries, which are normally directly available by your favorite package manager.

* [Boost](https://www.boost.org/)
* [TCLAP](https://tclap.sourceforge.net/)
* [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
* [libPNG](http://www.libpng.org/pub/png/libpng.html)

To ensure you have the right packages on a (Debian-type) of operating system,
you can run the following

```
sudo apt install build-essential cmake libboost-all-dev pkg-config libeigen3-dev \
libpng-dev libtclap-dev
```

To compile the program:
```
mkdir build
cd build
cmake ../src
make -j9
```

## Execution
```
./dftcxx -i ../molecules/h2.in
```
