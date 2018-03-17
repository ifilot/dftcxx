# DFTCXX

DFTCXX calculates the electronic structure of simple molecules within the framework of Density Functional Theory (DFT) at the LDA level of theory. It is mainly written for educational purposes. The source code has been documented (i.e. commented) relatively extensively to provide students the opportunity to read and understand the algorithm.

## Compilation

DFTCXX depends on a couple of libraries, which are normally directly available by your favorite package manager.

* Boost
* TCLAP
* Eigen3

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
