# xhdf5
**xhdf5** is a collection of modern Fortran wrappers to the HDF5 library. 
It allows for a straight forward interface to access HDF5 files in modern Fortran projects. 
The wrappers support serial and parallel (MPI) file writing with any striding.

## Supported features:
- Manage HDF5 groups
- Write and read datasets for **bool** and **string** data
- Write and read datasets of any rank for **integer**, **float**, **double**, **complex float** and **complex double**
- Serial file access
- Parallel (MPI) file access

# Install
The library can be easily installed with cmake (version >= 3.5)

## Supported compilers:
The library has been tested for **Intel comilers** (version >= 2021) and **GNU compilers** (version >= 8).

## Required libraries:
To install **xhdf5**, the **HDF5 library** is needed (version >= 1.12.0). 

There is an option to install the library without **HDF5**. In this case, calling the subroutines will have no effect.
