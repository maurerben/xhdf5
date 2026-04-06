---
project: solhdf5
summary: Fortran HDF5 Utility Library
author: Benedikt Maurer
src_dir: ./src
output_dir: ./doc
exclude_dir: ./build_parallel
project_github: https://github.com/benediktmaurer/solhdf5
docmark: !>
predocmark: 
docmark_alt: !!
predocmark_alt: 
warn: true
quiet: false
print_creation_date: true
predocmark_alt: !!
docmark_alt: !!
---

# solhdf5 Documentation

Fortran HDF5 utility library with MPI-aware parallel I/O support. This library provides easy-to-use wrappers for HDF5 file operations, group management, and dataset I/O.

## Features

- **Easy HDF5 Wrappers**: Type-safe interface for HDF5 operations
- **Parallel Support**: Optional MPI support for parallel I/O
- **Multiple Data Types**: Support for integer, real, and complex scalars, vectors, and matrices
- **Group Management**: Create and manage HDF5 groups and hierarchical structures
- **Comprehensive Tests**: Full test suite with CTest integration

## Quick Start

Build with documentation:

```bash
cd build_parallel
cmake ..
make
ford ford.md
```

The documentation is generated in the `doc/` directory.

## Module Structure

The main module `solhdf5` provides:

- `h5file_t` type for file handling
- Generic interfaces for writing/reading different data types
- Group and link management routines
- Support for HDF5 hyperslabs for complex data distributions
