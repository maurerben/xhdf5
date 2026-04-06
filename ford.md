---
project: solhdf5
summary: Fortran HDF5 Utility Library with MPI Support
author: Benedikt Maurer
output_dir: ./doc
exclude_dir: ./build*
project_github: https://github.com/benediktmaurer/solhdf5
predocmark: >  
docmark: !  
display: public  
graph: false  
warn: false  
search: false  
proc_internals: false  
src_dir: ./src/
         ./src/hdf5_wrappers
         ./src/utitilities
         ./src/templates 

[//]: # "Note, ford commands can not be separated by whitelines."  
[//]: # "More information on ford's project file options can be found at:"  
[//]: # "https://github.com/Fortran-FOSS-Programmers/ford/wiki/Project-File-Options"  

---

# solhdf5 Documentation

## Overview

**solhdf5** is a modern Fortran library that provides easy-to-use, type-safe wrappers for HDF5 file operations with comprehensive MPI parallel I/O support. The library is designed to simplify HDF5 usage in scientific computing applications while maintaining high performance and flexibility.

## Key Features

### **Easy HDF5 Wrappers**

- Type-safe interface for all HDF5 operations
- Automatic memory management and error handling
- Intuitive Fortran-style API

### **Parallel MPI Support**

- Full MPI-aware parallel I/O operations
- Flexible hyperslab management for complex data distributions
- Configurable communicator support

### **Comprehensive Data Type Support**

- **Scalars**: Integer (int32), Real (real32/real64), Complex (complex32/complex64), Strings, Booleans
- **Arrays**: 1D, 2D, and higher-dimensional arrays with automatic shape handling
- **Complex Data**: Native support for complex numbers with proper HDF5 storage

### **Advanced Group Management**

- Hierarchical group creation and navigation
- Link existence checking and deletion
- Dataset shape querying and metadata access

### **Robust Testing**

- Complete test suite with CTest integration
- MPI scaling tests (configurable process counts)
- File I/O validation across different data types and distributions

## Quick Start

### Basic Usage

```fortran
use solhdf5

type(h5file_t) :: h5file
real(real64), allocatable :: buffer(:,:)

! Initialize HDF5 file
call h5file%init("output.h5")

! Write 2D matrix
call h5file%write("group", "dataset", buffer)

! Read buffer back
call h5file%read("group", "dataset", buffer)

! Clean up
call h5file%delete()
```

### Parallel MPI Usage

```fortran
use solhdf5
use mpi_f08

type(h5file_t) :: h5file
integer, parameter :: global_size(2) = [100, 100]
integer :: buffer_start(2)
real(real64), allocatable :: local_buffer(:,:)


! Initialize with MPI communicator
call h5file%init("parallel_output.h5", MPI_COMM_WORLD)

! Write distributed data
! Obtain buffer_start and local_buffer based on MPI rank and global_size
call h5file%write("/", "distributed_data", local_data, datastart=buffer_start, datasize=global_size)

call h5file%delete()
```

## Architecture

The library consists of several key modules:

- **`solhdf5`**: Main user interface module with `h5file_t` type
- **`hdf5_file`**: Low-level file operations (create/open/close)
- **`hdf5_group`**: Group and link management
- **`hdf5_dataset`**: Dataset creation and metadata
- **`hdf5_datatype`**: HDF5 datatype mapping
- **`hdf5_write`/`hdf5_read`**: I/O operations with hyperslab support
- **`hyperslab`**: Complex data distribution utilities
- **`mpi_utils`**: MPI communicator and rank management

## Build System

### Prerequisites

- CMake 3.15+
- Fortran compiler (gfortran, ifort, etc.)
- HDF5 library with Fortran bindings
- MPI library (optional, for parallel support)
- FORD (optional, for documentation)

### Build Commands

```bash
# Configure
cmake .. -DENABLE_MPI=ON

# Build
make -j$(nproc)

# Test
ctest -j 2

# Generate docs (if FORD installed)
make doc
```

## API Reference

The main interface is provided through the `h5file_t` derived type, which supports:

### File Operations

- `init()`: Initialize/open HDF5 file
- `delete()`: Close file and finalize HDF5

### Data I/O

- `write()`: Generic write interface for all supported data types
- `read()`: Generic read interface for all supported data types

### Group Management

- `init_group()`: Create new groups
- `link_exists()`: Check link existence
- `delete_link()`: Remove links/datasets

### Metadata

- `dataset_shape()`: Query dataset dimensions

## Examples

See the `test/` directory for comprehensive usage examples covering:

- Scalar data I/O
- Vector operations with different distributions
- Matrix I/O with MPI scaling
- Group and link management
- Error handling patterns

## Contributing

Contributions are welcome! Please ensure:

- All tests pass with `ctest`
- Code follows Fortran 2008 standards
- Documentation is updated for new features
- MPI scaling is tested for parallel features

## License

MIT License - see LICENSE file for details.

## Acknowledgments

This library builds upon the HDF5 library and MPI standards, providing a modern Fortran interface for scientific data management.

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
