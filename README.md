# xhdf5

`xhdf5` is a Fortran HDF5 utility library with MPI-aware test coverage. It provides wrappers and utilities for easier HDF5 file, group, dataset, and datatype handling in Fortran.

## Repository layout

- `CMakeLists.txt` — root project build configuration
- `src/` — library sources and wrappers
- `test/` — test executables and CTest integration
- `build_parallel/` — recommended out-of-source build directory

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

## Prerequisites

- CMake 3.15+
- Fortran compiler (gfortran, ifort, etc.)
- HDF5 library with Fortran bindings (for parallel support, ensure HDF5 is built with MPI)
- MPI library (optional, for parallel support)
- FORD (optional, for documentation)

## Quick start

1. Create a build directory:

   ```bash
   mkdir -p build_parallel
   cd build_parallel
   ```

2. Configure the project with CMake:

   ```bash
   cmake ..
   ```

3. Build the library and tests:

   ```bash
   make -j$(nproc)
   ```

4. Run the full test suite:

   ```bash
   ctest -j 2
   ```

### Basic Usage

```fortran
use xhdf5

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
use xhdf5
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

## Notes

- MPI test executables are launched via `${MPIEXEC_EXECUTABLE}` when `ENABLE_MPI` is enabled.
- The test suite is configured to use `min(3, nproc)` MPI ranks by default.
- `ctest -j` controls the number of test executables run in parallel, not the number of MPI ranks inside each test.

## Documentation

API documentation is generated using [FORD](https://github.com/Fortran-FOSS-Programmers/ford). To generate docs:

```bash
ford ford.md
```

Documentation output is in `doc/index.html`.
