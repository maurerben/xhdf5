# hdf5fort

`hdf5fort` is a Fortran HDF5 utility library with MPI-aware test coverage. It provides wrappers and utilities for easier HDF5 file, group, dataset, and datatype handling in Fortran.

## Repository layout

- `CMakeLists.txt` — root project build configuration
- `src/` — library sources and wrappers
- `test/` — test executables and CTest integration
- `build_parallel/` — recommended out-of-source build directory

## Features

- HDF5 Fortran support
- Optional MPI support
- Unit tests for initialization, group/link operations, scalar/vector/matrix I/O
- CTest-based automated validation

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

## Notes

- MPI test executables are launched via `${MPIEXEC_EXECUTABLE}` when `ENABLE_MPI` is enabled.
- The test suite is configured to use `min(3, nproc)` MPI ranks by default.
- `ctest -j` controls the number of test executables run in parallel, not the number of MPI ranks inside each test.
