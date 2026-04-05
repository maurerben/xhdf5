# INSTALL

## Prerequisites

- CMake 3.15 or newer
- Fortran compiler (e.g. `gfortran`, `ifort`, `pgfortran`)
- HDF5 with Fortran support
- MPI libraries and compiler wrappers if using MPI support
- Linux utility `nproc` is used by the test CMake configuration to set default MPI ranks

## Build steps

1. Create and enter the build directory:

   ```bash
   mkdir -p build_parallel
   cd build_parallel
   ```

2. Configure the project:

   - Without MPI:

     ```bash
     cmake ..
     ```

   - With MPI enabled:

     ```bash
     cmake .. -DENABLE_MPI=ON
     ```

   If you want to explicitly choose the number of MPI ranks for tests, pass:

   ```bash
   cmake .. -DENABLE_MPI=ON -DNUM_PROCESSES=2
   ```

3. Build the library and test executables:

   ```bash
   make -j$(nproc)
   ```

## Running tests

From `build_parallel/`:

```bash
ctest -j 2
```

### Important behavior

- `ctest -j N` controls how many test programs run concurrently.
- Each MPI-enabled test still launches with `-np ${NUM_PROCESSES}` as configured in `test/CMakeLists.txt`.
- By default, `NUM_PROCESSES` is set to `min(3, nproc)`.

## Common options

- `-DENABLE_HDF5=OFF` to build without HDF5 support
- `-DENABLE_MPI=OFF` to build without MPI support
- `-DNUM_PROCESSES=<n>` to override the default MPI rank count used for tests

## Example full workflow

```bash
mkdir -p build_parallel
cd build_parallel
cmake .. -DENABLE_MPI=ON
make -j$(nproc)
ctest -j 2
```
