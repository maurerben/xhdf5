# INSTALL

## Prerequisites

### Build Requirements

- CMake 3.15 or newer
- Fortran compiler (e.g. `gfortran`, `ifort`, `pgfortran`)
- HDF5 with Fortran support
- MPI libraries and compiler wrappers if using MPI support
- [fypp](https://github.com/aradi/fypp) preprocessor, used to generate Fortran sources from `.fypp` templates (install with `pip install fypp`)
- Linux utility `nproc` is used by the test CMake configuration to set default MPI ranks

### Documentation (optional)

- FORD (Fortran Automatic Documentation) for API documentation generation
  - Install with: `pip install ford`

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
- `-DFYPP_EXECUTABLE=<path>` to use a specific fypp instead of the auto-detected one

## Specifying fypp

The build generates Fortran sources from `.fypp` templates and needs a fypp preprocessor. By default CMake locates one automatically, in order: a `fypp` executable on `PATH`, then the `fypp` Python module (`python3 -m fypp`). If neither is found, configuration fails with an error.

To use a particular fypp instead of relying on auto-detection, point `FYPP_EXECUTABLE` at it. This is also how you use a fypp script vendored in your source tree (it is executed directly, so it must be executable and carry a `#!/usr/bin/env python3` shebang, as the upstream `bin/fypp` does):

```bash
cmake .. -DFYPP_EXECUTABLE=/path/to/fypp
```

## Generating API Documentation

To generate HTML documentation with FORD (requires FORD to be installed):

```bash
ford ford.md
```

The documentation output will be placed in `doc/` directory. Open `doc/index.html` in a browser to view the API documentation.

## Example full workflow

```bash
mkdir -p build_parallel
cd build_parallel
cmake .. -DENABLE_MPI=ON
make -j$(nproc)
ctest -j 2

# Generate documentation
ford ../ford.md
```
