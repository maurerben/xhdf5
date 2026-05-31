# Get number of processors for MPI tests
if(NOT DEFINED NUM_PROCESSES)
  execute_process(COMMAND nproc OUTPUT_VARIABLE NPROC OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NPROC)
    if(NPROC GREATER 3)
      set(NUM_PROCESSES 3)
    else()
      set(NUM_PROCESSES ${NPROC})
    endif()
    message(STATUS "Setting NUM_PROCESSES to ${NUM_PROCESSES} (min(3, ${NPROC}))")
  else()
    set(NUM_PROCESSES 1)
    message(STATUS "nproc command failed, setting NUM_PROCESSES to 1")
  endif()
endif()

# Register a unit-test executable and its ctest entry.
#
# Usage:
#   add_xhdf5_test(<name> <source>)         # plain Fortran source
#   add_xhdf5_test(<name> <source.fypp> FYPP)  # generate the source from a fypp template first
function(add_xhdf5_test name source)
  cmake_parse_arguments(ARG "FYPP" "" "" ${ARGN})

  if(ARG_FYPP)
    set(generated "${CMAKE_CURRENT_BINARY_DIR}/${name}.F90")
    fypp_generate("${CMAKE_CURRENT_SOURCE_DIR}/${source}" "${generated}")
    set(source "${generated}")
  else()
    set(source "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
  endif()

  add_executable(${name} "${source}")
  target_link_libraries(${name} PRIVATE xhdf5)
  target_include_directories(${name} PRIVATE
    ${CMAKE_BINARY_DIR}/modules
    ${CMAKE_BINARY_DIR}/modules_internal
  )

  if(ENABLE_MPI)
    add_test(NAME ${name}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_PROCESSES} $<TARGET_FILE:${name}>)
  else()
    add_test(NAME ${name} COMMAND ${name})
  endif()
endfunction()

# Select which tests are built based on HDF5 availability. The test sources no
# longer guard themselves with #ifdef _HDF5_; instead the HDF5-requiring tests
# are only built when HDF5 is enabled, and the no-HDF5 behaviour test only when
# it is not.
if(ENABLE_HDF5)
  add_xhdf5_test(minimal_hdf5_test                  src/minimal_test.F90)
  add_xhdf5_test(init_delete_file                   src/init_delete_file.F90)
  add_xhdf5_test(init_group_path_exists_delete_link src/init_group_path_exists_delete_link.F90)
  add_xhdf5_test(write_read_data                    src/write_read_data.F90)
  add_xhdf5_test(write_read_numerical_scalar        src/write_read_numerical_scalar.F90.fypp FYPP)
  add_xhdf5_test(write_read_numerical_vector        src/write_read_numerical_vector.F90.fypp FYPP)
  add_xhdf5_test(write_read_numerical_matrix        src/write_read_numerical_matrix.F90.fypp FYPP)
else()
  # Without HDF5 the library routines must remain callable but do nothing / return false.
  add_xhdf5_test(no_hdf5_behavior src/no_hdf5_behavior.F90)
endif()
