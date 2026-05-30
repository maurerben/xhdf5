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

# Minimal test executable
add_executable(minimal_hdf5_test src/minimal_test.F90)
target_link_libraries(minimal_hdf5_test PRIVATE xhdf5)

# Set include directories for tests
target_include_directories(minimal_hdf5_test PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add test
if(ENABLE_MPI)
  add_test(NAME minimal_hdf5_test COMMAND ${MPIEXEC} -np ${NUM_PROCESSES} $<TARGET_FILE:minimal_hdf5_test>)
else()
  add_test(NAME minimal_hdf5_test COMMAND minimal_hdf5_test)
endif()


# Unit tests for init/delete file routines
add_executable(init_delete_file src/init_delete_file.F90)
target_link_libraries(init_delete_file PRIVATE xhdf5)

# Set include directories for tests
target_include_directories(init_delete_file PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add tests for initialize/delete group and checking links routins 
if(ENABLE_MPI)
  add_test(NAME init_delete_file COMMAND ${MPIEXEC_EXECUTABLE} -np ${NUM_PROCESSES} $<TARGET_FILE:init_delete_file>)
else()
  add_test(NAME init_delete_file COMMAND init_delete_file)
endif()


# Unit tests for initialize/finalize routines
add_executable(init_group_path_exists_delete_link src/init_group_path_exists_delete_link.F90)
target_link_libraries(init_group_path_exists_delete_link PRIVATE xhdf5)

# Set include directories for tests
target_include_directories(init_group_path_exists_delete_link PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add test
if(ENABLE_MPI)
  add_test(NAME init_group_path_exists_delete_link COMMAND ${MPIEXEC_EXECUTABLE} -np ${NUM_PROCESSES} $<TARGET_FILE:init_group_path_exists_delete_link>)
else()
  add_test(NAME init_group_path_exists_delete_link COMMAND init_group_path_exists_delete_link)
endif()


# Unit tests for write/read data routines (strings and booleans)
add_executable(write_read_data src/write_read_data.F90)
target_link_libraries(write_read_data PRIVATE xhdf5)

# Set include directories for tests
target_include_directories(write_read_data PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add test
if(ENABLE_MPI)
  add_test(NAME write_read_data COMMAND ${MPIEXEC_EXECUTABLE} -np ${NUM_PROCESSES} $<TARGET_FILE:write_read_data>)
else()
  add_test(NAME write_read_data COMMAND write_read_data)
endif()


# Unit tests for write/read numerical scalar routines
add_executable(write_read_numerical_scalar src/write_read_numerical_scalar.F90)
target_link_libraries(write_read_numerical_scalar PRIVATE xhdf5)

# Set include directories for tests
target_include_directories(write_read_numerical_scalar PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add test
if(ENABLE_MPI)
  add_test(NAME write_read_numerical_scalar COMMAND ${MPIEXEC_EXECUTABLE} -np ${NUM_PROCESSES} $<TARGET_FILE:write_read_numerical_scalar>)
else()
  add_test(NAME write_read_numerical_scalar COMMAND write_read_numerical_scalar)
endif()


# Unit tests for write/read numerical vector routines
add_executable(write_read_numerical_vector src/write_read_numerical_vector.F90)
target_link_libraries(write_read_numerical_vector PRIVATE xhdf5)

# Set include directories for tests
target_include_directories(write_read_numerical_vector PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add test
if(ENABLE_MPI)
  add_test(NAME write_read_numerical_vector COMMAND ${MPIEXEC_EXECUTABLE} -np ${NUM_PROCESSES} $<TARGET_FILE:write_read_numerical_vector>)
else()
  add_test(NAME write_read_numerical_vector COMMAND write_read_numerical_vector)
endif()


# Unit tests for write/read numerical matrix routines
add_executable(write_read_numerical_matrix src/write_read_numerical_matrix.F90)
target_link_libraries(write_read_numerical_matrix PRIVATE xhdf5)

# Set include directories for tests
target_include_directories(write_read_numerical_matrix PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add test
if(ENABLE_MPI)
  add_test(NAME write_read_numerical_matrix COMMAND ${MPIEXEC_EXECUTABLE} -np ${NUM_PROCESSES} $<TARGET_FILE:write_read_numerical_matrix>)
else()
  add_test(NAME write_read_numerical_matrix COMMAND write_read_numerical_matrix)
endif()
