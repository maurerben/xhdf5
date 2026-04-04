# Minimal test executable
add_executable(minimal_hdf5_test minimal_test.f90)
target_link_libraries(minimal_hdf5_test PRIVATE solhdf5)

# Set include directories for tests
target_include_directories(minimal_hdf5_test PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add test
add_test(NAME minimal_hdf5_test COMMAND minimal_hdf5_test)

# Unit tests for initialize/finalize routines
add_executable(solhdf5_file_utilities solhdf5_file_utilities.f90)
target_link_libraries(solhdf5_file_utilities PRIVATE solhdf5)

# Set include directories for tests
target_include_directories(solhdf5_file_utilities PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add test
add_test(NAME solhdf5_file_utilities COMMAND solhdf5_file_utilities)
