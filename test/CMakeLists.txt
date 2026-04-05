# Minimal test executable
add_executable(minimal_hdf5_test src/minimal_test.F90)
target_link_libraries(minimal_hdf5_test PRIVATE solhdf5)

# Set include directories for tests
target_include_directories(minimal_hdf5_test PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add test
add_test(NAME minimal_hdf5_test COMMAND minimal_hdf5_test)


# Unit tests for init/delete file routines
add_executable(init_delete_file src/init_delete_file.F90)
target_link_libraries(init_delete_file PRIVATE solhdf5)

# Set include directories for tests
target_include_directories(init_delete_file PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add tests for initialize/delete group and checking links routins 
add_test(NAME init_delete_file COMMAND init_delete_file)


# Unit tests for initialize/finalize routines
add_executable(init_group_path_exists_delete_link src/init_group_path_exists_delete_link.F90)
target_link_libraries(init_group_path_exists_delete_link PRIVATE solhdf5)

# Set include directories for tests
target_include_directories(init_group_path_exists_delete_link PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add test
add_test(NAME init_group_path_exists_delete_link COMMAND init_group_path_exists_delete_link)


# Unit tests for write/read string routines
add_executable(write_read_string src/write_read_string.F90)
target_link_libraries(write_read_string PRIVATE solhdf5)

# Set include directories for tests
target_include_directories(write_read_string PRIVATE
  ${CMAKE_BINARY_DIR}/modules
)

# Add test
add_test(NAME write_read_string COMMAND write_read_string)
