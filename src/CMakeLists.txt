# Automatically find all Fortran source files
file(GLOB_RECURSE SOLHDF5_SOURCES
  "*.F90"
  "*.f90"
)

# Create the library
add_library(solhdf5 STATIC ${SOLHDF5_SOURCES})

# Set module output directory
set_target_properties(solhdf5 PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
)

# Enable automatic Fortran module dependency scanning
set_target_properties(solhdf5 PROPERTIES
  Fortran_PREPROCESS ON
)

# Link against HDF5 if it was found
if(HDF5_FOUND)
  target_link_libraries(solhdf5 PUBLIC ${HDF5_LIBS})
  target_include_directories(solhdf5 PUBLIC ${HDF5_Fortran_INCLUDE_DIRS})
endif()

# Link against MPI if it was found
if(MPI_FOUND)
  target_link_libraries(solhdf5 PUBLIC MPI::MPI_Fortran)
endif()

# Expose modules directory for other targets
target_include_directories(solhdf5 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
)
