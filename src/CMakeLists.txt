# Automatically find all Fortran source files
file(GLOB_RECURSE xhdf5_SOURCES
  "*.F90"
  "*.f90"
)

set(xhdf5_FYPP_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/xhdf5.F90.fypp")
set(xhdf5_GENERATED_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/xhdf5.F90")

fypp_generate(${xhdf5_FYPP_SOURCE} ${xhdf5_GENERATED_SOURCE})

list(APPEND xhdf5_SOURCES ${xhdf5_GENERATED_SOURCE})

# Create the library
add_library(xhdf5 STATIC ${xhdf5_SOURCES})

# Set module output directory
set_target_properties(xhdf5 PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
)

# Enable automatic Fortran module dependency scanning
set_target_properties(xhdf5 PROPERTIES
  Fortran_PREPROCESS ON
)

# Link against HDF5 if it was found
if(HDF5_FOUND)
  target_link_libraries(xhdf5 PUBLIC ${HDF5_LIBS})
  target_include_directories(xhdf5 PUBLIC ${HDF5_Fortran_INCLUDE_DIRS})
endif()

# Link against MPI if it was found
if(MPI_FOUND)
  target_link_libraries(xhdf5 PUBLIC MPI::MPI_Fortran)
endif()

# Expose modules directory for other targets
target_include_directories(xhdf5 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
)
