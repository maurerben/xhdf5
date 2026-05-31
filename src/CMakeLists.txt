# Module directories: only the public one (holding xhdf5.mod) is exposed to consumers.
set(xhdf5_PUBLIC_MODULE_DIR  "${CMAKE_BINARY_DIR}/modules")
set(xhdf5_PRIVATE_MODULE_DIR "${CMAKE_BINARY_DIR}/modules_internal")

# ---------------------------------------------------------------------------
# Internal modules (utilities + HDF5 wrappers).
# Compiled into an object library whose .mod files land in a private module
# directory that is never added to the consumers' include path.
# ---------------------------------------------------------------------------
file(GLOB_RECURSE xhdf5_INTERNAL_SOURCES
  "*.F90"
  "*.f90"
)

add_library(xhdf5_internal OBJECT ${xhdf5_INTERNAL_SOURCES})

set_target_properties(xhdf5_internal PROPERTIES
  Fortran_MODULE_DIRECTORY "${xhdf5_PRIVATE_MODULE_DIR}"
  Fortran_PREPROCESS ON
)

# HDF5 / MPI are required to compile the internal sources.
if(HDF5_FOUND)
  target_link_libraries(xhdf5_internal PUBLIC ${HDF5_LIBS})
  target_include_directories(xhdf5_internal PUBLIC ${HDF5_Fortran_INCLUDE_DIRS})
endif()
if(MPI_FOUND)
  target_link_libraries(xhdf5_internal PUBLIC MPI::MPI_Fortran)
endif()

# ---------------------------------------------------------------------------
# Public library: the generated `xhdf5` module plus the internal objects.
# ---------------------------------------------------------------------------
set(xhdf5_FYPP_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/xhdf5.F90.fypp")
set(xhdf5_GENERATED_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/xhdf5.F90")

fypp_generate(${xhdf5_FYPP_SOURCE} ${xhdf5_GENERATED_SOURCE})

add_library(xhdf5 STATIC ${xhdf5_GENERATED_SOURCE})

# xhdf5.mod is written to the public module directory.
set_target_properties(xhdf5 PROPERTIES
  Fortran_MODULE_DIRECTORY "${xhdf5_PUBLIC_MODULE_DIR}"
  Fortran_PREPROCESS ON
)

# Pull the internal objects into libxhdf5.a. PRIVATE so the internal module
# directory is never propagated to consumers.
target_link_libraries(xhdf5 PRIVATE xhdf5_internal)

# xhdf5.F90 `use`s the internal modules, so it needs their .mod files at build time.
target_include_directories(xhdf5 PRIVATE "${xhdf5_PRIVATE_MODULE_DIR}")

# HDF5 / MPI remain part of the public link interface for downstream targets.
if(HDF5_FOUND)
  target_link_libraries(xhdf5 PUBLIC ${HDF5_LIBS})
  target_include_directories(xhdf5 PUBLIC ${HDF5_Fortran_INCLUDE_DIRS})
endif()
if(MPI_FOUND)
  target_link_libraries(xhdf5 PUBLIC MPI::MPI_Fortran)
endif()

# Only the public module directory (containing xhdf5.mod) is exposed to consumers.
target_include_directories(xhdf5 PUBLIC
  $<BUILD_INTERFACE:${xhdf5_PUBLIC_MODULE_DIR}>
)
