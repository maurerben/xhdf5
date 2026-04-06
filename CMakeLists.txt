cmake_minimum_required(VERSION 3.15)
project(solhdf5 Fortran)

# Set Fortran standard
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

# Options
option(ENABLE_HDF5 "Enable HDF5 support" ON)
option(ENABLE_MPI "Enable MPI support" OFF)

# Find HDF5 if enabled
if(ENABLE_HDF5)
  find_package(HDF5 COMPONENTS Fortran)
  if(HDF5_FOUND)
    add_compile_definitions(_HDF5_)
    message(STATUS "HDF5 found and enabled")
    # Add the HDF5 definitions to the compiler options
    add_definitions(${HDF5_DEFINITIONS})
    # Include directories for HDF5 header files
    include_directories(${HDF5_INCLUDE_DIRS})
    # Link directories for HDF5 libraries
    link_directories(${HDF5_LIBRARY_DIRS})

    # Set the libraries to link against, combining C, C++, and Fortran libraries
    set(HDF5_LIBS "")
    if(DEFINED HDF5_C_LIBRARIES)
      list(APPEND HDF5_LIBS ${HDF5_C_LIBRARIES})
    endif()
    if(DEFINED HDF5_CXX_LIBRARIES)
      list(APPEND HDF5_LIBS ${HDF5_CXX_LIBRARIES})
    endif()
    if(DEFINED HDF5_Fortran_LIBRARIES)
      list(APPEND HDF5_LIBS ${HDF5_Fortran_LIBRARIES})
    endif()
    if(DEFINED HDF5_hdf5_LIBRARY_RELEASE)
      if(NOT HDF5_hdf5_LIBRARY_RELEASE STREQUAL "HDF5_hdf5_LIBRARY_RELEASE-NOTFOUND")
        list(APPEND HDF5_LIBS ${HDF5_hdf5_LIBRARY_RELEASE})
      endif()
    endif()
    if(DEFINED HDF5_hdf5_fortran_LIBRARY_RELEASE)
      if(NOT HDF5_hdf5_fortran_LIBRARY_RELEASE STREQUAL "HDF5_hdf5_fortran_LIBRARY_RELEASE-NOTFOUND")
        list(APPEND HDF5_LIBS ${HDF5_hdf5_fortran_LIBRARY_RELEASE})
      endif()
    endif()
    if(DEFINED HDF5_hdf5_LIBRARY_DEBUG)
      if(NOT HDF5_hdf5_LIBRARY_DEBUG STREQUAL "HDF5_hdf5_LIBRARY_DEBUG-NOTFOUND")
        list(APPEND HDF5_LIBS ${HDF5_hdf5_LIBRARY_DEBUG})
      endif()
    endif()
    if(DEFINED HDF5_hdf5_fortran_LIBRARY_DEBUG)
      if(NOT HDF5_hdf5_fortran_LIBRARY_DEBUG STREQUAL "HDF5_hdf5_fortran_LIBRARY_DEBUG-NOTFOUND")
        list(APPEND HDF5_LIBS ${HDF5_hdf5_fortran_LIBRARY_DEBUG})
      endif()
    endif()
  else()
    message(WARNING "HDF5 not found, building without HDF5 support")
    set(ENABLE_HDF5 OFF)
  endif()
endif()

# Find MPI if enabled
if(ENABLE_MPI)
  find_package(MPI COMPONENTS Fortran)
  if(MPI_FOUND)
    add_compile_definitions(MPI)
    message(STATUS "MPI found and enabled")
  else()
    message(WARNING "MPI not found, building without MPI support")
    set(ENABLE_MPI OFF)
  endif()
endif()

# Enable testing
enable_testing()

# Find FORD for documentation generation (optional)
find_program(FORD_EXECUTABLE ford)
if(FORD_EXECUTABLE)
  add_custom_target(doc
    COMMAND ${FORD_EXECUTABLE} ${CMAKE_SOURCE_DIR}/ford.md
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating documentation with FORD..."
  )
  message(STATUS "FORD found - use 'make doc' to generate documentation")
else()
  message(STATUS "FORD not found - documentation generation disabled (install with: pip install ford)")
endif()

# Add subdirectories
add_subdirectory(src)
add_subdirectory(test)
