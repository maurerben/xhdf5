cmake_minimum_required(VERSION 3.15)
project(xhdf5 Fortran)

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

# Locate fypp, used to generate Fortran sources from .fypp templates.
#
# Resolution order (first match wins):
#   1. An explicit fypp passed via -DFYPP_EXECUTABLE=/path/to/fypp
#      (also the way to use a vendored script shipped in the source tree)
#   2. A locally installed `fypp` executable on PATH
#   3. The `fypp` Python module (`python3 -m fypp`)
set(FYPP_COMMAND "")

if(FYPP_EXECUTABLE)
  # 1. Explicit user override (e.g. a vendored script) - used as-is.
  message(STATUS "Using user-provided fypp: ${FYPP_EXECUTABLE}")
else()
  # 2. Locally installed fypp executable.
  find_program(FYPP_EXECUTABLE NAMES fypp)
  if(FYPP_EXECUTABLE)
    message(STATUS "Using installed fypp: ${FYPP_EXECUTABLE}")
  endif()
endif()

if(NOT FYPP_EXECUTABLE)
  # 3. fypp provided as a Python module.
  find_package(Python3 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "import fypp"
      RESULT_VARIABLE FYPP_PYTHON_RESULT
      OUTPUT_QUIET
      ERROR_QUIET
    )
    if(FYPP_PYTHON_RESULT EQUAL 0)
      set(FYPP_EXECUTABLE ${Python3_EXECUTABLE})
      set(FYPP_COMMAND -m fypp)
      message(STATUS "Using fypp Python module: ${Python3_EXECUTABLE} -m fypp")
    endif()
  endif()
endif()

if(NOT FYPP_EXECUTABLE)
  message(FATAL_ERROR "fypp is required to build xhdf5. Install the 'fypp' executable or "
    "Python module, or point -DFYPP_EXECUTABLE=... at a vendored script.")
endif()

# Generate a Fortran source file from a .fypp template.
# Usage: fypp_generate(<input.fypp> <output.F90>)
function(fypp_generate FYPP_INPUT FYPP_OUTPUT)
  add_custom_command(
    OUTPUT ${FYPP_OUTPUT}
    COMMAND ${FYPP_EXECUTABLE} ${FYPP_COMMAND} ${FYPP_INPUT} ${FYPP_OUTPUT}
    DEPENDS ${FYPP_INPUT}
    COMMENT "Generating ${FYPP_OUTPUT} from fypp"
    VERBATIM
  )
endfunction()

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


# Default install prefix to a writable location inside the build tree.
# Only applies when the user did NOT pass -DCMAKE_INSTALL_PREFIX=...
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install"
      CACHE PATH "Install prefix" FORCE)
endif()
message(STATUS "xhdf5 install prefix: ${CMAKE_INSTALL_PREFIX}")


install(TARGETS xhdf5
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
  DESTINATION include
  FILES_MATCHING PATTERN "*.mod")
