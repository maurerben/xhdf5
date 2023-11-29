cmake_minimum_required(VERSION 3.5)

project(xhdf5 VERSION 0.0.0 DESCRIPTION "Library for a modern Fortran HDF5 interface")
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-ffree-form -std=f2018 -fimplicit-none -ffree-line-length-0 -cpp")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f18 -free -implicitnone -O3 -ip -unroll -scalar-rep -cpp")
    set(bounds "-check bounds")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")


find_package(HDF5 REQUIRED COMPONENTS Fortran HL)

add_compile_options(-D_HDF5_)

#if(${HDF5_F})


add_subdirectory(lib)
#add_subdirectory(test)



