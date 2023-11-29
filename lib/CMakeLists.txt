file(GLOB_RECURSE sources ./*.f90 ./*.F90)
add_library(xhdf5 ${sources})
target_link_libraries(xhdf5 HDF5::HDF5)