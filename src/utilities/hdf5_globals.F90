!> Global definitions for the solhdf5 module
module hdf5_globals
   use iso_fortran_env, only: int32

#ifdef _HDF5_
   use hdf5, only: HID_T, HSIZE_T, HSSIZE_T
#endif

   implicit none

   private

   !> Root group of any HDF5 file
   character(*), parameter, public :: h5group_root = './'

   !> HDF5 integer kinds
#ifdef _HDF5_
   integer, parameter, public :: hdf5_id = HID_T
   integer, parameter, public :: hdf5_size = HSIZE_T
   integer, parameter, public :: hdf5_ssize = HSSIZE_T
#else
   integer, parameter, public :: hdf5_id = int32
   integer, parameter, public :: hdf5_size = int32
   integer, parameter, public :: hdf5_ssize = int32
#endif

   !> Value for undefined integer.
   integer(int32), parameter, public :: int32_undefined = -1
   !> Value for undefined HDF5 id.
   integer(hdf5_id), parameter, public  :: file_id_undefined = 0
   !> Value for undefined HDF5 size.
   integer(hdf5_size), parameter, public :: h5size_undefined = -1

   !> Blank, to fill non needed character length.
   character(1), parameter, public :: hdf5_blank = ' '
   !> String to represent a `.true.` value.
   character(*), parameter, public :: true_string = 'true'
   !> String to represent a `.false.` value.
   character(*), parameter, public :: false_string = 'false'

end module hdf5_globals
