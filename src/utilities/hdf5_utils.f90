!> Utilities for hdf5 wrapper
module hdf5_utils
   use iso_fortran_env, only: real32, real64, error_unit

#ifdef xhdf5
   use hdf5
#endif

   use mpi_utils, only: terminate
   use iso_fortran_env, only: real32, real64, int32

   implicit none

   !> Root group of any HDF5 file
   character(*), parameter :: hdf5_root = './'

   !> HDF5 integer kinds
#ifdef xhdf5
   integer, parameter :: hdf5_id = HID_T
   integer, parameter :: hdf5_size = HSIZE_T
   integer, parameter :: hdf5_ssize = HSSIZE_T
#else
   integer, parameter :: hdf5_id = real32
   integer, parameter :: hdf5_size = real32
   integer, parameter :: hdf5_ssize = real32
#endif


   !> Max. string length for write and read
   integer, parameter :: hdf5_max_string_len = 1024
   !> Blank, to fill non needed character length
   character(1), parameter :: hdf5_blank = ' '

contains


   !> Return HDF5 id for double
   integer(hdf5_id) function hdf5_double()
#ifdef xhdf5
      hdf5_double = H5T_NATIVE_DOUBLE
#else
      hdf5_double = real64
#endif
   end function

   !> Return HDF5 id for real (sp)
   integer(hdf5_id) function hdf5_float()
#ifdef xhdf5
      hdf5_float = H5T_NATIVE_REAL
#else
      hdf5_float = real32
#endif
   end function hdf5_float


   !> Return HDF5 id for integer
   integer(hdf5_id) function hdf5_integer_int32()
#ifdef xhdf5
      hdf5_integer_int32 = H5T_NATIVE_INTEGER
#else
      hdf5_integer_int32 = int32
#endif
   end function hdf5_integer_int32

   !> Return HDF5 id for character
   integer(hdf5_id) function hdf5_character()
#ifdef xhdf5
      hdf5_character = H5T_NATIVE_CHARACTER
#else
      hdf5_character = int32
#endif
   end function hdf5_character


   !> Assert if the HDF5 error flag `h5err` is zero or not. If it is not zero, the
   !> routine prints out a message that tells the routine name and the value of `h5err`.
   subroutine handle_hdf5_error(hdf5_routine_name, h5err, unit_out)
      !> Name of the hdf5 routine.
      character(*), intent(in) :: hdf5_routine_name
      !> HDF5 error flag.
      integer, intent(in) :: h5err
      !> Unit to write to. By default `error_unit`.
      integer, intent(in), optional :: unit_out

      integer :: unit_out_local

      unit_out_local = error_unit
      if (present(unit_out)) unit_out_local = unit_out

      if (h5err /= 0) then
         write(unit_out_local, '(A, A, A, I3)') 'HDF5 routine ', trim(adjustl(hdf5_routine_name)), ' failed with h5err = ', h5err
         call terminate()
      end if
   end subroutine handle_hdf5_error

end module hdf5_utils
