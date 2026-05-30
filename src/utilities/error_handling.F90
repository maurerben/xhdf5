! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Benedikt Maurer

!> Error handling for xhdf5.
module error_handling

   use mpi_utils, only: terminate, mpi_comm_type

   implicit none

   private
   public :: assert_true, handle_hdf5_error


contains


   !> Assert if a condition is true. If not the all processes of the MPI comminicator will be terminated.
   subroutine assert_true(condition, message)
      !> Condition that must be true to pass the assert_true.
      logical, intent(in) :: condition
      !> Message to print if assertion fails
      character(*), intent(in) :: message

      if (.not. condition) then
         call terminate(message=message)
      end if

   end subroutine assert_true


   !> Assert that the HDF5 error flag `h5err` is zero.
   !> If not, call [[terminate]] with an error message containing the name of `calling_routine`.
   subroutine handle_hdf5_error(mpi_comm, calling_routine, h5err)
      !> MPI communicator.
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Name of the hdf5 routine.
      character(*), intent(in) :: calling_routine
      !> HDF5 error flag.
      integer, intent(in) :: h5err

      character(3) :: h5err_char
      character(:), allocatable :: message

      if (h5err /= 0) then
         write(h5err_char, '(I3)') h5err
         message = 'HDF5 routine ' // calling_routine // ' failed with h5err = '// h5err_char
         call terminate(mpi_comm, message)
      end if
   end subroutine handle_hdf5_error

end module error_handling
