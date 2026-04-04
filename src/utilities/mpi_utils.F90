!> MPI utilities for the HDF5 wrappers
module mpi_utils

#ifdef MPI
   use mpi_f08
#endif

   use iso_fortran_env, only: error_unit

   implicit none

   private
   public :: mpi_comm_type, init_with_mpi_comm_world, terminate, &
      comm_to_rank, init_with_mpi_fint_comm, comm_to_fint, mpi_info_null_fint

   !> Rank of the root process.
   integer, parameter, public :: root_rank = 0

   !> Error code to send to `mpi_abort`.
   integer, parameter, private :: error_code = 101

   !> Type for handling the MPI communicator.
   type mpi_comm_type
#ifdef MPI
      type(MPI_Comm) :: handle = MPI_COMM_NULL
#else
      integer :: handle = 0
#endif
   end type



contains

   !> Initialize with MPI_COMM_WORLD
   subroutine init_with_mpi_comm_world(this)
      type(mpi_comm_type) :: this
#ifdef MPI
      this%handle = MPI_COMM_WORLD
#endif
   end subroutine


   !> Initialize communicator from a Fortran-integer communicator handle.
   subroutine init_with_mpi_fint_comm(this, comm_fint)
      type(mpi_comm_type), intent(out) :: this
      integer, intent(in) :: comm_fint
#ifdef MPI
      ! Bridge legacy integer communicator handles into mpi_f08 communicator type.
      this%handle = transfer(comm_fint, this%handle)
#else
      this%handle = comm_fint
#endif
   end subroutine init_with_mpi_fint_comm


   !> Convert communicator to a Fortran-integer communicator handle.
   integer function comm_to_fint(mpi_comm)
      type(mpi_comm_type), intent(in) :: mpi_comm
#ifdef MPI
      ! Bridge mpi_f08 communicator type back to legacy integer handle.
      comm_to_fint = transfer(mpi_comm%handle, comm_to_fint)
#else
      comm_to_fint = mpi_comm%handle
#endif
   end function comm_to_fint


   !> Return MPI_INFO_NULL as Fortran-integer handle.
   integer function mpi_info_null_fint()
#ifdef MPI
      mpi_info_null_fint = transfer(MPI_INFO_NULL, mpi_info_null_fint)
#else
      mpi_info_null_fint = 0
#endif
   end function mpi_info_null_fint

   !> Return the rank an MPI communicator handle belongs to.
   !>
   !> Wrapper for `mpi_comm_rank`.
   integer function comm_to_rank(mpi_comm)
      !> MPI communicator handle
      type(mpi_comm_type), intent(in) :: mpi_comm

      integer :: mpierr

#ifdef MPI
      call mpi_comm_rank(mpi_comm%handle, comm_to_rank, mpierr)
      call handle_mpi_error(mpi_comm, 'mpi_comm_rank', mpierr)
#else
      comm_to_rank = root_rank
#endif

   end function comm_to_rank


   !> Terminate an MPI environment.
   subroutine terminate(mpi_comm, message)
      !> MPI communicator.
      type(mpi_comm_type), intent(in), optional :: mpi_comm
      !> Error message
      character(len=*), intent(in), optional :: message

      integer :: mpierr, calling_rank

#ifdef MPI
      if(present(mpi_comm)) then
         calling_rank = comm_to_rank(mpi_comm)
      else
         call mpi_comm_rank(MPI_COMM_WORLD, calling_rank, mpierr)
         call handle_mpi_error(mpi_comm, 'mpi_comm_rank', mpierr)
      end if
      if(calling_rank == root_rank) then
         if(present(message)) write(error_unit, *) trim(adjustl(message))
      end if
      call mpi_abort(MPI_COMM_WORLD, error_code, mpierr)
#else
      if(present(message)) write(error_unit, *) trim(adjustl(message))
      stop
#endif

   end subroutine terminate


   !> Assert that the MPI error flag `mpierr` is zero.
   !> When not, `terminate_mpi_comm` is called to terminate the MPI environment.
   subroutine handle_mpi_error(mpi_comm, calling_routine, mpierr)
      !> MPI communicator.
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Name of the hdf5 routine.
      character(*), intent(in) :: calling_routine
      !> MPI error flag.
      integer, intent(in) :: mpierr

      character(3) :: mpierr_char
      character(:), allocatable :: message

      if (mpierr /= 0) then
         write(mpierr_char, '(I3)') mpierr
         message = 'MPI routine ' // calling_routine // ' failed with mpierr = '// mpierr_char
         call terminate(mpi_comm, message)
      end if

   end subroutine handle_mpi_error


end module
