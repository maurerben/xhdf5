! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Benedikt Maurer

!> Unit test verifying the behaviour of xhdf5 when built without HDF5 support.
!> In that configuration the routines must remain callable but simply do nothing
!> (no file is created) or return `.false.`, so that code can be linked and run
!> without HDF5 without crashing.
program no_hdf5_behavior_test
#ifdef MPI
   use mpi_f08
#endif
   use iso_fortran_env, only: output_unit, int32
   use xhdf5

   implicit none

   integer :: num_tests = 0
   integer :: num_passed = 0
   integer :: num_failed = 0
#ifdef MPI
   integer :: mpi_err
#endif

#ifdef MPI
   ! Initialize MPI
   call mpi_init(mpi_err)
#endif

   write(output_unit, '(a)') "=========================================="
   write(output_unit, '(a)') "  xhdf5 No-HDF5 Behavior Unit Tests"
   write(output_unit, '(a)') "=========================================="
   write(output_unit, *)

   call test_no_hdf5_noop()

   ! Summary
   write(output_unit, *)
   write(output_unit, '(a)') "=========================================="
   write(output_unit, '(a,i0,a,i0,a,i0)') "Results: ", num_passed, " passed, ", &
      num_failed, " failed out of ", num_tests, " tests"
   write(output_unit, '(a)') "=========================================="

   ! Exit with appropriate code
   if (num_failed > 0) then
#ifdef MPI
      call mpi_finalize(mpi_err)
#endif
      stop 1
   end if

#ifdef MPI
   ! Finalize MPI
   call mpi_finalize(mpi_err)
#endif

contains

   !> Without HDF5 every routine must be callable but must not do anything: no file is
   !> created on disk, the query routines return `.false.`, and nothing aborts the program.
   subroutine test_no_hdf5_noop()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: routines are callable and do nothing without HDF5"
      character(*), parameter :: filename = "no_hdf5_behavior.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_dataset"
      integer(int32) :: scalar_value = 42_int32
      logical :: file_exists

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Initialising must neither open anything nor create a file on disk.
      call h5file%init(filename)

      ! The file must report as not open.
      if (h5file%is_open()) then
         write(output_unit, '(a)') "    FAILED: is_open() returned .true. without HDF5"
         num_failed = num_failed + 1
         return
      end if

      ! Queries must return .false. rather than crashing.
      if (h5file%link_exists(h5path)) then
         write(output_unit, '(a)') "    FAILED: link_exists() returned .true. without HDF5"
         num_failed = num_failed + 1
         return
      end if

      ! Writing must be a safe no-op.
      call h5file%write(h5path, dataset, scalar_value)

      ! Closing must be a safe no-op.
      call h5file%delete()

      ! No file must have been created on disk.
      inquire(file=filename, exist=file_exists)
      if (file_exists) then
         write(output_unit, '(a)') "    FAILED: a file was created on disk without HDF5"
         num_failed = num_failed + 1
         return
      end if

      write(output_unit, '(a)') "    Routines were callable, returned .false., and created no file"
      num_passed = num_passed + 1
   end subroutine test_no_hdf5_noop

end program no_hdf5_behavior_test
