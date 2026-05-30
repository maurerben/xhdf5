! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Benedikt Maurer

!> Minimal test program - just open and close an HDF5 file
program minimal_hdf5_test
#ifdef MPI
   use mpi_f08
#endif
   use iso_fortran_env, only: output_unit
   use xhdf5

   implicit none

#ifdef MPI
   integer :: mpi_err
#endif
   type(h5file_t) :: h5file

#ifdef MPI
   call mpi_init(mpi_err)
#endif

   write(output_unit, '(a)') "Starting minimal HDF5 test..."

#ifdef _HDF5_
#ifdef MPI
   ! Open/Create file with MPI_COMM_WORLD
   write(output_unit, '(a)') "Opening HDF5 file with MPI_COMM_WORLD..."
#endif
   call h5file%init("test_minimal.h5", serial_access=.false.)

   if (h5file%is_open()) then
      write(output_unit, '(a)') "  File opened successfully"
   else
      write(output_unit, '(a)') "  Failed to open file"
   end if

   ! Close file
   write(output_unit, '(a)') "Closing HDF5 file..."
   call h5file%delete()

   if (h5file%file_id == 0) then
      write(output_unit, '(a)') "  File closed successfully"
   else
      write(output_unit, '(a)') "  Failed to close file"
   end if

#ifdef MPI
   call mpi_finalize(mpi_err)
#endif

   write(output_unit, '(a)') "  Test completed successfully"
#else
   write(output_unit, '(a)') "HDF5 not available - test skipped"
#endif

end program minimal_hdf5_test
