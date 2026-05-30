!> Unit tests for write and read numerical matrix (rank 2 array) routines
program write_read_numerical_matrix_test
#ifdef MPI
   use mpi_f08
#endif
   use iso_fortran_env, only: output_unit, real32, real64, int32
   use xhdf5
   use mpi_utils

   implicit none

   integer :: num_tests = 0
   integer :: num_passed = 0
   integer :: num_failed = 0
#ifdef MPI
   integer :: mpi_err, comm_size, my_rank
#endif

#ifdef MPI
   ! Initialize MPI
   call mpi_init(mpi_err)
   call mpi_comm_size(MPI_COMM_WORLD, comm_size, mpi_err)
   call mpi_comm_rank(MPI_COMM_WORLD, my_rank, mpi_err)
#endif

   write(output_unit, '(a)') "=========================================="
   write(output_unit, '(a)') "  xhdf5 Numerical Matrix Unit Tests"
   write(output_unit, '(a)') "=========================================="
   write(output_unit, *)

#ifdef _HDF5_
   ! Test integer(int32) matrices
   call test_write_read_integer_int32_matrix()
   call test_write_read_integer_int32_matrix_serial()

   ! Test real(real32) matrices
   call test_write_read_real_real32_matrix()
   call test_write_read_real_real32_matrix_serial()

   ! Test real(real64) matrices
   call test_write_read_real_real64_matrix()
   call test_write_read_real_real64_matrix_serial()

   ! Test complex(real32) matrices
   call test_write_read_complex_real32_matrix()
   call test_write_read_complex_real32_matrix_serial()

   ! Test complex(real64) matrices
   call test_write_read_complex_real64_matrix()
   call test_write_read_complex_real64_matrix_serial()

#ifdef MPI
   ! Test row-distributed matrices
   call test_write_read_integer_int32_matrix_row_distributed()

   !  ! Test column-distributed matrices
   call test_write_read_integer_int32_matrix_col_distributed()

   !  ! Test block-distributed matrices
   call test_write_read_integer_int32_matrix_block_distributed()
#endif

#else
   write(output_unit, '(a)') "HDF5 not available - tests skipped"
#endif

   ! Summary
#ifdef MPI
   write(output_unit, *)
   write(output_unit, '(a)') "=========================================="
   write(output_unit, '(a,i0,a,i0,a,i0,a,i0,a)') "Results: ", num_passed, " passed, ", &
      num_failed, " failed out of ", num_tests, " tests (MPI_COMM_WORLD size: ", comm_size, ")"
   write(output_unit, '(a)') "=========================================="
#else
   write(output_unit, *)
   write(output_unit, '(a)') "=========================================="
   write(output_unit, '(a,i0,a,i0,a,i0)') "Results: ", num_passed, " passed, ", &
      num_failed, " failed out of ", num_tests, " tests"
   write(output_unit, '(a)') "=========================================="
#endif

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

   ! Test for integer(int32) matrix - row-wise parallel distribution
   subroutine test_write_read_integer_int32_matrix()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read integer(int32) matrix (row-distributed, parallel mode)"
      character(*), parameter :: filename = "test_integer_int32_matrix_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: n_rows = 20, n_cols = 10
      integer(int32), allocatable :: local_matrix(:,:), read_matrix(:,:)
      integer(int32) :: local_rows, start_row, i, j
#ifdef MPI
      integer(int32) :: datastart(2), datasize(2)
#endif

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
#ifdef MPI
         ! Row-wise distribution: each rank gets n_rows/comm_size rows
         local_rows = n_rows / comm_size
         start_row = my_rank * local_rows + 1

         allocate(local_matrix(local_rows, n_cols))
         allocate(read_matrix(local_rows, n_cols))

         ! Fill matrix with rank-specific data
         do j = 1, n_cols
            do i = 1, local_rows
               local_matrix(i, j) = (start_row + i - 2) * n_cols + j
            end do
         end do

         ! Set hyperslab parameters
         datastart(1) = start_row
         datastart(2) = 1
         datasize(1) = n_rows
         datasize(2) = n_cols

         call h5file%write(h5path, dataset, local_matrix, datastart=datastart, datasize=datasize)
         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " wrote rows ", start_row, " to ", start_row + local_rows - 1

         call h5file%read(h5path, dataset, read_matrix, datastart=datastart)

         if (all(read_matrix == local_matrix)) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix)
#else
         ! Single process mode - write and read the full matrix
         allocate(local_matrix(n_rows, n_cols))
         allocate(read_matrix(n_rows, n_cols))
         do j = 1, n_cols
            do i = 1, n_rows
               local_matrix(i, j) = (i - 1) * n_cols + j
            end do
         end do

         call h5file%write(h5path, dataset, local_matrix)
         write(output_unit, '(a)') "    Full matrix written successfully"

         call h5file%read(h5path, dataset, read_matrix)

         if (all(read_matrix == local_matrix)) then
            write(output_unit, '(a)') "    Matrix read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Matrix read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix)
#endif
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_integer_int32_matrix

   subroutine test_write_read_integer_int32_matrix_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read integer(int32) matrix (serial mode)"
      character(*), parameter :: filename = "test_integer_int32_matrix_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: n_rows = 20, n_cols = 10
      integer(int32), allocatable :: full_matrix(:,:), read_matrix(:,:)
      integer(int32) :: i, j

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call init_with_mpi_comm_world(h5file%mpi_comm)

      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            allocate(full_matrix(n_rows, n_cols))
            allocate(read_matrix(n_rows, n_cols))
            do j = 1, n_cols
               do i = 1, n_rows
                  full_matrix(i, j) = (i - 1) * n_cols + j
               end do
            end do

            call h5file%write(h5path, dataset, full_matrix)
            write(output_unit, '(a)') "    Full matrix written successfully"

            call h5file%read(h5path, dataset, read_matrix)

            if (all(read_matrix == full_matrix)) then
               write(output_unit, '(a)') "    Matrix read matches written values"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Matrix read mismatch"
               num_failed = num_failed + 1
            end if

            deallocate(full_matrix, read_matrix)
            call h5file%delete()
         else
            write(output_unit, '(a)') "    File initialization failed"
            num_failed = num_failed + 1
         end if
      else
         write(output_unit, '(a)') "    Non-root process - skipping serial operations"
         num_passed = num_passed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_integer_int32_matrix_serial

   ! Test for real(real32) matrix
   subroutine test_write_read_real_real32_matrix()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real32) matrix (row-distributed, parallel mode)"
      character(*), parameter :: filename = "test_real_real32_matrix_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: n_rows = 20, n_cols = 10
      real(real32), allocatable :: local_matrix(:,:), read_matrix(:,:)
      integer(int32) :: local_rows, start_row, i, j
#ifdef MPI
      integer(int32) :: datastart(2), datasize(2)
#endif

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
#ifdef MPI
         local_rows = n_rows / comm_size
         start_row = my_rank * local_rows + 1

         allocate(local_matrix(local_rows, n_cols))
         allocate(read_matrix(local_rows, n_cols))

         do j = 1, n_cols
            do i = 1, local_rows
               local_matrix(i, j) = real((start_row + i - 2) * n_cols + j, real32) + 0.5_real32
            end do
         end do

         datastart(1) = start_row
         datastart(2) = 1
         datasize(1) = n_rows
         datasize(2) = n_cols

         call h5file%write(h5path, dataset, local_matrix, datastart=datastart, datasize=datasize)
         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " wrote rows ", start_row, " to ", start_row + local_rows - 1

         call h5file%read(h5path, dataset, read_matrix, datastart=datastart)

         if (all(abs(read_matrix - local_matrix) < epsilon(local_matrix))) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix)
#else
         allocate(local_matrix(n_rows, n_cols))
         allocate(read_matrix(n_rows, n_cols))
         do j = 1, n_cols
            do i = 1, n_rows
               local_matrix(i, j) = real((i - 1) * n_cols + j, real32) + 0.5_real32
            end do
         end do

         call h5file%write(h5path, dataset, local_matrix)
         write(output_unit, '(a)') "    Full matrix written successfully"

         call h5file%read(h5path, dataset, read_matrix)

         if (all(abs(read_matrix - local_matrix) < epsilon(local_matrix))) then
            write(output_unit, '(a)') "    Matrix read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Matrix read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix)
#endif
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_real_real32_matrix

   subroutine test_write_read_real_real32_matrix_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real32) matrix (serial mode)"
      character(*), parameter :: filename = "test_real_real32_matrix_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: n_rows = 20, n_cols = 10
      real(real32), allocatable :: full_matrix(:,:), read_matrix(:,:)
      integer(int32) :: i, j

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call init_with_mpi_comm_world(h5file%mpi_comm)

      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            allocate(full_matrix(n_rows, n_cols))
            allocate(read_matrix(n_rows, n_cols))
            do j = 1, n_cols
               do i = 1, n_rows
                  full_matrix(i, j) = real((i - 1) * n_cols + j, real32) + 0.5_real32
               end do
            end do

            call h5file%write(h5path, dataset, full_matrix)
            write(output_unit, '(a)') "    Full matrix written successfully"

            call h5file%read(h5path, dataset, read_matrix)

            if (all(abs(read_matrix - full_matrix) < epsilon(full_matrix))) then
               write(output_unit, '(a)') "    Matrix read matches written values"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Matrix read mismatch"
               num_failed = num_failed + 1
            end if

            deallocate(full_matrix, read_matrix)
            call h5file%delete()
         else
            write(output_unit, '(a)') "    File initialization failed"
            num_failed = num_failed + 1
         end if
      else
         write(output_unit, '(a)') "    Non-root process - skipping serial operations"
         num_passed = num_passed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_real_real32_matrix_serial

   ! Test for real(real64) matrix
   subroutine test_write_read_real_real64_matrix()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real64) matrix (row-distributed, parallel mode)"
      character(*), parameter :: filename = "test_real_real64_matrix_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: n_rows = 20, n_cols = 10
      real(real64), allocatable :: local_matrix(:,:), read_matrix(:,:)
      integer(int32) :: local_rows, start_row, i, j
#ifdef MPI
      integer(int32) :: datastart(2), datasize(2)
#endif

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
#ifdef MPI
         local_rows = n_rows / comm_size
         start_row = my_rank * local_rows + 1

         allocate(local_matrix(local_rows, n_cols))
         allocate(read_matrix(local_rows, n_cols))

         do j = 1, n_cols
            do i = 1, local_rows
               local_matrix(i, j) = real((start_row + i - 2) * n_cols + j, real64) + 0.5_real64
            end do
         end do

         datastart(1) = start_row
         datastart(2) = 1
         datasize(1) = n_rows
         datasize(2) = n_cols

         call h5file%write(h5path, dataset, local_matrix, datastart=datastart, datasize=datasize)
         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " wrote rows ", start_row, " to ", start_row + local_rows - 1

         call h5file%read(h5path, dataset, read_matrix, datastart=datastart)

         if (all(abs(read_matrix - local_matrix) < epsilon(local_matrix))) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix)
#else
         allocate(local_matrix(n_rows, n_cols))
         allocate(read_matrix(n_rows, n_cols))
         do j = 1, n_cols
            do i = 1, n_rows
               local_matrix(i, j) = real((i - 1) * n_cols + j, real64) + 0.5_real64
            end do
         end do

         call h5file%write(h5path, dataset, local_matrix)
         write(output_unit, '(a)') "    Full matrix written successfully"

         call h5file%read(h5path, dataset, read_matrix)

         if (all(abs(read_matrix - local_matrix) < epsilon(local_matrix))) then
            write(output_unit, '(a)') "    Matrix read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Matrix read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix)
#endif
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_real_real64_matrix

   subroutine test_write_read_real_real64_matrix_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real64) matrix (serial mode)"
      character(*), parameter :: filename = "test_real_real64_matrix_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: n_rows = 20, n_cols = 10
      real(real64), allocatable :: full_matrix(:,:), read_matrix(:,:)
      integer(int32) :: i, j

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call init_with_mpi_comm_world(h5file%mpi_comm)

      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            allocate(full_matrix(n_rows, n_cols))
            allocate(read_matrix(n_rows, n_cols))
            do j = 1, n_cols
               do i = 1, n_rows
                  full_matrix(i, j) = real((i - 1) * n_cols + j, real64) + 0.5_real64
               end do
            end do

            call h5file%write(h5path, dataset, full_matrix)
            write(output_unit, '(a)') "    Full matrix written successfully"

            call h5file%read(h5path, dataset, read_matrix)

            if (all(abs(read_matrix - full_matrix) < epsilon(full_matrix))) then
               write(output_unit, '(a)') "    Matrix read matches written values"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Matrix read mismatch"
               num_failed = num_failed + 1
            end if

            deallocate(full_matrix, read_matrix)
            call h5file%delete()
         else
            write(output_unit, '(a)') "    File initialization failed"
            num_failed = num_failed + 1
         end if
      else
         write(output_unit, '(a)') "    Non-root process - skipping serial operations"
         num_passed = num_passed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_real_real64_matrix_serial

   ! Test for complex(real32) matrix
   subroutine test_write_read_complex_real32_matrix()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real32) matrix (row-distributed, parallel mode)"
      character(*), parameter :: filename = "test_complex_real32_matrix_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: n_rows = 20, n_cols = 10
      complex(real32), allocatable :: local_matrix(:,:), read_matrix(:,:)
      integer(int32) :: local_rows, start_row, i, j
#ifdef MPI
      integer(int32) :: datastart(2), datasize(2)
#endif

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
#ifdef MPI
         local_rows = n_rows / comm_size
         start_row = my_rank * local_rows + 1

         allocate(local_matrix(local_rows, n_cols))
         allocate(read_matrix(local_rows, n_cols))

         do j = 1, n_cols
            do i = 1, local_rows
               local_matrix(i, j) = cmplx(real((start_row + i - 2) * n_cols + j, real32), &
                  real((start_row + i - 2) * n_cols + j, real32) + 0.5, real32)
            end do
         end do

         datastart(1) = start_row
         datastart(2) = 1
         datasize(1) = n_rows
         datasize(2) = n_cols

         call h5file%write(h5path, dataset, local_matrix, datastart=datastart, datasize=datasize)
         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " wrote rows ", start_row, " to ", start_row + local_rows - 1

         call h5file%read(h5path, dataset, read_matrix, datastart=datastart)

         if (all(abs(read_matrix - local_matrix) < epsilon(real(local_matrix)))) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix)
#else
         allocate(local_matrix(n_rows, n_cols))
         allocate(read_matrix(n_rows, n_cols))
         do j = 1, n_cols
            do i = 1, n_rows
               local_matrix(i, j) = cmplx(real((i - 1) * n_cols + j, real32), &
                  real((i - 1) * n_cols + j, real32) + 0.5, real32)
            end do
         end do

         call h5file%write(h5path, dataset, local_matrix)
         write(output_unit, '(a)') "    Full matrix written successfully"

         call h5file%read(h5path, dataset, read_matrix)

         if (all(abs(read_matrix - local_matrix) < epsilon(real(local_matrix)))) then
            write(output_unit, '(a)') "    Matrix read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Matrix read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix)
#endif
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_complex_real32_matrix

   subroutine test_write_read_complex_real32_matrix_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real32) matrix (serial mode)"
      character(*), parameter :: filename = "test_complex_real32_matrix_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: n_rows = 20, n_cols = 10
      complex(real32), allocatable :: full_matrix(:,:), read_matrix(:,:)
      integer(int32) :: i, j

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call init_with_mpi_comm_world(h5file%mpi_comm)

      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            allocate(full_matrix(n_rows, n_cols))
            allocate(read_matrix(n_rows, n_cols))
            do j = 1, n_cols
               do i = 1, n_rows
                  full_matrix(i, j) = cmplx(real((i - 1) * n_cols + j, real32), &
                     real((i - 1) * n_cols + j, real32) + 0.5, real32)
               end do
            end do

            call h5file%write(h5path, dataset, full_matrix)
            write(output_unit, '(a)') "    Full matrix written successfully"

            call h5file%read(h5path, dataset, read_matrix)

            if (all(abs(read_matrix - full_matrix) < epsilon(real(full_matrix)))) then
               write(output_unit, '(a)') "    Matrix read matches written values"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Matrix read mismatch"
               num_failed = num_failed + 1
            end if

            deallocate(full_matrix, read_matrix)
            call h5file%delete()
         else
            write(output_unit, '(a)') "    File initialization failed"
            num_failed = num_failed + 1
         end if
      else
         write(output_unit, '(a)') "    Non-root process - skipping serial operations"
         num_passed = num_passed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_complex_real32_matrix_serial

   ! Test for complex(real64) matrix
   subroutine test_write_read_complex_real64_matrix()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real64) matrix (row-distributed, parallel mode)"
      character(*), parameter :: filename = "test_complex_real64_matrix_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: n_rows = 20, n_cols = 10
      complex(real64), allocatable :: local_matrix(:,:), read_matrix(:,:)
      integer(int32) :: local_rows, start_row, i, j
#ifdef MPI
      integer(int32) :: datastart(2), datasize(2)
#endif

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
#ifdef MPI
         local_rows = n_rows / comm_size
         start_row = my_rank * local_rows + 1

         allocate(local_matrix(local_rows, n_cols))
         allocate(read_matrix(local_rows, n_cols))

         do j = 1, n_cols
            do i = 1, local_rows
               local_matrix(i, j) = cmplx(real((start_row + i - 2) * n_cols + j, real64), &
                  real((start_row + i - 2) * n_cols + j, real64) + 0.5, real64)
            end do
         end do

         datastart(1) = start_row
         datastart(2) = 1
         datasize(1) = n_rows
         datasize(2) = n_cols

         call h5file%write(h5path, dataset, local_matrix, datastart=datastart, datasize=datasize)
         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " wrote rows ", start_row, " to ", start_row + local_rows - 1

         call h5file%read(h5path, dataset, read_matrix, datastart=datastart)

         if (all(abs(read_matrix - local_matrix) < epsilon(real(local_matrix)))) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix)
#else
         allocate(local_matrix(n_rows, n_cols))
         allocate(read_matrix(n_rows, n_cols))
         do j = 1, n_cols
            do i = 1, n_rows
               local_matrix(i, j) = cmplx(real((i - 1) * n_cols + j, real64), &
                  real((i - 1) * n_cols + j, real64) + 0.5, real64)
            end do
         end do

         call h5file%write(h5path, dataset, local_matrix)
         write(output_unit, '(a)') "    Full matrix written successfully"

         call h5file%read(h5path, dataset, read_matrix)

         if (all(abs(read_matrix - local_matrix) < epsilon(real(local_matrix)))) then
            write(output_unit, '(a)') "    Matrix read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Matrix read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix)
#endif
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_complex_real64_matrix

   subroutine test_write_read_complex_real64_matrix_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real64) matrix (serial mode)"
      character(*), parameter :: filename = "test_complex_real64_matrix_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: n_rows = 20, n_cols = 10
      complex(real64), allocatable :: full_matrix(:,:), read_matrix(:,:)
      integer(int32) :: i, j

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call init_with_mpi_comm_world(h5file%mpi_comm)

      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            allocate(full_matrix(n_rows, n_cols))
            allocate(read_matrix(n_rows, n_cols))
            do j = 1, n_cols
               do i = 1, n_rows
                  full_matrix(i, j) = cmplx(real((i - 1) * n_cols + j, real64), &
                     real((i - 1) * n_cols + j, real64) + 0.5, real64)
               end do
            end do

            call h5file%write(h5path, dataset, full_matrix)
            write(output_unit, '(a)') "    Full matrix written successfully"

            call h5file%read(h5path, dataset, read_matrix)

            if (all(abs(read_matrix - full_matrix) < epsilon(real(full_matrix)))) then
               write(output_unit, '(a)') "    Matrix read matches written values"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Matrix read mismatch"
               num_failed = num_failed + 1
            end if

            deallocate(full_matrix, read_matrix)
            call h5file%delete()
         else
            write(output_unit, '(a)') "    File initialization failed"
            num_failed = num_failed + 1
         end if
      else
         write(output_unit, '(a)') "    Non-root process - skipping serial operations"
         num_passed = num_passed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_complex_real64_matrix_serial

   ! Test for row-distributed matrix
   subroutine test_write_read_integer_int32_matrix_row_distributed()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read integer(int32) matrix (row-distributed)"
      character(*), parameter :: filename = "test_integer_int32_matrix_row_distributed.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: rows_per_rank = 8, n_cols = 6
      integer(int32) :: n_rows
      integer(int32), allocatable :: local_matrix(:,:), read_matrix(:,:)
      integer(int32) :: local_rows, start_row, i, j
      type(hyperslab_type), allocatable :: hyperslabs(:)
      integer(int32) :: datastart(2), datasize(2)

#ifdef MPI
      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Total matrix size depends on number of processes
         n_rows = comm_size * rows_per_rank

         ! Simple uniform row-wise distribution
         local_rows = rows_per_rank
         start_row = my_rank * local_rows + 1

         allocate(local_matrix(local_rows, n_cols))
         allocate(read_matrix(local_rows, n_cols))
         allocate(hyperslabs(1))

         ! Fill matrix data
         do j = 1, n_cols
            do i = 1, local_rows
               local_matrix(i, j) = (start_row + i - 2) * n_cols + j
            end do
         end do

         ! Set up hyperslab
         datastart(1) = start_row
         datastart(2) = 1
         datasize(1) = n_rows
         datasize(2) = n_cols

         call hyperslabs(1)%init(shape(local_matrix, kind=int32), datasize, [-1], datastart, [-1], [-1], [-1], .false.)

         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " writing rows ", start_row, " to ", start_row + local_rows - 1

         call h5file%write(h5path, dataset, local_matrix, hyperslabs=hyperslabs)
         call h5file%read(h5path, dataset, read_matrix, hyperslabs=hyperslabs)

         if (all(read_matrix == local_matrix)) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix, hyperslabs)
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if
#else
      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name
      write(output_unit, '(a)') "    Single process mode - skipping row distribution test"
      num_passed = num_passed + 1
#endif

      write(output_unit, *)
   end subroutine test_write_read_integer_int32_matrix_row_distributed

   ! Test for column-distributed matrix
   subroutine test_write_read_integer_int32_matrix_col_distributed()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read integer(int32) matrix (column-distributed)"
      character(*), parameter :: filename = "test_integer_int32_matrix_col_distributed.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: n_rows = 10, cols_per_rank = 5
      integer(int32) :: n_cols
      integer(int32), allocatable :: local_matrix(:,:), read_matrix(:,:)
      integer(int32) :: local_cols, start_col, i, j
      type(hyperslab_type), allocatable :: hyperslabs(:)
      integer(int32) :: datastart(2), datasize(2)

#ifdef MPI
      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Total matrix size depends on number of processes
         n_cols = comm_size * cols_per_rank

         ! Simple uniform column-wise distribution
         local_cols = cols_per_rank
         start_col = my_rank * local_cols + 1

         allocate(local_matrix(n_rows, local_cols))
         allocate(read_matrix(n_rows, local_cols))
         allocate(hyperslabs(1))

         ! Fill matrix data
         do j = 1, local_cols
            do i = 1, n_rows
               local_matrix(i, j) = (i - 1) * n_cols + (start_col + j - 1)
            end do
         end do

         ! Set up hyperslab
         datastart(1) = 1
         datastart(2) = start_col
         datasize(1) = n_rows
         datasize(2) = n_cols

         call hyperslabs(1)%init(shape(local_matrix, kind=int32), datasize, [-1], datastart, [-1], [-1], [-1], .false.)

         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " writing columns ", start_col, " to ", start_col + local_cols - 1

         call h5file%write(h5path, dataset, local_matrix, hyperslabs=hyperslabs)
         call h5file%read(h5path, dataset, read_matrix, hyperslabs=hyperslabs)

         if (all(read_matrix == local_matrix)) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix, hyperslabs)
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if
#else
      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name
      write(output_unit, '(a)') "    Single process mode - skipping column distribution test"
      num_passed = num_passed + 1
#endif

      write(output_unit, *)
   end subroutine test_write_read_integer_int32_matrix_col_distributed

   ! Test for block-distributed matrix (2D domain decomposition)
   subroutine test_write_read_integer_int32_matrix_block_distributed()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read integer(int32) matrix (block-distributed 2x2 decomposition)"
      character(*), parameter :: filename = "test_integer_int32_matrix_block_distributed.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_matrix"
      integer(int32), parameter :: block_rows_per_rank = 4, block_cols_per_rank = 4
      integer(int32) :: n_rows, n_cols
      integer(int32), allocatable :: local_matrix(:,:), read_matrix(:,:)
      integer(int32) :: local_rows, local_cols, start_row, start_col, i, j
      integer(int32) :: block_row, block_col, grid_rows, grid_cols
      type(hyperslab_type), allocatable :: hyperslabs(:)
      integer(int32) :: datastart(2), datasize(2)

#ifdef MPI
      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Calculate grid dimensions for 2D decomposition
         ! Find factors of comm_size that are as square as possible
         grid_rows = int(sqrt(real(comm_size)))
         do while (grid_rows > 1 .and. mod(comm_size, grid_rows) /= 0)
            grid_rows = grid_rows - 1
         end do
         grid_cols = comm_size / grid_rows

         ! Total matrix size depends on number of processes and block size
         n_rows = grid_rows * block_rows_per_rank
         n_cols = grid_cols * block_cols_per_rank

         ! Calculate rank position in the grid
         block_row = my_rank / grid_cols
         block_col = mod(my_rank, grid_cols)

         ! Local block dimensions
         local_rows = block_rows_per_rank
         local_cols = block_cols_per_rank

         ! Starting position for this rank's block
         start_row = block_row * local_rows + 1
         start_col = block_col * local_cols + 1

         allocate(local_matrix(local_rows, local_cols))
         allocate(read_matrix(local_rows, local_cols))
         allocate(hyperslabs(1))

         ! Fill matrix data
         do j = 1, local_cols
            do i = 1, local_rows
               local_matrix(i, j) = (start_row + i - 2) * n_cols + (start_col + j - 1)
            end do
         end do

         ! Set up hyperslab
         datastart(1) = start_row
         datastart(2) = start_col
         datasize(1) = n_rows
         datasize(2) = n_cols

         call hyperslabs(1)%init(shape(local_matrix, kind=int32), datasize, [-1], datastart, [-1], [-1], [-1], .false.)

         write(output_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0)') "    Rank ", my_rank, " writing rows ", start_row, &
            " to ", start_row + local_rows - 1, " cols ", start_col, " to ", start_col + local_cols - 1

         call h5file%write(h5path, dataset, local_matrix, hyperslabs=hyperslabs)
         call h5file%read(h5path, dataset, read_matrix, hyperslabs=hyperslabs)

         if (all(read_matrix == local_matrix)) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_matrix, read_matrix, hyperslabs)
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if
#else
      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name
      write(output_unit, '(a)') "    Single process mode - skipping block distribution test"
      num_passed = num_passed + 1
#endif

      write(output_unit, *)
   end subroutine test_write_read_integer_int32_matrix_block_distributed

end program write_read_numerical_matrix_test
