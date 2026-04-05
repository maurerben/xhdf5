!> Unit tests for write and read numerical vector (rank 1 array) routines
program write_read_numerical_vector_test
#ifdef MPI
   use mpi_f08
#endif
   use iso_fortran_env, only: output_unit, real32, real64, int32
   use solhdf5
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
   write(output_unit, '(a)') "  solhdf5 Numerical Vector Unit Tests"
   write(output_unit, '(a)') "=========================================="
   write(output_unit, *)

#ifdef _HDF5_
   ! Test integer(int32) vectors
   call test_write_read_integer_int32_vector()
   call test_write_read_integer_int32_vector_serial()

   ! Test real(real32) vectors
   call test_write_read_real_real32_vector()
   call test_write_read_real_real32_vector_serial()

   ! Test real(real64) vectors
   call test_write_read_real_real64_vector()
   call test_write_read_real_real64_vector_serial()

   ! Test complex(real32) vectors
   call test_write_read_complex_real32_vector()
   call test_write_read_complex_real32_vector_serial()

   ! Test complex(real64) vectors
   call test_write_read_complex_real64_vector()
   call test_write_read_complex_real64_vector_serial()

#ifdef MPI   
   ! Test uneven data distribution
   call test_write_read_integer_int32_vector_uneven()

   ! Test non-contiguous data distribution
   call test_write_read_integer_int32_vector_noncontiguous()
#endif

#else
   write(output_unit, '(a)') "HDF5 not available - tests skipped"
#endif
   
   ! Summary
#ifdef MPI
   write(output_unit, *)
   write(output_unit, '(a)') "=========================================="
   write(output_unit, '(a,i0,a,i0,a,i0)') "Results: ", num_passed, " passed, ", &
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

   ! Test for integer(int32) vector
   subroutine test_write_read_integer_int32_vector()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read integer(int32) vector (parallel mode)"
      character(*), parameter :: filename = "test_integer_int32_vector_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), parameter :: elements_per_rank = 10
      integer(int32), parameter :: total_elements = 20  ! 2 ranks * 10 elements each
      integer(int32), allocatable :: full_vector(:), local_vector(:), read_vector(:)
      integer(int32) :: start_idx, end_idx, local_size, i
#ifdef MPI
      integer(int32) :: datastart(1), datasize(1)
#endif

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)
      print *, "File initialized with ID: ", h5file%file_id

      if (h5file%file_id /= 0) then
#ifdef MPI
         ! Set up data distribution: rank 0 gets 1-10, rank 1 gets 11-20, etc.
         start_idx = my_rank * elements_per_rank + 1
         end_idx = (my_rank + 1) * elements_per_rank
         local_size = elements_per_rank
         print *, "Rank ", my_rank, " will write elements ", start_idx, " to ", end_idx

         ! Create local data for this rank
         allocate(local_vector(local_size))
         allocate(read_vector(local_size))
         do i = 1, local_size
            local_vector(i) = start_idx + i - 1
         end do

         ! Set up hyperslab for this rank's portion
         datastart(1) = start_idx
         datasize(1) = local_size * comm_size  ! Total size of the dataset is total_elements, but we specify the full size for collective I/O
         print *, "Rank ", my_rank, " datastart: ", datastart(1), " datasize: ", datasize(1)

         ! Write local portion
         call h5file%write(h5path, dataset, local_vector, datastart=datastart, datasize=datasize)
         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " wrote elements ", start_idx, " to ", end_idx

         ! Read back local portion
         call h5file%read(h5path, dataset, read_vector, datastart=datastart)

         ! Check if values match
         if (all(read_vector == local_vector)) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_vector, read_vector)
#else
         ! Single process mode - write and read the full vector
         allocate(full_vector(total_elements))
         allocate(read_vector(total_elements))
         do i = 1, total_elements
            full_vector(i) = i
         end do

         call h5file%write(h5path, dataset, full_vector)
         write(output_unit, '(a)') "    Full vector written successfully"

         call h5file%read(h5path, dataset, read_vector)

         if (all(read_vector == full_vector)) then
            write(output_unit, '(a)') "    Vector read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Vector read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(full_vector, read_vector)
#endif
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_integer_int32_vector

   subroutine test_write_read_integer_int32_vector_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read integer(int32) vector (serial mode)"
      character(*), parameter :: filename = "test_integer_int32_vector_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), parameter :: total_elements = 20
      integer(int32), allocatable :: full_vector(:), read_vector(:)
      integer(int32) :: i

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Initialize MPI communicator first to check rank
      call init_with_mpi_comm_world(h5file%mpi_comm)

      ! Only root process performs serial HDF5 operations
      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         ! Init file in serial mode (only root process allowed)
         call h5file%init(filename, serial_access=.true.)

         if (h5file%file_id /= 0) then
            ! Create and write full vector
            allocate(full_vector(total_elements))
            allocate(read_vector(total_elements))
            do i = 1, total_elements
               full_vector(i) = i
            end do

            call h5file%write(h5path, dataset, full_vector)
            write(output_unit, '(a)') "    Full vector written successfully"

            call h5file%read(h5path, dataset, read_vector)

            if (all(read_vector == full_vector)) then
               write(output_unit, '(a)') "    Vector read matches written values"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Vector read mismatch"
               num_failed = num_failed + 1
            end if

            deallocate(full_vector, read_vector)
            call h5file%delete()
         else
            write(output_unit, '(a)') "    File initialization failed"
            num_failed = num_failed + 1
         end if
      else
         ! Non-root processes should not participate in the test
         write(output_unit, '(a)') "    Non-root process - skipping serial operations"
         num_passed = num_passed + 1  ! Count as passed since this is expected behavior
      end if

      write(output_unit, *)
   end subroutine test_write_read_integer_int32_vector_serial

   ! Test for real(real32) vector
   subroutine test_write_read_real_real32_vector()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real32) vector (parallel mode)"
      character(*), parameter :: filename = "test_real_real32_vector_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), parameter :: elements_per_rank = 10
      integer(int32), parameter :: total_elements = 20
      real(real32), allocatable :: local_vector(:), read_vector(:)
      integer(int32) :: start_idx, local_size, i
#ifdef MPI
      integer(int32) :: datastart(1), datasize(1)
#endif

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%file_id /= 0) then
#ifdef MPI
         start_idx = my_rank * elements_per_rank + 1
         local_size = elements_per_rank

         allocate(local_vector(local_size))
         allocate(read_vector(local_size))
         do i = 1, local_size
            local_vector(i) = real(start_idx + i - 1, real32) + 0.5_real32
         end do

         datastart(1) = start_idx
         datasize(1) = local_size * comm_size

         call h5file%write(h5path, dataset, local_vector, datastart=datastart, datasize=datasize)
         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " wrote elements ", start_idx, " to ", start_idx + local_size - 1

         call h5file%read(h5path, dataset, read_vector, datastart=datastart)

         if (all(abs(read_vector - local_vector) < epsilon(local_vector))) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_vector, read_vector)
#else
         allocate(local_vector(total_elements))
         allocate(read_vector(total_elements))
         do i = 1, total_elements
            local_vector(i) = real(i, real32) + 0.5_real32
         end do

         call h5file%write(h5path, dataset, local_vector)
         write(output_unit, '(a)') "    Full vector written successfully"

         call h5file%read(h5path, dataset, read_vector)

         if (all(abs(read_vector - local_vector) < epsilon(local_vector))) then
            write(output_unit, '(a)') "    Vector read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Vector read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_vector, read_vector)
#endif
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_real_real32_vector

   subroutine test_write_read_real_real32_vector_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real32) vector (serial mode)"
      character(*), parameter :: filename = "test_real_real32_vector_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), parameter :: total_elements = 20
      real(real32), allocatable :: full_vector(:), read_vector(:)
      integer(int32) :: i

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call init_with_mpi_comm_world(h5file%mpi_comm)

      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         call h5file%init(filename, serial_access=.true.)

         if (h5file%file_id /= 0) then
            allocate(full_vector(total_elements))
            allocate(read_vector(total_elements))
            do i = 1, total_elements
               full_vector(i) = real(i, real32) + 0.5_real32
            end do

            call h5file%write(h5path, dataset, full_vector)
            write(output_unit, '(a)') "    Full vector written successfully"

            call h5file%read(h5path, dataset, read_vector)

            if (all(abs(read_vector - full_vector) < epsilon(full_vector))) then
               write(output_unit, '(a)') "    Vector read matches written values"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Vector read mismatch"
               num_failed = num_failed + 1
            end if

            deallocate(full_vector, read_vector)
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
   end subroutine test_write_read_real_real32_vector_serial

   ! Test for real(real64) vector
   subroutine test_write_read_real_real64_vector()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real64) vector (parallel mode)"
      character(*), parameter :: filename = "test_real_real64_vector_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), parameter :: elements_per_rank = 10
      integer(int32), parameter :: total_elements = 20
      real(real64), allocatable :: local_vector(:), read_vector(:)
      integer(int32) :: start_idx, local_size, i
#ifdef MPI
      integer(int32) :: datastart(1), datasize(1)
#endif

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call h5file%init(filename, serial_access=.false.)

      if (h5file%file_id /= 0) then
#ifdef MPI
         start_idx = my_rank * elements_per_rank + 1
         local_size = elements_per_rank

         allocate(local_vector(local_size))
         allocate(read_vector(local_size))
         do i = 1, local_size
            local_vector(i) = real(start_idx + i - 1, real64) + 0.5_real64
         end do

         datastart(1) = start_idx
         datasize(1) = local_size * comm_size

         call h5file%write(h5path, dataset, local_vector, datastart=datastart, datasize=datasize)
         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " wrote elements ", start_idx, " to ", start_idx + local_size - 1

         call h5file%read(h5path, dataset, read_vector, datastart=datastart)

         if (all(abs(read_vector - local_vector) < epsilon(local_vector))) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_vector, read_vector)
#else
         allocate(local_vector(total_elements))
         allocate(read_vector(total_elements))
         do i = 1, total_elements
            local_vector(i) = real(i, real64) + 0.5_real64
         end do

         call h5file%write(h5path, dataset, local_vector)
         write(output_unit, '(a)') "    Full vector written successfully"

         call h5file%read(h5path, dataset, read_vector)

         if (all(abs(read_vector - local_vector) < epsilon(local_vector))) then
            write(output_unit, '(a)') "    Vector read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Vector read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_vector, read_vector)
#endif
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_real_real64_vector

   subroutine test_write_read_real_real64_vector_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real64) vector (serial mode)"
      character(*), parameter :: filename = "test_real_real64_vector_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), parameter :: total_elements = 20
      real(real64), allocatable :: full_vector(:), read_vector(:)
      integer(int32) :: i

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call init_with_mpi_comm_world(h5file%mpi_comm)

      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         call h5file%init(filename, serial_access=.true.)

         if (h5file%file_id /= 0) then
            allocate(full_vector(total_elements))
            allocate(read_vector(total_elements))
            do i = 1, total_elements
               full_vector(i) = real(i, real64) + 0.5_real64
            end do

            call h5file%write(h5path, dataset, full_vector)
            write(output_unit, '(a)') "    Full vector written successfully"

            call h5file%read(h5path, dataset, read_vector)

            if (all(abs(read_vector - full_vector) < epsilon(full_vector))) then
               write(output_unit, '(a)') "    Vector read matches written values"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Vector read mismatch"
               num_failed = num_failed + 1
            end if

            deallocate(full_vector, read_vector)
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
   end subroutine test_write_read_real_real64_vector_serial

   ! Test for complex(real32) vector
   subroutine test_write_read_complex_real32_vector()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real32) vector (parallel mode)"
      character(*), parameter :: filename = "test_complex_real32_vector_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), parameter :: elements_per_rank = 10
      integer(int32), parameter :: total_elements = 20
      complex(real32), allocatable :: local_vector(:), read_vector(:)
      integer(int32) :: start_idx, local_size, i
#ifdef MPI
      integer(int32) :: datastart(1), datasize(1)
#endif

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call h5file%init(filename, serial_access=.false.)

      if (h5file%file_id /= 0) then
#ifdef MPI
         start_idx = my_rank * elements_per_rank + 1
         local_size = elements_per_rank

         allocate(local_vector(local_size))
         allocate(read_vector(local_size))
         do i = 1, local_size
            local_vector(i) = cmplx(real(start_idx + i - 1, real32), real(start_idx + i - 1, real32) + 0.5, real32)
         end do

         datastart(1) = start_idx
         datasize(1) = local_size * comm_size

         call h5file%write(h5path, dataset, local_vector, datastart=datastart, datasize=datasize)
         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " wrote elements ", start_idx, " to ", start_idx + local_size - 1

         call h5file%read(h5path, dataset, read_vector, datastart=datastart)

         if (all(abs(read_vector - local_vector) < epsilon(real(local_vector)))) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_vector, read_vector)
#else
         allocate(local_vector(total_elements))
         allocate(read_vector(total_elements))
         do i = 1, total_elements
            local_vector(i) = cmplx(real(i, real32), real(i, real32) + 0.5, real32)
         end do

         call h5file%write(h5path, dataset, local_vector)
         write(output_unit, '(a)') "    Full vector written successfully"

         call h5file%read(h5path, dataset, read_vector)

         if (all(abs(read_vector - local_vector) < epsilon(real(local_vector)))) then
            write(output_unit, '(a)') "    Vector read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Vector read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_vector, read_vector)
#endif
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_complex_real32_vector

   subroutine test_write_read_complex_real32_vector_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real32) vector (serial mode)"
      character(*), parameter :: filename = "test_complex_real32_vector_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), parameter :: total_elements = 20
      complex(real32), allocatable :: full_vector(:), read_vector(:)
      integer(int32) :: i

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call init_with_mpi_comm_world(h5file%mpi_comm)

      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         call h5file%init(filename, serial_access=.true.)

         if (h5file%file_id /= 0) then
            allocate(full_vector(total_elements))
            allocate(read_vector(total_elements))
            do i = 1, total_elements
               full_vector(i) = cmplx(real(i, real32), real(i, real32) + 0.5, real32)
            end do

            call h5file%write(h5path, dataset, full_vector)
            write(output_unit, '(a)') "    Full vector written successfully"

            call h5file%read(h5path, dataset, read_vector)

            if (all(abs(read_vector - full_vector) < epsilon(real(full_vector)))) then
               write(output_unit, '(a)') "    Vector read matches written values"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Vector read mismatch"
               num_failed = num_failed + 1
            end if

            deallocate(full_vector, read_vector)
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
   end subroutine test_write_read_complex_real32_vector_serial

   ! Test for complex(real64) vector
   subroutine test_write_read_complex_real64_vector()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real64) vector (parallel mode)"
      character(*), parameter :: filename = "test_complex_real64_vector_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), parameter :: elements_per_rank = 10
      integer(int32), parameter :: total_elements = 20
      complex(real64), allocatable :: local_vector(:), read_vector(:)
      integer(int32) :: start_idx, local_size, i
#ifdef MPI
      integer(int32) :: datastart(1), datasize(1)
#endif

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call h5file%init(filename, serial_access=.false.)

      if (h5file%file_id /= 0) then
#ifdef MPI
         start_idx = my_rank * elements_per_rank + 1
         local_size = elements_per_rank

         allocate(local_vector(local_size))
         allocate(read_vector(local_size))
         do i = 1, local_size
            local_vector(i) = cmplx(real(start_idx + i - 1, real64), real(start_idx + i - 1, real64) + 0.5, real64)
         end do

         datastart(1) = start_idx
         datasize(1) = local_size * comm_size

         call h5file%write(h5path, dataset, local_vector, datastart=datastart, datasize=datasize)
         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " wrote elements ", start_idx, " to ", start_idx + local_size - 1

         call h5file%read(h5path, dataset, read_vector, datastart=datastart)

         if (all(abs(read_vector - local_vector) < epsilon(real(local_vector)))) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_vector, read_vector)
#else
         allocate(local_vector(total_elements))
         allocate(read_vector(total_elements))
         do i = 1, total_elements
            local_vector(i) = cmplx(real(i, real64), real(i, real64) + 0.5, real64)
         end do

         call h5file%write(h5path, dataset, local_vector)
         write(output_unit, '(a)') "    Full vector written successfully"

         call h5file%read(h5path, dataset, read_vector)

         if (all(abs(read_vector - local_vector) < epsilon(real(local_vector)))) then
            write(output_unit, '(a)') "    Vector read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Vector read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_vector, read_vector)
#endif
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_complex_real64_vector

   subroutine test_write_read_complex_real64_vector_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real64) vector (serial mode)"
      character(*), parameter :: filename = "test_complex_real64_vector_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), parameter :: total_elements = 20
      complex(real64), allocatable :: full_vector(:), read_vector(:)
      integer(int32) :: i

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      call init_with_mpi_comm_world(h5file%mpi_comm)

      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         call h5file%init(filename, serial_access=.true.)

         if (h5file%file_id /= 0) then
            allocate(full_vector(total_elements))
            allocate(read_vector(total_elements))
            do i = 1, total_elements
               full_vector(i) = cmplx(real(i, real64), real(i, real64) + 0.5, real64)
            end do

            call h5file%write(h5path, dataset, full_vector)
            write(output_unit, '(a)') "    Full vector written successfully"

            call h5file%read(h5path, dataset, read_vector)

            if (all(abs(read_vector - full_vector) < epsilon(real(full_vector)))) then
               write(output_unit, '(a)') "    Vector read matches written values"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Vector read mismatch"
               num_failed = num_failed + 1
            end if

            deallocate(full_vector, read_vector)
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
   end subroutine test_write_read_complex_real64_vector_serial

   ! Test for uneven data distribution (different ranks write different amounts)
   subroutine test_write_read_integer_int32_vector_uneven()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read integer(int32) vector (uneven distribution)"
      character(*), parameter :: filename = "test_integer_int32_vector_uneven.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), allocatable :: local_vector(:), read_vector(:)
      integer(int32) :: local_size, start_idx, i
      type(hyperslab_type), allocatable :: hyperslabs(:)
      integer(int32) :: datastart(1), datasize(1)
      
#ifdef MPI

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%file_id /= 0) then
         ! Uneven distribution: rank 0 writes 7 elements, rank 1 writes 18 elements
         if (my_rank == 0) then
            local_size = 7
            start_idx = 1
         else if (my_rank == 1) then
            local_size = 18
            start_idx = 8
         else
            ! For more than 2 ranks, distribute remaining elements
            local_size = 25
            start_idx = 26 + (my_rank - 2) * local_size
         end if

         allocate(local_vector(local_size))
         allocate(read_vector(local_size))
         allocate(hyperslabs(1))

         ! Create local data
         do i = 1, local_size
            local_vector(i) = start_idx + i - 1
         end do

         ! Set up hyperslab for this rank's portion
         datastart(1) = start_idx
         if(comm_size == 1) then
            datasize(1) = local_size  ! Single rank writes all elements
         else
            datasize(1) = 25 + (comm_size - 2) * 25  ! Total size: 7 + 18 + 25*(comm_size-1)
         end if

         call hyperslabs(1)%init(shape(local_vector, kind=int32), datasize, [-1], datastart, [-1], [-1], [-1], .false.)

         write(output_unit, '(a,i0,a,i0,a,i0)') "    Rank ", my_rank, " writing ", local_size, " elements starting at ", start_idx

         ! Write using hyperslabs
         call h5file%write(h5path, dataset, local_vector, hyperslabs=hyperslabs)

         ! Read back using the same hyperslab
         call h5file%read(h5path, dataset, read_vector, hyperslabs=hyperslabs)

         ! Check if values match
         if (all(read_vector == local_vector)) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_vector, read_vector, hyperslabs)

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if
#else 
      ! Single process mode - do nothing 
      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name
      write(output_unit, '(a)') "    Single process mode - skipping uneven distribution test"
      num_passed = num_passed + 1  ! Count as passed since this is expected behavior
#endif


      write(output_unit, *)
   end subroutine test_write_read_integer_int32_vector_uneven

   ! Test for non-contiguous data distribution (multiple hyperslab selections per rank)
   subroutine test_write_read_integer_int32_vector_noncontiguous()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read integer(int32) vector (non-contiguous distribution)"
      character(*), parameter :: filename = "test_integer_int32_vector_noncontiguous.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_vector"
      integer(int32), parameter :: elements_per_block = 3  ! Elements per contiguous block
      integer(int32), parameter :: num_blocks = 2  ! Number of blocks per rank
      integer(int32), allocatable :: local_vector(:), read_vector(:)
      integer(int32) :: local_size, block_start, i, block_idx
      type(hyperslab_type), allocatable :: hyperslabs(:)
      integer(int32) :: datastart(1), datasize(1), memstart(1), memsize(1)
#ifdef MPI
      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      if (comm_size < 2) then
         write(output_unit, '(a)') "    Test requires at least 2 MPI ranks - skipping"
         num_passed = num_passed + 1
         write(output_unit, *)
         return
      end if

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%file_id /= 0) then
         ! Non-contiguous distribution: each rank writes to 2 separate blocks
         ! Rank 0: blocks at positions 1-3 and 7-9
         ! Rank 1: blocks at positions 4-6 and 10-12
         local_size = elements_per_block * num_blocks
         allocate(local_vector(local_size))
         allocate(read_vector(local_size))
         allocate(hyperslabs(num_blocks))

         ! Create local data and hyperslabs
         do block_idx = 1, num_blocks
            if (my_rank == 0) then
               block_start = (block_idx - 1) * 6 + 1  ! 1 and 7
            else if (my_rank == 1) then
               block_start = (block_idx - 1) * 6 + 4  ! 4 and 10
            else
               ! For additional ranks, distribute in pattern
               block_start = 13 + (my_rank - 2) * 6 + (block_idx - 1) * 3
            end if

            ! Fill data for this block
            do i = 1, elements_per_block
               local_vector((block_idx-1)*elements_per_block + i) = block_start + i - 1
            end do

            ! Set up hyperslab for this block
            datastart(1) = block_start
            datasize(1) = local_size * comm_size  ! Total dataset size in file
            memstart(1) = (block_idx - 1) * elements_per_block + 1
            memsize(1) = local_size

            call hyperslabs(block_idx)%init([elements_per_block], datasize, memsize, datastart, memstart, [-1], [-1], .false.)
         end do

         write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " writing non-contiguous blocks"

         ! Write using multiple hyperslabs
         call h5file%write(h5path, dataset, local_vector, hyperslabs=hyperslabs)

         ! Read back using the same hyperslabs
         call h5file%read(h5path, dataset, read_vector, hyperslabs=hyperslabs)

         ! Check if values match
         if (all(read_vector == local_vector)) then
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read matches written values"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,i0,a)') "    Rank ", my_rank, " read mismatch"
            num_failed = num_failed + 1
         end if

         deallocate(local_vector, read_vector, hyperslabs)
         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

#else 
      ! Single process mode - do nothing 
      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name
      write(output_unit, '(a)') "    Single process mode - skipping uneven distribution test"
      num_passed = num_passed + 1  ! Count as passed since this is expected behavior
#endif
      write(output_unit, *)
   end subroutine test_write_read_integer_int32_vector_noncontiguous

end program write_read_numerical_vector_test
