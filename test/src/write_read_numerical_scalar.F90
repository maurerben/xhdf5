!> Unit tests for write and read numerical scalar routines
program write_read_numerical_scalar_test
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
   integer :: mpi_err
#endif

#ifdef MPI
   ! Initialize MPI
   call mpi_init(mpi_err)
#endif

   write(output_unit, '(a)') "=========================================="
   write(output_unit, '(a)') "  solhdf5 Numerical Scalar Unit Tests"
   write(output_unit, '(a)') "=========================================="
   write(output_unit, *)

#ifdef _HDF5_
   ! Test integer(int32) scalars
   call test_write_read_integer_int32_scalar()
   call test_write_read_integer_int32_scalar_serial()

   ! Test real(real32) scalars
   call test_write_read_real_real32_scalar()
   call test_write_read_real_real32_scalar_serial()

   ! Test real(real64) scalars
   call test_write_read_real_real64_scalar()
   call test_write_read_real_real64_scalar_serial()

   ! Test complex(real32) scalars
   call test_write_read_complex_real32_scalar()
   call test_write_read_complex_real32_scalar_serial()

   ! Test complex(real64) scalars
   call test_write_read_complex_real64_scalar()
   call test_write_read_complex_real64_scalar_serial()
#else
   write(output_unit, '(a)') "HDF5 not available - tests skipped"
#endif

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

   ! Test for integer(int32) scalar
   subroutine test_write_read_integer_int32_scalar()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read integer(int32) scalar (parallel mode)"
      character(*), parameter :: filename = "test_integer_int32_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_scalar"
      integer(int32) :: test_value = 42_int32
      integer(int32) :: read_value

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Write scalar
         call h5file%write(h5path, dataset, test_value)
         write(output_unit, '(a)') "    Scalar written successfully"

         ! Read scalar back
         call h5file%read(h5path, dataset, read_value)

         ! Check if values match
         if (read_value == test_value) then
            write(output_unit, '(a)') "    Scalar read matches written value"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Scalar mismatch: expected ", test_value, ", got ", read_value
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_integer_int32_scalar

   subroutine test_write_read_integer_int32_scalar_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read integer(int32) scalar (serial mode)"
      character(*), parameter :: filename = "test_integer_int32_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_scalar"
      integer(int32) :: test_value = 42_int32
      integer(int32) :: read_value

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Initialize MPI communicator first to check rank
      call init_with_mpi_comm_world(h5file%mpi_comm)

      ! Only root process performs serial HDF5 operations
      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         ! Init file in serial mode (only root process allowed)
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            ! Write scalar
            call h5file%write(h5path, dataset, test_value)
            write(output_unit, '(a)') "    Scalar written successfully"

            ! Read scalar back
            call h5file%read(h5path, dataset, read_value)

            ! Check if values match
            if (read_value == test_value) then
               write(output_unit, '(a)') "    Scalar read matches written value"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Scalar mismatch: expected ", test_value, ", got ", read_value
               num_failed = num_failed + 1
            end if

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
   end subroutine test_write_read_integer_int32_scalar_serial

   ! Test for real(real32) scalar
   subroutine test_write_read_real_real32_scalar()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real32) scalar (parallel mode)"
      character(*), parameter :: filename = "test_real_real32_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_scalar"
      real(real32) :: test_value = 3.14159_real32
      real(real32) :: read_value

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Write scalar
         call h5file%write(h5path, dataset, test_value)
         write(output_unit, '(a)') "    Scalar written successfully"

         ! Read scalar back
         call h5file%read(h5path, dataset, read_value)

         ! Check if values match
         if (read_value == test_value) then
            write(output_unit, '(a)') "    Scalar read matches written value"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Scalar mismatch: expected ", test_value, ", got ", read_value
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_real_real32_scalar

   subroutine test_write_read_real_real32_scalar_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real32) scalar (serial mode)"
      character(*), parameter :: filename = "test_real_real32_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_scalar"
      real(real32) :: test_value = 3.14159_real32
      real(real32) :: read_value

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Initialize MPI communicator first to check rank
      call init_with_mpi_comm_world(h5file%mpi_comm)

      ! Only root process performs serial HDF5 operations
      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         ! Init file in serial mode (only root process allowed)
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            ! Write scalar
            call h5file%write(h5path, dataset, test_value)
            write(output_unit, '(a)') "    Scalar written successfully"

            ! Read scalar back
            call h5file%read(h5path, dataset, read_value)

            ! Check if values match
            if (read_value == test_value) then
               write(output_unit, '(a)') "    Scalar read matches written value"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Scalar mismatch: expected ", test_value, ", got ", read_value
               num_failed = num_failed + 1
            end if

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
   end subroutine test_write_read_real_real32_scalar_serial

   ! Test for real(real64) scalar
   subroutine test_write_read_real_real64_scalar()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real64) scalar (parallel mode)"
      character(*), parameter :: filename = "test_real_real64_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_scalar"
      real(real64) :: test_value = 2.718281828459045_real64
      real(real64) :: read_value

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Write scalar
         call h5file%write(h5path, dataset, test_value)
         write(output_unit, '(a)') "    Scalar written successfully"

         ! Read scalar back
         call h5file%read(h5path, dataset, read_value)

         ! Check if values match
         if (read_value == test_value) then
            write(output_unit, '(a)') "    Scalar read matches written value"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Scalar mismatch: expected ", test_value, ", got ", read_value
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_real_real64_scalar

   subroutine test_write_read_real_real64_scalar_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read real(real64) scalar (serial mode)"
      character(*), parameter :: filename = "test_real_real64_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_scalar"
      real(real64) :: test_value = 2.718281828459045_real64
      real(real64) :: read_value

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Initialize MPI communicator first to check rank
      call init_with_mpi_comm_world(h5file%mpi_comm)

      ! Only root process performs serial HDF5 operations
      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         ! Init file in serial mode (only root process allowed)
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            ! Write scalar
            call h5file%write(h5path, dataset, test_value)
            write(output_unit, '(a)') "    Scalar written successfully"

            ! Read scalar back
            call h5file%read(h5path, dataset, read_value)

            ! Check if values match
            if (read_value == test_value) then
               write(output_unit, '(a)') "    Scalar read matches written value"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Scalar mismatch: expected ", test_value, ", got ", read_value
               num_failed = num_failed + 1
            end if

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
   end subroutine test_write_read_real_real64_scalar_serial

   ! Test for complex(real32) scalar
   subroutine test_write_read_complex_real32_scalar()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real32) scalar (parallel mode)"
      character(*), parameter :: filename = "test_complex_real32_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_scalar"
      complex(real32) :: test_value = (1.0_real32, 2.0_real32)
      complex(real32) :: read_value

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Write scalar
         call h5file%write(h5path, dataset, test_value)
         write(output_unit, '(a)') "    Scalar written successfully"

         ! Read scalar back
         call h5file%read(h5path, dataset, read_value)

         ! Check if values match
         if (read_value == test_value) then
            write(output_unit, '(a)') "    Scalar read matches written value"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Scalar mismatch: expected ", test_value, ", got ", read_value
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_complex_real32_scalar

   subroutine test_write_read_complex_real32_scalar_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real32) scalar (serial mode)"
      character(*), parameter :: filename = "test_complex_real32_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_scalar"
      complex(real32) :: test_value = (1.0_real32, 2.0_real32)
      complex(real32) :: read_value

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Initialize MPI communicator first to check rank
      call init_with_mpi_comm_world(h5file%mpi_comm)

      ! Only root process performs serial HDF5 operations
      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         ! Init file in serial mode (only root process allowed)
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            ! Write scalar
            call h5file%write(h5path, dataset, test_value)
            write(output_unit, '(a)') "    Scalar written successfully"

            ! Read scalar back
            call h5file%read(h5path, dataset, read_value)

            ! Check if values match
            if (read_value == test_value) then
               write(output_unit, '(a)') "    Scalar read matches written value"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Scalar mismatch: expected ", test_value, ", got ", read_value
               num_failed = num_failed + 1
            end if

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
   end subroutine test_write_read_complex_real32_scalar_serial

   ! Test for complex(real64) scalar
   subroutine test_write_read_complex_real64_scalar()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real64) scalar (parallel mode)"
      character(*), parameter :: filename = "test_complex_real64_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_scalar"
      complex(real64) :: test_value = (3.14159_real64, 2.71828_real64)
      complex(real64) :: read_value

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Write scalar
         call h5file%write(h5path, dataset, test_value)
         write(output_unit, '(a)') "    Scalar written successfully"

         ! Read scalar back
         call h5file%read(h5path, dataset, read_value)

         ! Check if values match
         if (read_value == test_value) then
            write(output_unit, '(a)') "    Scalar read matches written value"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Scalar mismatch: expected ", test_value, ", got ", read_value
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_complex_real64_scalar

   subroutine test_write_read_complex_real64_scalar_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test: Write and read complex(real64) scalar (serial mode)"
      character(*), parameter :: filename = "test_complex_real64_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_scalar"
      complex(real64) :: test_value = (3.14159_real64, 2.71828_real64)
      complex(real64) :: read_value

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Initialize MPI communicator first to check rank
      call init_with_mpi_comm_world(h5file%mpi_comm)

      ! Only root process performs serial HDF5 operations
      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         ! Init file in serial mode (only root process allowed)
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            ! Write scalar
            call h5file%write(h5path, dataset, test_value)
            write(output_unit, '(a)') "    Scalar written successfully"

            ! Read scalar back
            call h5file%read(h5path, dataset, read_value)

            ! Check if values match
            if (read_value == test_value) then
               write(output_unit, '(a)') "    Scalar read matches written value"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Scalar mismatch: expected ", test_value, ", got ", read_value
               num_failed = num_failed + 1
            end if

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
   end subroutine test_write_read_complex_real64_scalar_serial

end program write_read_numerical_scalar_test
