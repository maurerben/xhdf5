!> Unit tests for write and read string routines
program write_read_string_test
#ifdef MPI
   use mpi_f08
#endif
   use iso_fortran_env, only: output_unit
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
   write(output_unit, '(a)') "  solhdf5 Write/Read Unit Tests"
   write(output_unit, '(a)') "=========================================="
   write(output_unit, *)

#ifdef _HDF5_
   ! Test 1: Write and read string
   call test_write_read_string()

   ! Test 2: Write and read string in serial mode
   call test_write_read_string_serial()

   ! Test 3: Write and read empty string
   call test_write_read_empty_string()

   ! Test 4: Write and read long string
   call test_write_read_long_string()

   ! Test 5: Overwrite string dataset
   call test_overwrite_string_dataset()

   ! Test 6: Write and read boolean
   call test_write_read_bool()

   ! Test 7: Write and read boolean in serial mode
   call test_write_read_bool_serial()

   ! Test 8: Overwrite boolean dataset
   call test_overwrite_bool_dataset()
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

   subroutine test_write_read_string()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 1: Write and read string (parallel mode)"
      character(*), parameter :: filename = "test_write_read_string_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_string"
      character(*), parameter :: test_string = "Hello, HDF5 World!"
      character(:), allocatable :: read_string

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Write string
         call h5file%write(h5path, dataset, test_string)
         write(output_unit, '(a)') "    String written successfully"

         ! Read string back
         call h5file%read(h5path, dataset, read_string)

         ! Check if strings match
         if (allocated(read_string)) then
            if (trim(read_string) == trim(test_string)) then
               write(output_unit, '(a)') "    String read matches written string"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a,a,a)') "    String mismatch: expected '", trim(test_string), &
                  "', got '", trim(read_string), "'"
               num_failed = num_failed + 1
            end if
         else
            write(output_unit, '(a)') "    Read string not allocated"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_string


   subroutine test_write_read_string_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 2: Write and read string (serial mode)"
      character(*), parameter :: filename = "test_write_read_string_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_string"
      character(*), parameter :: test_string = "Hello, HDF5 World!"
      character(:), allocatable :: read_string

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Initialize MPI communicator first to check rank
      call init_with_mpi_comm_world(h5file%mpi_comm)

      ! Only root process performs serial HDF5 operations
      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         ! Init file in serial mode (only root process allowed)
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            ! Write string
            call h5file%write(h5path, dataset, test_string)
            write(output_unit, '(a)') "    String written successfully"

            ! Read string back
            call h5file%read(h5path, dataset, read_string)

            ! Check if strings match
            if (allocated(read_string)) then
               if (trim(read_string) == trim(test_string)) then
                  write(output_unit, '(a)') "    String read matches written string"
                  num_passed = num_passed + 1
               else
                  write(output_unit, '(a,a,a)') "    String mismatch: expected '", trim(test_string), &
                     "', got '", trim(read_string), "'"
                  num_failed = num_failed + 1
               end if
            else
               write(output_unit, '(a)') "    Read string not allocated"
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
   end subroutine test_write_read_string_serial


   subroutine test_write_read_empty_string()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 3: Write and read empty string"
      character(*), parameter :: filename = "test_write_read_empty_string.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "empty_string"
      character(*), parameter :: test_string = ""
      character(:), allocatable :: read_string

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Write empty string
         call h5file%write(h5path, dataset, test_string)
         write(output_unit, '(a)') "    Empty string written successfully"

         ! Read string back
         call h5file%read(h5path, dataset, read_string)

         ! Check if strings match
         if (allocated(read_string)) then
            if (len_trim(read_string) == 0) then
               write(output_unit, '(a)') "    Empty string read correctly"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a,a)') "    Expected empty string, got '", trim(read_string), "'"
               num_failed = num_failed + 1
            end if
         else
            write(output_unit, '(a)') "    Read string not allocated"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_empty_string


   subroutine test_write_read_long_string()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 4: Write and read long string"
      character(*), parameter :: filename = "test_write_read_long_string.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "long_string"
      character(500) :: test_string
      character(:), allocatable :: read_string
      integer :: i

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Create a long test string
      test_string = ""
      do i = 1, 50
         test_string = trim(test_string) // "This is a test string. "
      end do

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Write long string
         call h5file%write(h5path, dataset, trim(test_string))
         write(output_unit, '(a)') "    Long string written successfully"

         ! Read string back
         call h5file%read(h5path, dataset, read_string)

         ! Check if strings match
         if (allocated(read_string)) then
            if (trim(read_string) == trim(test_string)) then
               write(output_unit, '(a)') "    Long string read matches written string"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a,i0,a,i0)') "    String length mismatch: expected ", len_trim(test_string), &
                  ", got ", len_trim(read_string)
               num_failed = num_failed + 1
            end if
         else
            write(output_unit, '(a)') "    Read string not allocated"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_long_string


   subroutine test_overwrite_string_dataset()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 5: Overwrite string dataset"
      character(*), parameter :: filename = "test_overwrite_string_dataset.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "overwrite_string"
      character(*), parameter :: original_string = "Original string"
      character(*), parameter :: new_string = "New overwritten string"
      character(:), allocatable :: read_string

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Write original string
         call h5file%write(h5path, dataset, original_string)
         write(output_unit, '(a)') "    Original string written"

         ! Overwrite with new string
         call h5file%write(h5path, dataset, new_string)
         write(output_unit, '(a)') "    String overwritten"

         ! Read string back
         call h5file%read(h5path, dataset, read_string)

         ! Check if the new string is read
         if (allocated(read_string)) then
            if (trim(read_string) == trim(new_string)) then
               write(output_unit, '(a)') "    Overwritten string read correctly"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a,a,a)') "    Expected '", trim(new_string), &
                  "', got '", trim(read_string), "'"
               num_failed = num_failed + 1
            end if
         else
            write(output_unit, '(a)') "    Read string not allocated"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_overwrite_string_dataset


   subroutine test_write_read_bool()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 6: Write and read boolean"
      character(*), parameter :: filename = "test_write_read_bool_parallel.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_bool"
      logical :: test_bool = .true.
      logical :: read_bool

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Write boolean
         call h5file%write(h5path, dataset, test_bool)
         write(output_unit, '(a)') "    Boolean written successfully"

         ! Read boolean back
         call h5file%read(h5path, dataset, read_bool)

         ! Check if booleans match
         if (read_bool .eqv. test_bool) then
            write(output_unit, '(a)') "    Boolean read matches written boolean"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,l1,a,l1)') "    Boolean mismatch: expected ", test_bool, &
               ", got ", read_bool
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_write_read_bool


   subroutine test_write_read_bool_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 7: Write and read boolean (serial mode)"
      character(*), parameter :: filename = "test_write_read_bool_serial.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "test_bool"
      logical :: test_bool = .false.
      logical :: read_bool

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Initialize MPI communicator first to check rank
      call init_with_mpi_comm_world(h5file%mpi_comm)

      ! Only root process performs serial HDF5 operations
      if (comm_to_rank(h5file%mpi_comm) == root_rank) then
         ! Init file in serial mode (only root process allowed)
         call h5file%init(filename, serial_access=.true.)

         if (h5file%is_open()) then
            ! Write boolean
            call h5file%write(h5path, dataset, test_bool)
            write(output_unit, '(a)') "    Boolean written successfully"

            ! Read boolean back
            call h5file%read(h5path, dataset, read_bool)

            ! Check if booleans match
            if (read_bool .eqv. test_bool) then
               write(output_unit, '(a)') "    Boolean read matches written boolean"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a,l1,a,l1)') "    Boolean mismatch: expected ", test_bool, &
                  ", got ", read_bool
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
   end subroutine test_write_read_bool_serial


   subroutine test_overwrite_bool_dataset()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 8: Overwrite boolean dataset"
      character(*), parameter :: filename = "test_overwrite_bool_dataset.h5"
      character(*), parameter :: h5path = "/"
      character(*), parameter :: dataset = "overwrite_bool"
      logical :: original_bool = .true.
      logical :: new_bool = .false.
      logical :: read_bool

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Write original boolean
         call h5file%write(h5path, dataset, original_bool)
         write(output_unit, '(a)') "    Original boolean written"

         ! Overwrite with new boolean
         call h5file%write(h5path, dataset, new_bool)
         write(output_unit, '(a)') "    Boolean overwritten"

         ! Read boolean back
         call h5file%read(h5path, dataset, read_bool)

         ! Check if the new boolean is read
         if (read_bool .eqv. new_bool) then
            write(output_unit, '(a)') "    Overwritten boolean read correctly"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a,l1,a,l1)') "    Expected ", new_bool, &
               ", got ", read_bool
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File initialization failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_overwrite_bool_dataset

end program write_read_string_test
