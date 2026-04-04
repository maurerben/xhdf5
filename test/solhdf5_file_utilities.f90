!> Unit tests for solhdf5 init and delete routines
program solhdf5_file_utilities
#ifdef MPI
   use mpi_f08
#endif
   use iso_fortran_env, only: output_unit
   use solhdf5
   use os_utils, only: join_paths

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
   write(output_unit, '(a)') "  solhdf5 Init/Delete Unit Tests"
   write(output_unit, '(a)') "=========================================="
   write(output_unit, *)

#ifdef _HDF5_
   ! Test 1: Init and delete with serial mode
   call test_init_delete_serial()

   ! Test 2: Init creates new file
   call test_init_creates_file()

   ! Test 3: Init opens existing file
   call test_init_opens_existing()

   ! Test 4: Multiple init/delete cycles
   call test_multiple_cycles()

#ifdef MPI
   ! Tests 5-6 only make sense with real MPI communicator
   ! Test 5: Init with non-serial mode (parallel mode with MPI_COMM_WORLD)
   call test_init_nonserial_mode()

   ! Test 6: Non-serial mode with existing file
   call test_init_nonserial_opens_existing()

   ! Test 7: initialize_group in parallel mode
   call test_initialize_group()

   ! Test 8: initialize_group_update_groupname in parallel mode
   call test_initialize_group_update_groupname()

   ! Test 9: link_exists in parallel mode
   call test_link_exists()

   ! Test 10: delete_link in parallel mode
   call test_delete_link()
#endif

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

   subroutine test_init_delete_serial()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 1: Init and delete (serial mode)"
      character(*), parameter :: filename = "test_init_delete_serial.h5"

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in serial mode
      call h5file%init(filename, serial_access=.true.)

      ! Check file was initialized
      if (h5file%file_id /= 0) then
         write(output_unit, '(a)') "    Init succeeded"

         ! Delete file
         call h5file%delete()

         ! Check file was deleted
         if (h5file%file_id == 0) then
            write(output_unit, '(a)') "    Delete succeeded"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    File ID not reset after delete"
            num_failed = num_failed + 1
         end if
      else
         write(output_unit, '(a)') "    Init failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_init_delete_serial


   subroutine test_init_creates_file()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 2: Init creates new file"
      character(*), parameter :: filename = "test_create_new.h5"
      logical :: file_exists_after

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init - should create new file
      call h5file%init(filename, serial_access=.true.)

      if (h5file%file_id /= 0) then
         write(output_unit, '(a)') "    File opened/created"

         ! Check root exists (file is valid)
         file_exists_after = h5file%link_exists("/")

         if (file_exists_after) then
            write(output_unit, '(a)') "    File is valid (root exists)"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    File root does not exist"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    Failed to create file"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_init_creates_file


   subroutine test_init_opens_existing()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 3: Init opens existing file"
      character(*), parameter :: filename = "test_open_existing.h5"

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! First, create a file
      call h5file%init(filename, serial_access=.true.)

      if (h5file%file_id /= 0) then
         call h5file%delete()

         ! Now open the same file
         call h5file%init(filename, serial_access=.true.)

         if (h5file%file_id /= 0) then
            write(output_unit, '(a)') "    Successfully opened existing file"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Failed to open existing file"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    Failed to create initial file"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_init_opens_existing


   subroutine test_multiple_cycles()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 4: Multiple init/delete cycles"
      character(*), parameter :: filename = "test_multiple_cycles.h5"
      integer :: i, cycle_count

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      cycle_count = 3
      do i = 1, cycle_count
         ! Init
         call h5file%init(filename, serial_access=.true.)

         if (h5file%file_id == 0) then
            write(output_unit, '(a,i0,a)') "    Init failed in cycle ", i, ""
            num_failed = num_failed + 1
            return
         end if

         ! Delete
         call h5file%delete()

         if (h5file%file_id /= 0) then
            write(output_unit, '(a,i0,a)') "    Delete failed in cycle ", i, ""
            num_failed = num_failed + 1
            return
         end if
      end do

      write(output_unit, '(a,i0,a)') "    Successfully completed ", cycle_count, " init/delete cycles"
      num_passed = num_passed + 1
      write(output_unit, *)
   end subroutine test_multiple_cycles


   subroutine test_init_nonserial_mode()
#ifdef MPI
      use mpi_f08, only: MPI_COMM_WORLD
#endif
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 5: Init with non-serial mode (MPI_COMM_WORLD)"
      character(*), parameter :: filename = "test_nonserial_mode.h5"

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

#ifdef MPI
      ! Init file in non-serial mode with MPI_COMM_WORLD
      call h5file%init(filename, MPI_COMM_WORLD, serial_access=.false.)

      if (h5file%file_id /= 0) then
         write(output_unit, '(a)') "    Init succeeded in non-serial mode"

         ! Verify serial_access flag is false
         if (.not. h5file%serial_access) then
            write(output_unit, '(a)') "    Serial_access flag correctly set to false"
            num_passed = num_passed + 1
         else
            write(output_unit, '(a)') "    Serial_access flag incorrectly set"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    Init failed in non-serial mode"
         num_failed = num_failed + 1
      end if
#else
      write(output_unit, '(a)') "    Skipped (MPI not available)"
      num_passed = num_passed + 1
#endif

      write(output_unit, *)
   end subroutine test_init_nonserial_mode


   subroutine test_init_nonserial_opens_existing()
#ifdef MPI
      use mpi_f08, only: MPI_COMM_WORLD
#endif
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 6: Non-serial mode with existing file (MPI_COMM_WORLD)"
      character(*), parameter :: filename = "test_nonserial_existing.h5"

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

#ifdef MPI
      ! First, create file in serial mode
      call h5file%init(filename, MPI_COMM_WORLD, serial_access=.true.)

      if (h5file%file_id /= 0) then
         call h5file%delete()

         ! Open same file in non-serial (parallel) mode with MPI_COMM_WORLD
         call h5file%init(filename, MPI_COMM_WORLD, serial_access=.false.)

         if (h5file%file_id /= 0) then
            write(output_unit, '(a)') "    Successfully opened existing file in non-serial mode"

            if (.not. h5file%serial_access) then
               write(output_unit, '(a)') "    Serial_access flag correctly set to false"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Serial_access flag incorrectly set"
               num_failed = num_failed + 1
            end if
         else
            write(output_unit, '(a)') "    Failed to open existing file in non-serial mode"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    Failed to create initial file"
         num_failed = num_failed + 1
      end if
#else
      write(output_unit, '(a)') "    Skipped (MPI not available)"
      num_passed = num_passed + 1
#endif

      write(output_unit, *)
   end subroutine test_init_nonserial_opens_existing


   subroutine test_initialize_group()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 7: initialize_group in parallel mode"
      character(*), parameter :: filename = "test_initialize_group.h5"
      character(*), parameter :: group_path = "/"
      character(*), parameter :: group_name = "test_group"

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%file_id /= 0) then
         ! Check group doesn't exist initially
         if (.not. h5file%link_exists(join_paths(group_path, group_name))) then
            write(output_unit, '(a)') "    Group does not exist initially"

            ! Create group
            call h5file%initialize_group(group_path, group_name)

            ! Check group exists now
            if (h5file%link_exists(join_paths(group_path, group_name))) then
               write(output_unit, '(a)') "    Group created successfully"

               ! Try to create again (should not fail)
               call h5file%initialize_group(group_path, group_name)

               if (h5file%link_exists(join_paths(group_path, group_name))) then
                  write(output_unit, '(a)') "    Group still exists after second create"
                  num_passed = num_passed + 1
               else
                  write(output_unit, '(a)') "    Group disappeared after second create"
                  num_failed = num_failed + 1
               end if
            else
               write(output_unit, '(a)') "    Group creation failed"
               num_failed = num_failed + 1
            end if
         else
            write(output_unit, '(a)') "    Group already exists unexpectedly"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File init failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_initialize_group


   subroutine test_initialize_group_update_groupname()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 8: initialize_group_update_groupname in parallel mode"
      character(*), parameter :: filename = "test_initialize_group_update.h5"
      character(*), parameter :: group_path = "/"
      character(:), allocatable :: group_name

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, MPI_COMM_WORLD, serial_access=.false.)

      if (h5file%file_id /= 0) then
         group_name = "test_group_update"

         ! Create group and update groupname
         call h5file%initialize_group_update_groupname(group_path, group_name)

         ! Check group exists
         if (h5file%link_exists(group_name)) then
            write(output_unit, '(a)') "    Group created successfully"

            ! Check groupname was updated to full path
            if (group_name == join_paths(group_path, "test_group_update")) then
               write(output_unit, '(a)') "    Group name updated to full path correctly"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Group name not updated correctly"
               num_failed = num_failed + 1
            end if
         else
            write(output_unit, '(a)') "    Group creation failed"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File init failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_initialize_group_update_groupname


   subroutine test_link_exists()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 9: link_exists in parallel mode"
      character(*), parameter :: filename = "test_link_exists.h5"
      character(*), parameter :: group_path = "/"
      character(*), parameter :: group_name = "test_link_group"

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, serial_access=.false.)

      if (h5file%file_id /= 0) then
         ! Check root exists
         if (h5file%link_exists("/")) then
            write(output_unit, '(a)') "    Root link exists"

            ! Check non-existing group doesn't exist
            if (.not. h5file%link_exists(join_paths(group_path, group_name))) then
               write(output_unit, '(a)') "    Non-existing group correctly reported as not existing"

               ! Create group
               call h5file%initialize_group(group_path, group_name)

               ! Check group exists now
               if (h5file%link_exists(join_paths(group_path, group_name))) then
                  write(output_unit, '(a)') "    Created group exists"
                  num_passed = num_passed + 1
               else
                  write(output_unit, '(a)') "    Created group does not exist"
                  num_failed = num_failed + 1
               end if
            else
               write(output_unit, '(a)') "    Non-existing group incorrectly reported as existing"
               num_failed = num_failed + 1
            end if
         else
            write(output_unit, '(a)') "    Root link does not exist"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File init failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_link_exists


   subroutine test_delete_link()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 10: delete_link in parallel mode"
      character(*), parameter :: filename = "test_delete_link.h5"
      character(*), parameter :: group_path = "/"
      character(*), parameter :: group_name = "test_delete_group"

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Init file in parallel mode
      call h5file%init(filename, MPI_COMM_WORLD, serial_access=.false.)

      if (h5file%file_id /= 0) then
         ! Create group
         call h5file%initialize_group(group_path, group_name)

         if (h5file%link_exists(join_paths(group_path, group_name))) then
            write(output_unit, '(a)') "    Group created successfully"

            ! Delete group
            call h5file%delete_link(join_paths(group_path, group_name))

            ! Check group no longer exists
            if (.not. h5file%link_exists(join_paths(group_path, group_name))) then
               write(output_unit, '(a)') "    Group deleted successfully"
               num_passed = num_passed + 1
            else
               write(output_unit, '(a)') "    Group still exists after delete"
               num_failed = num_failed + 1
            end if
         else
            write(output_unit, '(a)') "    Group creation failed"
            num_failed = num_failed + 1
         end if

         call h5file%delete()
      else
         write(output_unit, '(a)') "    File init failed"
         num_failed = num_failed + 1
      end if

      write(output_unit, *)
   end subroutine test_delete_link

end program solhdf5_file_utilities
