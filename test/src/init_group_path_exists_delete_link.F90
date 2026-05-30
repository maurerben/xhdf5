! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Benedikt Maurer

program init_group_path_exists_delete_link
#ifdef MPI
   use mpi_f08
#endif
   use iso_fortran_env, only: output_unit
   use xhdf5
   use os_utils, only: join_paths

   implicit none

   integer :: num_tests = 0
   integer :: num_passed = 0
   integer :: num_failed = 0
#ifdef MPI
   integer :: mpi_err
   logical :: mpi_initialized_flag, mpi_init_by_me
#endif

#ifdef MPI
   ! Check if MPI is already initialized (e.g., when running with mpiexec)
   call mpi_initialized(mpi_initialized_flag, mpi_err)
   if (.not. mpi_initialized_flag) then
      call mpi_init(mpi_err)
      mpi_init_by_me = .true.
   else
      mpi_init_by_me = .false.
   end if
#endif

   write(output_unit, '(a)') "=========================================="
   write(output_unit, '(a)') "  xhdf5 Init/Delete Unit Tests"
   write(output_unit, '(a)') "=========================================="
   write(output_unit, *)

#ifdef _HDF5_

   ! Test 1: init_group
   call test_init_group()

   ! Test 2: init_group_update_groupname
   call test_init_group_update_groupname()

   ! Test 3: link_exists
   call test_link_exists()

   ! Test 4: delete_link
   call test_delete_link()
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
      if (mpi_init_by_me) call mpi_finalize(mpi_err)
#endif
      stop 1
   end if

#ifdef MPI
   ! Finalize MPI
   if (mpi_init_by_me) call mpi_finalize(mpi_err)
#endif

contains


   subroutine test_init_group()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 1: init_group "
      character(*), parameter :: filename = "test_init_group.h5"
      character(*), parameter :: group_path = "/"
      character(*), parameter :: group_name = "test_group"

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Clean up any existing file
      call execute_command_line("rm -f " // filename)

      ! Init file
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Check group doesn't exist initially
         if (.not. h5file%link_exists(join_paths(group_path, group_name))) then
            write(output_unit, '(a)') "    Group does not exist initially"

            ! Create group
            call h5file%init_group(group_path, group_name)

            ! Check group exists now
            if (h5file%link_exists(join_paths(group_path, group_name))) then
               write(output_unit, '(a)') "    Group created successfully"

               ! Try to create again (should not fail)
               call h5file%init_group(group_path, group_name)

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
   end subroutine test_init_group


   subroutine test_init_group_update_groupname()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 2: init_group_update_groupname "
      character(*), parameter :: filename = "test_init_group_update.h5"
      character(*), parameter :: group_path = "/"
      character(:), allocatable :: group_name

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Clean up any existing file
      call execute_command_line("rm -f " // filename)

      ! Init file
#ifdef MPI
      call h5file%init(filename, MPI_COMM_WORLD, serial_access=.false.)
#else
      call h5file%init(filename, serial_access=.false.)
#endif

      if (h5file%is_open()) then
         group_name = "test_group_update"

         ! Create group and update groupname
         call h5file%init_group_update_groupname(group_path, group_name)

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
   end subroutine test_init_group_update_groupname


   subroutine test_link_exists()
      type(h5file_t) :: h5file
      character(*), parameter :: test_name = "Test 3: link_exists "
      character(*), parameter :: filename = "test_link_exists.h5"
      character(*), parameter :: group_path = "/"
      character(*), parameter :: group_name = "test_link_group"

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Clean up any existing file
      call execute_command_line("rm -f " // filename)

      ! Init file
      call h5file%init(filename, serial_access=.false.)

      if (h5file%is_open()) then
         ! Check root exists
         if (h5file%link_exists("/")) then
            write(output_unit, '(a)') "    Root link exists"

            ! Check non-existing group doesn't exist
            if (.not. h5file%link_exists(join_paths(group_path, group_name))) then
               write(output_unit, '(a)') "    Non-existing group correctly reported as not existing"

               ! Create group
               call h5file%init_group(group_path, group_name)

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
      character(*), parameter :: test_name = "Test 4: delete_link "
      character(*), parameter :: filename = "test_delete_link.h5"
      character(*), parameter :: group_path = "/"
      character(*), parameter :: group_name = "test_delete_group"

      num_tests = num_tests + 1
      write(output_unit, '(a)') test_name

      ! Clean up any existing file
      call execute_command_line("rm -f " // filename)

      ! Init file
#ifdef MPI
      call h5file%init(filename, MPI_COMM_WORLD, serial_access=.false.)
#else
      call h5file%init(filename, serial_access=.false.)
#endif

      if (h5file%is_open()) then
         ! Create group
         call h5file%init_group(group_path, group_name)

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


end program init_group_path_exists_delete_link
