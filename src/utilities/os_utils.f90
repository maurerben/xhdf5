! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Benedikt Maurer

module os_utils

   implicit none

   private
   public :: path_exists, join_paths, separate_path_and_filename

contains

   !> check if path exists
   function path_exists(path, ierr) result(exists)
      !> path name
      character(*), intent(in) :: path
      !> system dependent error code; `0` on success
      integer, intent(out), optional :: ierr
      logical :: exists

      integer :: ierr_local
      character(:), allocatable :: dir

#ifdef __INTEL_COMPILER
      dir = trim(path)
      if( dir(len(dir):len(dir)) == '/' ) then
         inquire(directory=trim(path), exist=exists, iostat=ierr_local)
      else
         inquire(file=trim(path), exist=exists, iostat=ierr_local)
      end if
#else
      inquire(file=trim(path), exist=exists, iostat=ierr_local)
#endif
      if (present(ierr)) ierr = ierr_local
   end function path_exists

   !> Join two paths following UNIX convention:
   !>
   !> - `join_paths('path/one/', 'path/two') --> 'path/one/path/two'
   !>
   !> - `join_paths('path/one', 'path/two') --> 'path/one/path/two'
   !>
   !> - `join_paths('path/one', '/path/two') --> 'path/one/path/two'
   !>
   !> - `join_paths('path/one/', '/path/two') --> 'path/one/path/two'
   !>
   !> - `join_paths('path/one////', '///path/two') --> 'path/one/path/two'
   function join_paths(path1, path2) result(joined_path)
      !> Paths to join.
      character(*), intent(in) :: path1, path2
      character(:), allocatable :: joined_path

      character(:), allocatable :: path1_cleaned, path2_cleaned

      character(1), parameter :: delimiter = '/'

      integer :: last_char_index

      ! Remove delimiter from the end of path1
      path1_cleaned = adjustl(trim(path1))
      last_char_index = len(path1_cleaned)
      if (path1_cleaned == delimiter) then
         path1_cleaned = ""
      else if (path1_cleaned(last_char_index:last_char_index) == delimiter) then
         path1_cleaned = path1_cleaned(:last_char_index - 1)
      end if

      ! Remove delimiter from the beginning of path2
      path2_cleaned = adjustl(trim(path2))
      if (path2_cleaned == delimiter) then
         path2_cleaned = ""
      else if (path2_cleaned(1:1) == delimiter) then
         path2_cleaned = path2_cleaned(2:)
      end if

      ! Join paths
      joined_path = path1_cleaned // delimiter // path2_cleaned
   end function join_paths

   !> Separate a full path into path and filename parts. For example:
   !> `separate_path_and_filename('path/to/file.txt', path, filename)` will set `path` to `path/to` and `filename` to `file.txt`.
   subroutine separate_path_and_filename(full_path, path, filename)
      !> Full path to separate
      character(*), intent(in) :: full_path
      !> Path part of the full path
      character(:), allocatable :: path
      !> Filename part of the full path
      character(:), allocatable :: filename

      integer :: last_delimiter_index

      last_delimiter_index = max(0, scan(full_path, '/\'//char(92), back=.true.))
      if (last_delimiter_index > 0) then
         path = full_path(:last_delimiter_index - 1)
         filename = full_path(last_delimiter_index + 1:)
      else
         path = ""
         filename = full_path
      end if
   end subroutine separate_path_and_filename

   !> Split a full path into its components. For example:
   !> `split_path('path/to/file.txt', components)` will set `components` to `['path', 'to', 'file.txt']`.
   subroutine split_path(full_path, components)
      !> Full path to split
      character(*), intent(in) :: full_path
      !> Components of the full path
      character(:), allocatable :: components(:)

      integer :: num_components, i, start_index, end_index
      character(1), parameter :: delimiter = '/'
      character(:), allocatable :: path_cleaned

      ! Remove leading and trailing delimiters
      path_cleaned = adjustl(trim(full_path))
      if (path_cleaned == delimiter) then
         path_cleaned = ""
      else
         if (path_cleaned(1:1) == delimiter) then
            path_cleaned = path_cleaned(2:)
         end if
         if (path_cleaned(len(path_cleaned):len(path_cleaned)) == delimiter) then
            path_cleaned = path_cleaned(:len(path_cleaned) - 1)
         end if
      end if

      ! Count number of components
      num_components = 1
      do i = 1, len(path_cleaned)
         if (path_cleaned(i:i) == delimiter) num_components = num_components + 1
      end do

      allocate(character(len=1) :: components(num_components))
      start_index = 1
      do i = 1, num_components
         end_index = scan(path_cleaned(start_index:), delimiter)
         if (end_index == 0) then
            components(i) = path_cleaned(start_index:)
         else
            components(i) = path_cleaned(start_index:start_index + end_index - 2)
            start_index = start_index + end_index
         end if
      end do
   end subroutine split_path


end module
