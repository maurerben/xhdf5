!> Module with utility functions for simple os commands.
module os_utils
  implicit none
  private

  public :: path_exists, join_paths

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

end module os_utils