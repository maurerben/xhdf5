module hdf5_group

#ifdef _HDF5_
   use hdf5
#endif

   use error_handling, only: handle_hdf5_error, assert_true
   use hdf5_globals, only: hdf5_id, file_id_undefined
   use mpi_utils, only: mpi_comm_type


   implicit none


   private
   public :: hdf5_create_group, hdf5_link_exists, hdf5_delete_link


contains


   !> Create a new `group` at `h5path`, associated to a link in an HDF5 file or group associated with `file_id`.
   subroutine hdf5_create_group(mpi_comm, file_id, h5path, group)
      !> MPI communicator.
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the file or group, used by HDF5.
      integer(hdf5_id), intent(in) :: file_id
      !> Path to an object in an HDF5 file. Can be given as `path/to/object`.
      character(*), intent(in) :: h5path
      !> Name of the group to create.
      character(*), intent(in) :: group
#ifdef _HDF5_
      integer :: h5err
      integer(hdf5_id) :: group_id, file_id_newgroup

      call assert_true(mpi_comm, file_id /= file_id_undefined, &
         'Error(xhdf5%write): HDF5 file is not initialized.')

      call h5gopen_f(file_id, h5path, group_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gopen_f' ,h5err)

      call h5gcreate_f(group_id, group, file_id_newgroup, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gcreate_f' ,h5err)

      call h5gclose_f(file_id_newgroup, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gclose_f', h5err)

      call h5gclose_f(group_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gclose_f', h5err)
#endif
   end subroutine hdf5_create_group


   !> Check if `h5path` is associated with a link in an HDF5 file or group corresponding to `file_id`.
   logical function hdf5_link_exists(mpi_comm, file_id, h5path)
      !> MPI communicator.
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the file or group, used by HDF5.
      integer(hdf5_id), intent(in) :: file_id
      !> Path to an object in an HDF5 file. Can be given as `path/to/object`.
      character(*), intent(in) :: h5path
#ifdef _HDF5_
      integer :: h5err

      call assert_true(mpi_comm, file_id /= file_id_undefined, &
         'Error(xhdf5%write): HDF5 file is not initialized.')

      call h5lexists_f(file_id, h5path, hdf5_link_exists, h5err)
      call handle_hdf5_error(mpi_comm, 'h5lexists_f' ,h5err)
#else
      hdf5_link_exists = .false.
#endif
   end function hdf5_link_exists


   !> Delete object associated with `h5path` in an HDF5 file or group corresponding to `file_id`.
   subroutine hdf5_delete_link(mpi_comm, file_id, h5path)
      !> MPI communicator.
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the file or group, used by HDF5.
      integer(hdf5_id), intent(in) :: file_id
      !> Path to an object in an HDF5 file. Can be given as `path/to/object`.
      character(*), intent(in) :: h5path
#ifdef _HDF5_
      integer :: h5err

      call assert_true(mpi_comm, file_id /= file_id_undefined, &
         'Error(xhdf5%write): HDF5 file is not initialized.')

      call h5ldelete_f(file_id, h5path, h5err)
      call handle_hdf5_error(mpi_comm, 'h5ldelete_f' ,h5err)
#endif
   end subroutine hdf5_delete_link

end module hdf5_group
