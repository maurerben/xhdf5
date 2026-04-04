!> Wrappers for initializing/finalizing HDF5 dependent globals and creating/opening/closing HDF5 files.
module hdf5_file

#ifdef _HDF5_
   use hdf5
#endif

   use error_handling, only: handle_hdf5_error
   use hdf5_globals, only: hdf5_id
   use mpi_utils, only: mpi_comm_type, comm_to_fint, mpi_info_null_fint


   implicit none


   private
   public :: hdf5_initialize, hdf5_finalize, hdf5_create_file, hdf5_open_file, hdf5_close_file


contains


   !> init global variables, used by HDF5 library functions.
   subroutine hdf5_initialize(mpi_comm)
      !> MPI communicator handle
      type(mpi_comm_type), intent(in) :: mpi_comm
#ifdef _HDF5_
      integer :: h5err
      call h5open_f(h5err)
      call handle_hdf5_error(mpi_comm, 'h5open_f', h5err)
#endif
   end subroutine hdf5_initialize

   !> Finalize global variables, used by HDF5 library functions.
   subroutine hdf5_finalize(mpi_comm)
      !> MPI communicator handle
      type(mpi_comm_type), intent(in) :: mpi_comm
#ifdef _HDF5_
      integer :: h5err

      call h5close_f(h5err)
      call handle_hdf5_error(mpi_comm, 'h5_close_f', h5err)
#endif
   end subroutine

   !> Create an HDF5 file at `path`.
   subroutine hdf5_create_file(path, mpi_comm, file_id, serial_access)
      !> Relative path to the HDF5 file. Is expected to end with the name of the file,
      !> _e.g._, `path = 'path/to/file.h5'. HDF5 files are expected to have the `.h5` suffix.
      character(*), intent(in) :: path
      !> MPI communicator.
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the file, used by HDF5.
      integer(hdf5_id), intent(out) :: file_id
      !> Set to `.true.` if only serial access is possible.
      logical, intent(in) :: serial_access

#ifdef _HDF5_
      integer :: h5err
      integer(hdf5_id) :: file_id_plist

      call h5pcreate_f(H5P_FILE_ACCESS_F, file_id_plist, h5err)
      call handle_hdf5_error(mpi_comm, 'h5pcreate_f', h5err)

#ifdef MPI
      if(.not. serial_access) then
         call h5pset_fapl_mpio_f(file_id_plist, comm_to_fint(mpi_comm), mpi_info_null_fint(), h5err)
         call handle_hdf5_error(mpi_comm, 'h5pset_fapl_mpio_f', h5err)
      end if
#endif

      call h5fcreate_f(path, H5F_ACC_TRUNC_F, file_id, h5err, access_prp = file_id_plist)
      call handle_hdf5_error(mpi_comm, 'h5fopen_f', h5err)

      call h5pclose_f(file_id_plist, h5err)
      call handle_hdf5_error(mpi_comm, 'h5pclose_f', h5err)
#endif
   end subroutine hdf5_create_file


   !> Open an HDF5 file at `path`.
   subroutine hdf5_open_file(path, mpi_comm, file_id, serial_access)
      !> Relative path to the HDF5 file. Is expected to end with the name of the file,
      !> _e.g._, `path = 'path/to/file.h5'. HDF5 files are expected to have the `.h5` suffix.
      character(*), intent(in) :: path
      !> MPI communicator.
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the file, used by HDF5.
      integer(hdf5_id), intent(out) :: file_id
      !> Set to `.true.` if only serial access is possible.
      logical, intent(in) :: serial_access

#ifdef _HDF5_
      integer :: h5err
      integer(hdf5_id) :: file_id_plist

      call h5pcreate_f(H5P_FILE_ACCESS_F, file_id_plist, h5err)
      call handle_hdf5_error(mpi_comm, 'h5pcreate_f', h5err)

#ifdef MPI
      if(.not. serial_access) then
         call h5pset_fapl_mpio_f(file_id_plist, comm_to_fint(mpi_comm), mpi_info_null_fint(), h5err)
         call handle_hdf5_error(mpi_comm, 'h5pset_fapl_mpio_f', h5err)
      end if
#endif

      call h5fopen_f(path, H5F_ACC_RDWR_F, file_id, h5err, access_prp = file_id_plist)
      call handle_hdf5_error(mpi_comm, 'h5fopen_f', h5err)

      call h5pclose_f(file_id_plist, h5err)
      call handle_hdf5_error(mpi_comm, 'h5pclose_f', h5err)
#endif
   end subroutine hdf5_open_file

   !> Close an HDF5 file with the corresponding identifier.
   subroutine hdf5_close_file(mpi_comm, file_id)
      !> MPI communicator handle
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the file, used by HDF5.
      integer(hdf5_id), intent(in) :: file_id
#ifdef _HDF5_
      integer :: h5err

      call h5fclose_f(file_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5fclose_f', h5err)
#endif
   end subroutine hdf5_close_file

end module hdf5_file
