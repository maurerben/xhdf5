module hdf5_dataset

#ifdef _HDF5_
   use hdf5
#endif

   use mpi_utils, only: mpi_comm_type
   use error_handling, only: handle_hdf5_error, assert_true
   use hdf5_globals, only: hdf5_id, hdf5_size, file_id_undefined

   implicit none

   private
   public :: hdf5_get_dataset_shape, hdf5_create_dataset


contains


   !> Get the shape of a dataset.
   subroutine hdf5_get_dataset_shape(mpi_comm, file_id, h5path, dataset, dims)
      !> MPI communicator handle
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the file or group, used by HDF5.
      integer(hdf5_id), intent(in) :: file_id
      !> Path to an object in an HDF5 file. Can be given as `path/to/object`.
      character(*), intent(in) :: h5path
      !> Name of the dataset to obtain shape of.
      character(*), intent(in) :: dataset
      !> Shape of the dataset.
      integer(hdf5_size), allocatable, intent(out) :: dims(:)

      integer :: rank
      integer(hdf5_size), allocatable :: dataset_shape_max(:)

#ifdef _HDF5_
      integer :: h5err
      integer(hdf5_id) :: group_id, dataset_id, dataspace_id, datatype_id

      call assert_true(mpi_comm, file_id /= file_id_undefined, &
         'Error(solhdf5%write): HDF5 file is not initialized.')

      call h5gopen_f(file_id, h5path, group_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gopen_f' ,h5err)

      call h5dopen_f(group_id, dataset, dataset_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5dopen_f', h5err)

      call h5dget_space_f(dataset_id, dataspace_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5dget_space_f', h5err)


      call h5sget_simple_extent_ndims_f(dataspace_id, rank, h5err)
      call handle_hdf5_error(mpi_comm, 'h5dget_space_f', h5err)

      allocate(dims(rank))
      allocate(dataset_shape_max(rank))
      call H5Sget_simple_extent_dims_f(dataspace_id, dims, dataset_shape_max, h5err)
      if (h5err == rank) h5err = 0
      call handle_hdf5_error(mpi_comm, 'H5Sget_simple_extent_dims', h5err)


      call h5sclose_f(dataspace_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5sclose_f', h5err)

      call h5dclose_f(dataset_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5dclose_f', h5err)

      call h5gclose_f(group_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gclose_f', h5err)
#endif
   end subroutine hdf5_get_dataset_shape


   !> Create a new data set without writing to it.
   subroutine hdf5_create_dataset(mpi_comm, file_id, h5path, dataset, dataset_type, dims)
      !> MPI communicator handle
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the file or group, used by HDF5.
      integer(hdf5_id), intent(in) :: file_id
      !> Path to an object in an HDF5 file. Can be given as `path/to/object`.
      character(*), intent(in) :: h5path
      !> Name of the dataset to obtain shape of.
      character(*), intent(in) :: dataset
      !> Type id of the data set.
      integer(hdf5_id), intent(in) :: dataset_type
      !> Shape of the dataset.
      integer(hdf5_size), intent(in) :: dims(:)

      integer :: rank

#ifdef _HDF5_
      integer :: h5err
      integer(hdf5_id) :: group_id, dataset_id, dataspace_id

      rank = size(dims)

      call assert_true(mpi_comm, file_id /= file_id_undefined, &
         'Error(solhdf5%write): HDF5 file is not initialized.')

      call h5gopen_f(file_id, h5path, group_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gopen_f' ,h5err)

      call h5screate_simple_f(rank, dims, dataspace_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5screate_simple_f', h5err)

      call h5dcreate_f(group_id, dataset, dataset_type, dataspace_id, dataset_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5dcreate', h5err)

      call h5sclose_f(dataspace_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5sclose_f', h5err)

      call h5dclose_f(dataset_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5dclose_f', h5err)

      call h5gclose_f(group_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gclose_f', h5err)
#endif
   end subroutine

end module hdf5_dataset
