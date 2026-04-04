!> Wrapper for writing any data set to an HDF5 file
module hdf5_write
   use iso_c_binding, only: c_ptr

#ifdef _HDF5_
   use hdf5
#endif

   use error_handling, only: handle_hdf5_error, assert_true
   use hdf5_globals, only: hdf5_id, hdf5_size, hdf5_ssize, file_id_undefined
   use mpi_utils, only: mpi_comm_type
   use hyperslab, only: hyperslab_type, check_hyperslabs


   implicit none


   private
   public :: hdf5_write_dataset

   ! interface hdf5_write_dataset
   !   module procedure :: hdf5_write_dataset_memspace, hdf5_write_dataset_hyperslab
   ! end interface


contains

   subroutine hdf5_write_dataset(mpi_comm, file_id, h5path, dataset, dataset_type, buffer_ptr, hyperslabs, serial_access)
      !> MPI communicator handle
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the file or group, used by HDF5.
      integer(hdf5_id), intent(in) :: file_id
      !> Absolute path in the HDF5 file to the group to write the data set to.
      character(*), intent(in) :: h5path
      !> Name of the data set to write.
      character(*), intent(in) :: dataset
      !> Type id of the data set.
      integer(hdf5_id), intent(in) :: dataset_type
      !> Pointer to the first element of the data set to be written.
      type(c_ptr), intent(inout) :: buffer_ptr
      !> Type that holds information about the global and local data set shape and datastart.
      type(hyperslab_type) :: hyperslabs(:)
      !> Set to `.true.` if only serial access is possible.
      logical, intent(in) :: serial_access

#ifdef _HDF5_
      type(hyperslab_type) :: hyperslab
      logical :: overwrite_local
      integer(hdf5_id) :: group_id, dataspace_id, dataset_id, memspace_id, file_id_plist
      integer :: h5err, hsdx, hs_operator
      integer :: rank
      integer(hdf5_size), allocatable :: memsize(:), datasize(:), memsize_use(:)
      logical :: has_selection

      call assert_true(mpi_comm, file_id /= file_id_undefined, &
         'Error(hdf5_write_dataset_memspace): HDF5 file is not initialized.')

      call check_hyperslabs(hyperslabs, 'hdf5_write_dataset', mpi_comm)

      rank = hyperslabs(1)%rank
      memsize = hyperslabs(1)%memsize
      datasize = hyperslabs(1)%datasize
      memsize_use = max(memsize, int(1, hdf5_size))
      ! Todo(Bene): assert that hyperslabs makes sense

      ! Create dataspace
      call h5screate_simple_f(rank, datasize, dataspace_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5screate_simple_f', h5err)

      ! Create memory space
      call h5screate_simple_f(rank, memsize_use, memspace_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5screate_simple_f', h5err)

      has_selection = .false.
      do hsdx=1, size(hyperslabs)
         if(any(hyperslabs(hsdx)%count == 0)) cycle

         if(hsdx==1) then
            hs_operator = H5S_SELECT_SET_F
         else
            hs_operator = H5S_SELECT_OR_F ! Join hyperslabs
         end if

         ! dataspace
         call h5sselect_hyperslab_f(dataspace_id, hs_operator, hyperslabs(hsdx)%datastart, hyperslabs(hsdx)%count, h5err, &
            hyperslabs(hsdx)%stride, hyperslabs(hsdx)%datablock)
         call handle_hdf5_error(mpi_comm, 'h5sselect_hyperslab_f:dataspace', h5err)

         ! memspace
         call h5sselect_hyperslab_f(memspace_id, hs_operator, hyperslabs(hsdx)%memstart, hyperslabs(hsdx)%count, h5err, &
            hyperslabs(hsdx)%stride, hyperslabs(hsdx)%datablock)
         call handle_hdf5_error(mpi_comm, 'h5sselect_hyperslab_f:memspace', h5err)
         has_selection = .true.
      end do

      if(.not. has_selection) then
         call h5sselect_none_f(dataspace_id, h5err)
         call handle_hdf5_error(mpi_comm, 'h5sselect_none_f:dataspace', h5err)
         call h5sselect_none_f(memspace_id, h5err)
         call handle_hdf5_error(mpi_comm, 'h5sselect_none_f:memspace', h5err)
      end if

      ! Create a property list to enable parallel writing, if MPI is used.
      call h5pcreate_f(H5P_DATASET_XFER_F, file_id_plist, h5err)
      call handle_hdf5_error(mpi_comm, 'h5pcreate_f', h5err)

#ifdef MPI
      if(.not. serial_access) then
         call h5pset_dxpl_mpio_f(file_id_plist, H5FD_MPIO_COLLECTIVE_F, h5err)
         call handle_hdf5_error(mpi_comm, 'h5pset_dxpl_mpio_f', h5err)
      end if
#endif

      ! Open group to write the data set in
      call h5gopen_f(file_id, h5path, group_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gopen_f' ,h5err)

      ! Create dataset in file
      call h5dcreate_f(group_id, dataset, dataset_type, dataspace_id, dataset_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5dcreate', h5err)

      ! Write data to file
      call h5dwrite_f(dataset_id, dataset_type, buffer_ptr, h5err, &
         file_space_id = dataspace_id, mem_space_id = memspace_id, xfer_prp = file_id_plist)
      call handle_hdf5_error(mpi_comm, 'h5dwrite_f', h5err)

      ! Close open objects
      call h5dclose_f(dataset_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5dclose_f', h5err)

      call h5gclose_f(group_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gclose_f', h5err)

      call h5pclose_f(file_id_plist, h5err)
      call handle_hdf5_error(mpi_comm, 'h5pclose_f', h5err)

      call h5sclose_f(dataspace_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5sclose_f', h5err)

      call h5sclose_f(memspace_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5sclose_f', h5err)
#endif
   end subroutine hdf5_write_dataset

end module hdf5_write
