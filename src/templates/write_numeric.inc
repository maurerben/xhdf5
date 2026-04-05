
! This file expects these macros to be defined before inclusion:
!   FT_TYPE        : Fortran type, e.g. real(sp)
!   FT_SUFFIX      : suffix token, e.g. sp
!   h5fYPE_EXPR() : expression yielding the HDF5 type id/tag
!   PROCNAME       : subroutine name, e.g. write_real_sp
!   COMPLEX_FLAG   : .true. / .false. (if you need it anywhere)

subroutine PROCNAME(this, h5path, dataset, buffer, datastart, datasize, hyperslabs)
  class(h5file_t), intent(inout) :: this
  character(*), intent(in) :: h5path
  character(*), intent(in) :: dataset
  FT_TYPE, intent(in), contiguous, target :: buffer(..)
  integer(int32), intent(in), optional :: datastart(:)
  integer(int32), intent(in), optional :: datasize(:)
  type(hyperslab_type), intent(in), optional :: hyperslabs(:)

  integer(int32), allocatable :: datasize_local(:), datastart_local(:)
  type(hyperslab_type), allocatable :: hyperslabs_use(:)
  type(c_ptr) :: buffer_ptr
  FT_TYPE, target :: dummy_buffer

  if(this%serial_access) then
     call assert_true(this%mpi_comm, comm_to_rank(this%mpi_comm) == root_rank, &
        'Error(solhdf5%write): serial_access is set to .true., only the root process &
        is allowed to call solhdf5 routines.')
  end if

  if(present(hyperslabs)) then ! Use predefined hyperslabs
    hyperslabs_use = hyperslabs
  else ! Define a hyperslab 
    datastart_local = [int32_undefined]
    if (present(datastart)) datastart_local = datastart
    datasize_local = [int32_undefined]
    if (present(datasize)) datasize_local = datasize
    allocate(hyperslabs_use(1))
    call hyperslabs_use(1)%init(shape(buffer, kind=int32), datasize_local, [-1], datastart_local, [-1], [-1], [-1], COMPLEX_FLAG)
  endif

  call this%handle_if_dataset_exists(h5path, dataset)
  if(product(shape(buffer, kind=int32)) == 0) then
    buffer_ptr = c_loc(dummy_buffer)
  else
    buffer_ptr = c_loc(buffer)
  end if
  call hdf5_write_dataset(this%mpi_comm, this%file_id, h5path, dataset, h5fYPE_EXPR(), buffer_ptr, &
                          hyperslabs_use, this%serial_access)
end subroutine PROCNAME
