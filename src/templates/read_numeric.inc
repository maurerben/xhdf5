! Expects macros:
!   FT_TYPE        : Fortran type, e.g. real(sp)
!   PROCNAME       : subroutine name, e.g. read_real_sp
!   h5fYPE_EXPR() : HDF5 type expression, e.g. hdf5_double()
!   PROCNAME_STR   : quoted string with name, e.g. "read_real_sp"

subroutine PROCNAME(this, h5path, dataset, buffer, datastart, hyperslabs)
  class(hdf5f_t), intent(inout) :: this
  character(*), intent(in) :: h5path
  character(*), intent(in) :: dataset
  FT_TYPE, intent(inout), contiguous, target :: buffer(..)
  integer(int32), intent(in), optional :: datastart(:)
  type(hyperslab_type), intent(in), optional :: hyperslabs(:)

  integer(int32), allocatable :: datasize_file(:), datastart_local(:)
  type(c_ptr) :: buffer_ptr
  type(hyperslab_type), allocatable :: hyperslabs_use(:)
  FT_TYPE, target :: dummy_buffer



  call this%dataset_shape(h5path, dataset, datasize_file, COMPLEX_FLAG) ! We need to compare the datasize HDF5 sees (no complex numbers)
  if(present(hyperslabs)) then ! Use predefined hyperslabs
    hyperslabs_use = hyperslabs
    call assert_true(this%mpi_comm, hyperslabs(1)%check_datasize(datasize_file), &
      'Error(solhdf5%read): hyperslabs have not the same size as dataset ' // h5path // '/' // dataset // '.' )
  else ! Define a hyperslab 
    datastart_local = [int32_undefined]
    if (present(datastart)) datastart_local = datastart

    allocate(hyperslabs_use(1))
    call hyperslabs_use(1)%initialize(shape(buffer, kind=int32), datasize_file, [-1], datastart_local, [-1], [-1], [-1], COMPLEX_FLAG)
  endif

  if(product(shape(buffer, kind=int32)) == 0) then
    buffer_ptr = c_loc(dummy_buffer)
  else
    buffer_ptr = c_loc(buffer)
  end if
  call hdf5_read_dataset(this%mpi_comm, this%file_id, h5path, dataset, h5fYPE_EXPR(), buffer_ptr, &
                         hyperslabs_use, this%serial_access)
end subroutine PROCNAME
