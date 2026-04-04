module hdf5_datatype

#ifdef _HDF5_
   use hdf5
#endif

   use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
   use error_handling, only: handle_hdf5_error, assert_true
   use hdf5_globals, only: hdf5_id, hdf5_size, file_id_undefined, h5size_undefined
   use mpi_utils, only: mpi_comm_type

   implicit none

   private
   public :: hdf5_double, hdf5_float, hdf5_integer_int32, hdf5_string, hdf5_get_type, hdf5_get_type_size, hdf5_datatypes_equal

   integer(int32) :: str_16 = 16

contains


   !> Return HDF5 precision for double
   integer(hdf5_id) function hdf5_double()
#ifdef _HDF5_
      hdf5_double = H5T_NATIVE_DOUBLE
#else
      hdf5_double = real64
#endif
   end function

   !> Return HDF5 precision for real (real32)
   integer(hdf5_id) function hdf5_float()
#ifdef _HDF5_
      hdf5_float = H5T_NATIVE_REAL
#else
      hdf5_float = real32
#endif
   end function hdf5_float


   !> Return HDF5 precision for integer
   integer(hdf5_id) function hdf5_integer_int32()
#ifdef _HDF5_
      hdf5_integer_int32 = H5T_NATIVE_INTEGER
#else
      hdf5_integer_int32 = int32
#endif
   end function hdf5_integer_int32


   !> Return HDF5 precision for character
   integer(hdf5_id) function hdf5_string(mpi_comm, string_len)
      !> MPI communicator handle
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Length of the string.
      integer(hdf5_size), intent(in) :: string_len

      integer(real64) :: string_len_local
      integer :: h5err

      string_len_local = string_len
#ifdef _HDF5_
      call h5tcopy_f(H5T_NATIVE_CHARACTER, hdf5_string, h5err)
      call handle_hdf5_error(mpi_comm, 'h5tcopy_f', h5err)

      call h5tset_size_f(hdf5_string, string_len_local, h5err)
      call handle_hdf5_error(mpi_comm, 'h5tset_size_f', h5err)
#else
      hdf5_string = str_16
#endif
   end function hdf5_string


   !> Get the datatype of a dataset
   integer(hdf5_id) function hdf5_get_type(mpi_comm, file_id, h5path, dataset)
      !> MPI communicator handle
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the file or group, used by HDF5.
      integer(hdf5_id), intent(in) :: file_id
      !> Path to an object in an HDF5 file. Can be given as `path/to/object`.
      character(*), intent(in) :: h5path
      !> Name of the dataset to obtain shape of.
      character(*), intent(in) :: dataset

#ifdef _HDF5_
      integer :: h5err
      integer(hdf5_id) :: group_id, dataset_id, file_id_dtype

      call assert_true(mpi_comm, file_id /= file_id_undefined, &
         'Error(solhdf5%write): HDF5 file is not initialized.')

      call h5gopen_f(file_id, h5path, group_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gopen_f' ,h5err)

      call h5dopen_f(group_id, dataset, dataset_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5dopen_f', h5err)

      call h5dget_type_f(dataset_id, hdf5_get_type, h5err)
      call handle_hdf5_error(mpi_comm, 'h5dget_type_f', h5err)

      call h5dclose_f(dataset_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5dclose_f', h5err)

      call h5gclose_f(group_id, h5err)
      call handle_hdf5_error(mpi_comm, 'h5gclose_f', h5err)
#else
      hdf5_get_type = file_id_undefined
#endif
   end function hdf5_get_type


   !> Get the size of a data precision.
   integer(hdf5_size) function hdf5_get_type_size(mpi_comm, file_id_type)
      !> MPI communicator handle
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the datatype.
      integer(hdf5_id), intent(in) :: file_id_type

#ifdef _HDF5_
      integer :: h5err

      call h5tget_size_f(file_id_type, hdf5_get_type_size, h5err)
      call handle_hdf5_error(mpi_comm, 'h5tget_size_f', h5err)
#else
      hdf5_get_type_size = h5size_undefined
#endif
   end function hdf5_get_type_size


   !> Determine if two datatypes are equal.
   logical function hdf5_datatypes_equal(mpi_comm, file_id_type1, file_id_type2)
      !> MPI communicator handle
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Identifier of the datatypes that are compared.
      integer(hdf5_id), intent(in) :: file_id_type1, file_id_type2

#ifdef _HDF5_
      integer :: h5err

      call h5tequal_f(file_id_type1, file_id_type1, hdf5_datatypes_equal, h5err)
      call handle_hdf5_error(mpi_comm, 'h5tequal_f', h5err)
#else
      hdf5_datatypes_equal = .false.
#endif
   end function hdf5_datatypes_equal

end module
