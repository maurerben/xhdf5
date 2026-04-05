!> Module for HDF5 wrappers.
!> Easy to use fortran wrappers for the HDF5 library.
!> Parallel MPI writing of array strides is supported.
!> Use these routines for writing large binary files that can be easily reused on other mashines.
!>
!> The use of this module requires to build exciting with HDF5.
!> If exciting is compiled without HDF5 the routines will do nothing. If your feature requires HDF5,
!> Make sure to call `[[abort_if_not_hdf5]]` in the beginning of the code to make sure that exciting gets
!> safely terminated if built without HDF5.
!>
!> As already mentioned, HDF5 is a library to write large binary files, especially writing large arrays.
!> Therefore it does absolutely make no sense to use them in loops.
!> If you are tempted to do this, you should seriously rethink about your data layout!
!>
!> For examples how to use the library we refer kindly to the unit tests (`[[solhdf5_test.f90]]).
!> These show all examples of usage.
!>
!> Be aware of that the writing routines will overwrite an existing data set with the same name.
module solhdf5
   use iso_c_binding, only: c_ptr, c_loc
   use iso_fortran_env, only: real32, real64, int32

   use os_utils, only: path_exists, join_paths

   use hdf5_globals
   use error_handling
#ifdef MPI
   use mpi_f08, only: MPI_Comm
#endif
   use mpi_utils
   use hdf5_file
   use hdf5_group
   use hdf5_write
   use hdf5_read
   use hdf5_dataset
   use hdf5_datatype
   use hyperslab


   implicit none


   private
   public :: h5file_t, h5group_root, hyperslab_type


   !> Type for handeling an HDF5 file.
   type h5file_t
      !> OS file path.
      character(:), allocatable :: path
      !> HDF5 file id.
      integer(hdf5_id) :: file_id = file_id_undefined
      !> MPI communicator handle
      type(mpi_comm_type) :: mpi_comm
      !> Flag to init HDF5 file in serial mode, even if compiled with MPI.
      !> By default this is set to false.
      logical :: serial_access = .false.
      !> Number of overwritten datasets
      integer :: n_overwritten_datasets = 0

   contains

#ifdef MPI
      generic :: init => initialize_mpi_comm_world,initialize_mpi_comm_type, initialize_mpi_comm, initialize_mpi_comm_f08
      procedure :: initialize_mpi_comm_world, initialize_mpi_comm_type, initialize_mpi_comm, initialize_mpi_comm_f08
#else
      generic :: init => initialize_mpi_comm_type, initialize_mpi_comm
      procedure :: initialize_mpi_comm_type, initialize_mpi_comm
#endif

      procedure :: delete, initialize_group, initialize_group_update_groupname, &
         dataset_shape, link_exists, delete_link

      procedure :: handle_if_dataset_exists
      procedure :: evaluate_overwritten_datasets

      generic :: write => write_string, write_bool, write_integer_int32, write_real_real32, &
         write_real_real64, write_complex_real32, write_complex_real64

      procedure :: write_string, write_bool, write_integer_int32, write_real_real32, &
         write_real_real64, write_complex_real32, write_complex_real64

      generic :: read => read_string, read_bool, read_integer_int32, read_real_real32, &
         read_real_real64, read_complex_real32, read_complex_real64

      procedure :: read_string, read_bool, read_integer_int32, read_real_real32, &
         read_real_real64, read_complex_real32, read_complex_real64

   end type h5file_t


contains

   !> init HDF5 library and, if `path` exists, open the file, else create a file at `path`.
   subroutine initialize_mpi_comm_world(this, path, serial_access)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Relative path to HDF5 file.
      character(*), intent(in) :: path
      !> Set to `.true.` for serial access in an MPI environment.
      !> With this setting enabled, only the root process is allowed to call solhdf5 routines.
      logical, intent(in), optional :: serial_access

      this%path = trim(path)

      call init_with_mpi_comm_world(this%mpi_comm)

      if(present(serial_access)) then
         this%serial_access = serial_access
      end if

      if(this%serial_access) then
         call assert_true(this%mpi_comm, comm_to_rank(this%mpi_comm) == root_rank, &
            'Error(solhdf5%initialize_mpi_comm_type): serial_access is set to .true., only the root process &
            is allowed to call solhdf5 routines.')
      end if

      call hdf5_initialize(this%mpi_comm)

      if(path_exists(this%path)) then
         call hdf5_open_file(this%path, this%mpi_comm, this%file_id, this%serial_access)
      else
         call hdf5_create_file(this%path, this%mpi_comm, this%file_id, this%serial_access)
      end if
   end subroutine

   !> init HDF5 library and, if `path` exists, open the file, else create a file at `path`.
   subroutine initialize_mpi_comm_type(this, path, mpi_comm, serial_access)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Relative path to HDF5 file.
      character(*), intent(in) :: path
      !> MPI communicator.
      type(mpi_comm_type), intent(in) :: mpi_comm
      !> Set to `.true.` for serial access in an MPI environment.
      !> With this setting enabled, only the root process is allowed to call solhdf5 routines.
      logical, intent(in), optional :: serial_access

      this%path = trim(path)
      this%mpi_comm = mpi_comm

      if(present(serial_access)) then
         this%serial_access = serial_access
      end if

      if(this%serial_access) then
         call assert_true(this%mpi_comm, comm_to_rank(this%mpi_comm) == root_rank, &
            'Error(solhdf5%initialize_mpi_comm_type): serial_access is set to .true., only the root process &
            is allowed to call solhdf5 routines.')
      end if

      call hdf5_initialize(this%mpi_comm)

      if(path_exists(this%path)) then
         call hdf5_open_file(this%path, this%mpi_comm, this%file_id, this%serial_access)
      else
         call hdf5_create_file(this%path, this%mpi_comm, this%file_id, this%serial_access)
      end if
   end subroutine

   !> init HDF5 library with MPI_Comm from mpi_f08.
   !> If `path` exists, open the file, else create a file at `path`.
   subroutine initialize_mpi_comm(this, path, mpi_comm, serial_access)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Relative path to HDF5 file.
      character(*), intent(in) :: path
      !> MPI integer communicator
      integer, intent(in) :: mpi_comm
      !> Set to `.true.` for serial access in an MPI environment.
      !> With this setting enabled, only the root process is allowed to call solhdf5 routines.
      logical, intent(in), optional :: serial_access

      type(mpi_comm_type) :: mpi_comm_use 

      ! Create mpi_comm_type from MPI_Comm
      call init_with_mpi_fint_comm(mpi_comm_use, mpi_comm)

      if (present(serial_access)) then
         call this%initialize_mpi_comm_type(path, mpi_comm_use, serial_access)
      else
         call this%initialize_mpi_comm_type(path, mpi_comm_use)
      end if
   end subroutine initialize_mpi_comm

#ifdef MPI
   !> init HDF5 library with MPI_Comm from mpi_f08.
   !> If `path` exists, open the file, else create a file at `path`.
   subroutine initialize_mpi_comm_f08(this, path, mpi_comm_f08, serial_access)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Relative path to HDF5 file.
      character(*), intent(in) :: path
      !> MPI communicator from mpi_f08.
      type(MPI_Comm), intent(in) :: mpi_comm_f08
      !> Set to `.true.` for serial access in an MPI environment.
      !> With this setting enabled, only the root process is allowed to call solhdf5 routines.
      logical, intent(in), optional :: serial_access

      type(mpi_comm_type) :: mpi_comm_use

      ! Create mpi_comm_type from MPI_Comm
      mpi_comm_use%handle = mpi_comm_f08

      if (present(serial_access)) then
         call this%initialize_mpi_comm_type(path, mpi_comm_use, serial_access)
      else
         call this%initialize_mpi_comm_type(path, mpi_comm_use)
      end if
   end subroutine initialize_mpi_comm_f08
#endif

   !> Close file and finalize HDF5 library.
   subroutine delete(this)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this

      if (this%file_id == file_id_undefined) return ! File was never initialized

      call hdf5_close_file(this%mpi_comm, this%file_id)
      call hdf5_finalize(this%mpi_comm)
      call this%evaluate_overwritten_datasets()
      this%file_id = file_id_undefined
   end subroutine delete


   !> Create a new group with name `group` at `h5path`. If the group already exists, the routine does nothing.
   subroutine initialize_group(this, h5path, groupname)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Absolute path in the hdf5 file.
      character(*), intent(in) :: h5path
      !> Group name.
      character(*), intent(in) :: groupname

      if (this%link_exists(join_paths(h5path, groupname))) return
      call hdf5_create_group(this%mpi_comm, this%file_id, trim(h5path), groupname)
   end subroutine initialize_group


   !> Create a new group with name `group` at `h5path`. If the group already exists, the routine does nothing.
   subroutine initialize_group_update_groupname(this, h5path, groupname)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Absolute path in the hdf5 file.
      character(*), intent(in) :: h5path
      !> Group name.
      character(:), allocatable, intent(inout) :: groupname

      call this%initialize_group(h5path, groupname)
      groupname = join_paths(h5path, groupname)
   end subroutine initialize_group_update_groupname

   !> Return true or false, whether the link to `h5path` exists or not.
   logical function link_exists(this, h5path)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Absolute path in the hdf5 file.
      character(*), intent(in) :: h5path

      link_exists = hdf5_link_exists(this%mpi_comm, this%file_id, trim(h5path))
   end function link_exists

   !> Delete link to `h5path`. This action frees the space on the drive allocated by the object `h5path` points to.
   subroutine delete_link(this, h5path)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Absolute path in the hdf5 file.
      character(*), intent(in) :: h5path

      call hdf5_delete_link(this%mpi_comm, this%file_id, trim(h5path))
   end subroutine delete_link


   !> Get the shape of a data set.
   subroutine dataset_shape(this, h5path, datasetname, shape, complex_dataset)
      !> HDF5 file handler.
      class(h5file_t), intent(in) :: this
      !> Absolute path in the hdf5 file.
      character(*), intent(in) :: h5path
      !> Dataset name.
      character(*), intent(in) :: datasetname
      !> Dataset shape
      integer, allocatable :: shape(:)
      !> Set to `.true.` if the data set is complex. If not, the first dimension will be
      !> always 2 for complex datasets.
      logical, optional :: complex_dataset

      logical :: complex_dataset_local
      integer(hdf5_size), allocatable :: shape_local(:)

      complex_dataset_local = .false.
      if(present(complex_dataset)) complex_dataset_local = complex_dataset

      call hdf5_get_dataset_shape(this%mpi_comm, this%file_id, h5path, datasetname, shape_local)
      if(complex_dataset_local) then
         shape = shape_local(2:)
      else
         shape = shape_local
      end if
   end subroutine dataset_shape


   !> Warn about overwritten datasets.
   subroutine evaluate_overwritten_datasets(this)
      !> HDF5 file handler.
      class(h5file_t), intent(in) :: this

#ifdef USE_ASSERT
      if(this%n_overwritten_datasets > 0) then
         write(error_unit, '(1x, A, I3, A, A, A)') &
            "Warning(solhdf5): ", this%n_overwritten_datasets, " dataset(s) were overwritten in file ", this%path, "."
      end if
#endif
   end subroutine evaluate_overwritten_datasets


   !> Handles what to do if a data set already exists.
   !> If this is the case, the data set will be overwritten with new data.
   !> In debug mode, a warning is printed to the terminal.
   subroutine handle_if_dataset_exists(this, h5path, datasetname)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Absolute path in the hdf5 file.
      character(*), intent(in) :: h5path
      !> Dataset name.
      character(*), intent(in) :: datasetname

      if (this%link_exists(join_paths(h5path, datasetname))) then
         this%n_overwritten_datasets = this%n_overwritten_datasets + 1
         call this%delete_link(join_paths(h5path, datasetname))

#ifdef USE_ASSERT
         write(error_unit, '(1x, A, A, A)') &
            "Warning(solhdf5): Dataset ", datasetname, " already existed and was deleted. &
            It will be overwritten with new data."
#endif
      end if

   end subroutine

!-----------------------------------------------------------------------------------------------------------------------
! WRITE DATASET

   !> Write a string to an HDF5 file.
   subroutine write_string(this, h5path, dataset, string)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Absolute path in the HDF5 file to the group to write the data set to.
      character(*), intent(in) :: h5path
      !> Name of the data set to write.
      character(*), intent(in) :: dataset
      !> String to write.
      character(*), intent(in), target :: string

      integer(hdf5_id) :: string_type
      type(c_ptr) :: buffer_ptr
      type(hyperslab_type) :: hyperslab
      character(len=7), target :: empty_marker

      call hyperslab%init([1], [-1], [-1], [-1], [-1], [-1], [-1], .false.)
      call this%handle_if_dataset_exists(h5path, dataset)

      ! Handle empty strings specially
      if (len(string) == 0) then
         empty_marker = "<EMPTY>"
         string_type = hdf5_string(this%mpi_comm, len(empty_marker, kind=hdf5_size))
         buffer_ptr = c_loc(empty_marker)
      else
         string_type = hdf5_string(this%mpi_comm, len(string, kind=hdf5_size))
         buffer_ptr = c_loc(string)
      end if

      call hdf5_write_dataset(this%mpi_comm, this%file_id, trim(h5path), dataset, string_type, &
         buffer_ptr, [hyperslab], this%serial_access)
   end subroutine write_string


   !> Write a boolean to an HDF5 file.
   !> The routine does not really writes a boolean but a string representing the boolean.
   !> For '.true.' it writes `'true'` and for `.false.` `'false'`.
   subroutine write_bool(this, h5path, dataset, bool)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Absolute path in the HDF5 file to the group to write the data set to.
      character(*), intent(in) :: h5path
      !> Name of the data set to write.
      character(*), intent(in) :: dataset
      !> Boolean to write.
      logical, intent(in), target :: bool

      character(:), allocatable :: bool_string

      if(bool) then
         bool_string = true_string
      else
         bool_string = false_string
      end if

      call write_string(this, h5path, dataset, bool_string)
   end subroutine write_bool

   !> Read a character scalar from an HDF5 file.
   subroutine read_string(this, h5path, dataset, string)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Absolute path in the HDF5 file to the group to read the data set from.
      character(*), intent(in) :: h5path
      !> Name of the data set to read from.
      character(*), intent(in) :: dataset
      !> String to read.
      character(:), allocatable, intent(out), target :: string

      integer, parameter :: dataset_rank = 1

      integer(hdf5_size) :: string_len
      integer(hdf5_id) :: dataset_type
      type(c_ptr) :: buffer_ptr
      type(hyperslab_type) :: hyperslab

      dataset_type = hdf5_get_type(this%mpi_comm, this%file_id, h5path, dataset)
      allocate(character(len=hdf5_get_type_size(this%mpi_comm, dataset_type)) :: string)

      call hyperslab%init([1], [-1], [-1], [-1], [-1], [-1], [-1], .false.)

      buffer_ptr = c_loc(string)
      call hdf5_read_dataset(this%mpi_comm, this%file_id, h5path, dataset, dataset_type, buffer_ptr, [hyperslab], this%serial_access)

      ! Handle empty strings specially
      if (string == "<EMPTY>") then
         deallocate(string)
         allocate(character(len=0) :: string)
      end if
   end subroutine read_string

   !> Read a boolean to an HDF5 file.
   !> The boolean is not saved as a boolean in the file but as a string representing the state
   !> of the boolean. See [[read_bool]].
   subroutine read_bool(this, h5path, dataset, bool)
      !> HDF5 file handler.
      class(h5file_t), intent(inout) :: this
      !> Absolute path in the HDF5 file to the group to write the data set to.
      character(*), intent(in) :: h5path
      !> Name of the data set to write.
      character(*), intent(in) :: dataset
      !> Boolean to read.
      logical, intent(out), target :: bool

      character(:), allocatable :: bool_string

      call read_string(this, h5path, dataset, bool_string)

      if (bool_string == true_string) then
         bool = .true.
      else if(bool_string == false_string) then
         bool = .false.
      else
         call assert_true(this%mpi_comm, .false., 'dataset does not contain a string representing the state of a bool.')
      end if
   end subroutine read_bool

   ! Numerical interfaces

   ! ---- integer(int32)
#define FT_TYPE        integer(int32)
#define h5fYPE_EXPR() hdf5_integer_int32()
#define COMPLEX_FLAG   .false.

#define PROCNAME       write_integer_int32
#include "templates/write_numeric.inc"
#undef PROCNAME
#define PROCNAME       read_integer_int32
#include "templates/read_numeric.inc"
#undef PROCNAME

#undef FT_TYPE
#undef h5fYPE_EXPR
#undef COMPLEX_FLAG

! ---- real(real32)
#define FT_TYPE        real(real32)
#define h5fYPE_EXPR() hdf5_float()
#define COMPLEX_FLAG   .false.

#define PROCNAME       write_real_real32
#include "templates/write_numeric.inc"
#undef PROCNAME
#define PROCNAME       read_real_real32
#include "templates/read_numeric.inc"
#undef PROCNAME

#undef FT_TYPE
#undef h5fYPE_EXPR
#undef COMPLEX_FLAG

! ---- real(real64)
#define FT_TYPE        real(real64)
#define h5fYPE_EXPR() hdf5_double()
#define COMPLEX_FLAG   .false.

#define PROCNAME       write_real_real64
#include "templates/write_numeric.inc"
#undef PROCNAME

#define PROCNAME       read_real_real64
#include "templates/read_numeric.inc"
#undef PROCNAME

#undef FT_TYPE
#undef h5fYPE_EXPR
#undef COMPLEX_FLAG

! ---- complex(real32)
#define FT_TYPE        complex(real32)
#define h5fYPE_EXPR() hdf5_float()
#define COMPLEX_FLAG   .true.

#define PROCNAME       write_complex_real32
#include "templates/write_numeric.inc"
#undef PROCNAME

#define PROCNAME       read_complex_real32
#include "templates/read_numeric.inc"
#undef PROCNAME

#undef FT_TYPE
#undef h5fYPE_EXPR
#undef COMPLEX_FLAG

! ---- complex(real64)
#define FT_TYPE        complex(real64)
#define h5fYPE_EXPR() hdf5_double()
#define COMPLEX_FLAG   .true.

#define PROCNAME       write_complex_real64
#include "templates/write_numeric.inc"
#undef PROCNAME

#define PROCNAME       read_complex_real64
#include "templates/read_numeric.inc"
#undef PROCNAME

#undef FT_TYPE
#undef h5fYPE_EXPR
#undef COMPLEX_FLAG

end module solhdf5
