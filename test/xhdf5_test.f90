module solhdf5_test
   use iso_fortran_env, only: real32, real64
   use modmpi, only: mpiinfo, barrier, nofset, firstofset, lastofset
   use unit_test_framework, only : unit_test_type
   use solhdf5
   use os_utils, only: system_cmd
   use math_utils, only: mod1, all_close

   implicit none

   private
   public :: solhdf5_test_driver

contains

   !> Run tests for math tools
   subroutine solhdf5_test_driver(mpiglobal, kill_on_failure)
      !> mpi information
      type(mpiinfo), intent(in) :: mpiglobal
      !> Kill the program before the test driver finishes
      !> if an assertion fails
      logical, optional :: kill_on_failure
      !> test object
      type(unit_test_type) :: test_report
      !> Number of assertions
      integer, parameter :: n_assertions = 47

#ifdef _HDF5_
      ! init test object
      call test_report%init(mpiglobal)

      ! Run and assert tests

      call test_solhdf5_file(test_report, mpiglobal)
      call barrier(mpiglobal)

      ! CHARACTER

      call test_solhdf5_write_and_read_character(test_report, mpiglobal)
      call barrier(mpiglobal)

      ! INTEGER

      call test_solhdf5_write_and_read_integer_rank_0(test_report, mpiglobal)
      call test_solhdf5_write_and_read_integer_rank_1(test_report, mpiglobal)
      call test_solhdf5_write_and_read_integer_rank_2(test_report, mpiglobal)
      call test_solhdf5_write_and_read_integer_rank_3(test_report, mpiglobal)

      ! REAL(real32)

      call test_solhdf5_write_and_read_real_real32_rank_0(test_report, mpiglobal)
      call test_solhdf5_write_and_read_real_real32_rank_1(test_report, mpiglobal)
      call barrier(mpiglobal)

      call test_solhdf5_write_and_read_real_real32_rank_2(test_report, mpiglobal)
      call barrier(mpiglobal)

      call test_solhdf5_write_and_read_real_real32_rank_3(test_report, mpiglobal)
      call barrier(mpiglobal)

      call test_solhdf5_write_and_read_real_real32_rank_4(test_report, mpiglobal)
      call barrier(mpiglobal)

      ! REAL(real64)

      call test_solhdf5_write_and_read_real_real64_rank_0(test_report, mpiglobal)
      call barrier(mpiglobal)

      call test_solhdf5_write_and_read_real_real64_rank_1(test_report, mpiglobal)
      call barrier(mpiglobal)

      call test_solhdf5_write_and_read_real_real64_rank_2(test_report, mpiglobal)
      call barrier(mpiglobal)

      call test_solhdf5_write_and_read_real_real64_rank_3(test_report, mpiglobal)
      call barrier(mpiglobal)

      call test_solhdf5_write_and_read_real_real64_rank_4(test_report, mpiglobal)
      call barrier(mpiglobal)

      ! COMPLEX(real64)

      call test_solhdf5_write_and_read_cmplx_real64_rank_1(test_report, mpiglobal)
      call barrier(mpiglobal)

      call test_solhdf5_write_and_read_cmplx_real64_rank_2(test_report, mpiglobal)
      call barrier(mpiglobal)

      call test_solhdf5_write_and_read_cmplx_real64_rank_3(test_report, mpiglobal)
      call barrier(mpiglobal)

      call test_solhdf5_write_and_read_cmplx_real64_rank_4(test_report, mpiglobal)
      call barrier(mpiglobal)
#else
      print*, 'Built without HDF5. Nothing to do here.'
      call test_report%init(0, mpiglobal)
#endif

      ! report results
      if (present(kill_on_failure)) then
         call test_report%report('solhdf5', kill_on_failure)
      else
         call test_report%report('solhdf5')
      end if

      ! Finalise test object
      call test_report%finalise()
   end subroutine solhdf5_test_driver


   !> Test hdf5 file utilities
   subroutine test_solhdf5_file(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr

      real(real64), parameter :: real_vec(5) = [1._real64, -23.0123_real64, 0._real64, 123.3245_real64, 3142._real64]
      complex(real64), parameter :: cmplx_vec(5) = cmplx(real_vec, -real_vec, real64)
      integer, allocatable :: dset_shape(:), dset_shape_ref(:)

      ! Create a new HDF5 file
      call h5f%init('test_file.h5', mpiglobal)
      call h5f%finalize()

      ! Open an existing HDF5 file
      call h5f%init('test_file.h5', mpiglobal)

      ! Create a new group
      call h5f%initialize_group('./', 'new_group')
      call h5f%initialize_group('./', 'new_group') ! Can initialize_group handle an existing group?

      ! Create a sub group
      call h5f%initialize_group('new_group', 'sub_group')
      call h5f%initialize_group('new_group', 'sub_group') ! Can initialize_group handle an existing group?

      ! Test link_exists
      call test_report%assert(h5f%exists('new_group'), &
         "Failed to create 'new_group'.")

      call h5f%delete('new_group')
      call test_report%assert(.not. h5f%exists('new_group'), &
         "Failed to delete 'new_group'.")

      call h5f%write('/', 'real_vec', real_vec, [1], [5])
      call h5f%write('/', 'cmplx_vec', cmplx_vec, [1], [5])

      call h5f%dataset_shape('/', 'real_vec', dset_shape)
      call test_report%assert(all(dset_shape == [5]), &
         'dataset shape for real vector not correctly detected.')

      deallocate(dset_shape)
      call h5f%dataset_shape('/', 'cmplx_vec', dset_shape)
      call test_report%assert(all(dset_shape == [2, 5]), &
         'dataset shape for complex vector not correctly detected (without complex specification).')

      deallocate(dset_shape)
      call h5f%dataset_shape('/', 'cmplx_vec', dset_shape, complex_dataset=.true.)
      call test_report%assert(all(dset_shape == [5]), &
         'dataset shape for complex vector not correctly detected (with complex specification).')

      deallocate(dset_shape)
      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_file


! CHARACTER


   !> Test read and write character.
   subroutine test_solhdf5_write_and_read_character(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr

      character(:), allocatable :: char_write, char_read
      logical :: bool_write, bool_read

      char_write = 'aloifhSADpi$%^flmdsa'

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')

      call h5f%write('datasets', 'character', char_write)
      call h5f%read('datasets', 'character', char_read)
      call test_report%assert(char_write .eq. char_read, &
         'character: Read array differs from written array.')

      bool_write = .true.
      call h5f%write('datasets', 'bool_true', bool_write)
      call h5f%read('datasets', 'bool_true', bool_read)
      call test_report%assert(bool_read .eqv. bool_write, &
         'bool: Read and write .true. failed.')

      bool_write = .false.
      call h5f%write('datasets', 'bool_true', bool_write)
      call h5f%read('datasets', 'bool_true', bool_read)
      call test_report%assert(bool_read .eqv. bool_write, &
         'bool: Read and write .false. failed.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_character


! INTEGER


   !> Test read and write for integer scalars.
   subroutine test_solhdf5_write_and_read_integer_rank_0(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr

      integer :: integer_write, integer_read

      integer_write = 15

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')

      call h5f%write('datasets', 'integer_r1', integer_write)
      call h5f%read('datasets', 'integer_r1', integer_read)
      call test_report%assert(integer_write == integer_read, &
         'integer rank 0: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_integer_rank_0


   !> Test read and write for integer rank 1 arrays.
   subroutine test_solhdf5_write_and_read_integer_rank_1(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, l

      integer :: integer_r1(8)
      integer, allocatable :: datachunk(:), datachunk_read(:)
      integer :: offset(1), dataset_shape(1), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(integer_r1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(integer_r1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(integer_r1), mpiglobal%procs)

      do l = 1, size(integer_r1)

         integer_r1(l) =  l

      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      datachunk = integer_r1(first : last)
      datachunk_read = spread(0, 1, size(datachunk))
      offset = first
      dataset_shape = shape(integer_r1)

      call h5f%write('datasets', 'integer_r1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'integer_r1', datachunk_read, offset)
      call test_report%assert(all(datachunk == datachunk_read), &
         'integer rank 1: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_integer_rank_1


   !> Test read and write for integer rank 2 arrays.
   subroutine test_solhdf5_write_and_read_integer_rank_2(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      integer :: integer_r2(8, 3)
      integer, allocatable :: datachunk(:, :), datachunk_read(:, :)
      integer :: offset(2), dataset_shape(2), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(integer_r2, 2), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(integer_r2, 2), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(integer_r2, 2), mpiglobal%procs)

      do l = 1, size(integer_r2, 2)
         do k = 1, size(integer_r2, 1)

            integer_r2(k, l) =  (l - k)

         end do
      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      ! Distribute writing over second rank
      datachunk = integer_r2(:, first : last)
      datachunk_read = datachunk
      datachunk_read = 0.
      offset = [1, first]
      dataset_shape = shape(integer_r2)

      call h5f%write('datasets', 'integer_r2-2', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'integer_r2-2', datachunk_read, offset)
      call test_report%assert(all(datachunk == datachunk_read), &
         'integer rank 2, dim 2 distributed: Read array differs from written array.')


      ! Distribute writing over first rank
      chunk_size = nofset(mpiglobal%rank, size(integer_r2, 1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(integer_r2, 1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(integer_r2, 1), mpiglobal%procs)

      datachunk = integer_r2(first : last, :)
      datachunk_read = datachunk
      datachunk_read = 0.
      offset = [first, 1]
      dataset_shape = shape(integer_r2)

      call h5f%write('datasets', 'integer_r2-1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'integer_r2-1', datachunk_read, offset)
      call test_report%assert(all(datachunk == datachunk_read), &
         'integer rank 2, dim 1 distributed: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_integer_rank_2


   !> Test read and write for integer rank 3 arrays.
   subroutine test_solhdf5_write_and_read_integer_rank_3(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      real(real32) :: integer_r3(8, 3, 5)
      real(real32), allocatable :: datachunk(:, :, :), datachunk_read(:, :, :)
      integer :: offset(3), dataset_shape(3), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(integer_r3, 3), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(integer_r3, 3), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(integer_r3, 3), mpiglobal%procs)

      do l = 1, size(integer_r3, 3)
         do k = 1, size(integer_r3, 2)
            do j = 1, size(integer_r3, 1)

               integer_r3(j, k, l) =  (j - k) + l

            end do
         end do
      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      ! Distribute writing over third rank
      datachunk = integer_r3(:, :, first : last)
      datachunk_read = datachunk
      datachunk_read = 0
      offset = [1, 1, first]
      dataset_shape = shape(integer_r3)

      call h5f%write('datasets', 'integer_r3-3', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'integer_r3-3', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 3, dim 3 distributed: Read array differs from written array.')


      ! Distribute writing over second rank
      chunk_size = nofset(mpiglobal%rank, size(integer_r3, 2), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(integer_r3, 2), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(integer_r3, 2), mpiglobal%procs)

      datachunk = integer_r3(:, first : last, :)
      datachunk_read = datachunk
      datachunk_read = 0
      offset = [1, first, 1]
      dataset_shape = shape(integer_r3)

      call h5f%write('datasets', 'integer_r3-2', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'integer_r3-2', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 3, dim 2 distributed: Read array differs from written array.')


      ! Distribute writing over first rank
      chunk_size = nofset(mpiglobal%rank, size(integer_r3, 1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(integer_r3, 1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(integer_r3, 1), mpiglobal%procs)

      datachunk = integer_r3(first : last, :, :)
      datachunk_read = datachunk
      datachunk_read = 0
      offset = [first, 1, 1]
      dataset_shape = shape(integer_r3)

      call h5f%write('datasets', 'integer_r3-1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'integer_r3-1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 3, dim 1 distributed: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_integer_rank_3


   ! REAL(real32)


   !> Test read and write for real(real32) scalars.
   subroutine test_solhdf5_write_and_read_real_real32_rank_0(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      real(real32) :: real_write, real_read

      real_write = 0.213_real32

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')

      call h5f%write('datasets', 'real_r1', real_write)
      call h5f%read('datasets', 'real_r1', real_read)
      call test_report%assert(all_close(real_write, real_read), &
         'real(real32) rank 0: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_real_real32_rank_0


   !> Test read and write for real(real32) rank 1 arrays.
   subroutine test_solhdf5_write_and_read_real_real32_rank_1(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      real(real32) :: real_r1(8)
      real(real32), allocatable :: datachunk(:), datachunk_read(:)
      integer :: offset(1), dataset_shape(1), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)

      do l = 1, size(real_r1)

         real_r1(l) =  real(l, real32)

      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      datachunk = real_r1(first : last)
      datachunk_read = datachunk
      datachunk_read = 0._real32
      offset = first
      dataset_shape = shape(real_r1)

      call h5f%write('datasets', 'real_r1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 1: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_real_real32_rank_1


   !> Test read and write for real(real32) rank 2 arrays.
   subroutine test_solhdf5_write_and_read_real_real32_rank_2(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      real(real32) :: real_r2(8, 3)
      real(real32), allocatable :: datachunk(:, :), datachunk_read(:, :)
      integer :: offset(2), dataset_shape(2), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)

      do l = 1, size(real_r2, 2)
         do k = 1, size(real_r2, 1)

            real_r2(k, l) =  real((l - k), real32)

         end do
      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      ! Distribute writing over second rank
      datachunk = real_r2(:, first : last)
      datachunk_read = datachunk
      datachunk_read = 0._real32
      offset = [1, first]
      dataset_shape = shape(real_r2)

      call h5f%write('datasets', 'real_r2-2', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r2-2', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 2, dim 2 distributed: Read array differs from written array.')


      ! Distribute writing over first rank
      chunk_size = nofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)

      datachunk = real_r2(first : last, :)
      datachunk_read = datachunk
      datachunk_read = 0._real32
      offset = [first, 1]
      dataset_shape = shape(real_r2)

      call h5f%write('datasets', 'real_r2-1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r2-1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 2, dim 1 distributed: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_real_real32_rank_2


   !> Test read and write for real(real32) rank 3 arrays.
   subroutine test_solhdf5_write_and_read_real_real32_rank_3(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      real(real32) :: real_r3(8, 3, 5)
      real(real32), allocatable :: datachunk(:, :, :), datachunk_read(:, :, :)
      integer :: offset(3), dataset_shape(3), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)

      do l = 1, size(real_r3, 3)
         do k = 1, size(real_r3, 2)
            do j = 1, size(real_r3, 1)

               real_r3(j, k, l) =  real((j - k), real32) / l

            end do
         end do
      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      ! Distribute writing over third rank
      datachunk = real_r3(:, :, first : last)
      datachunk_read = datachunk
      datachunk_read = 0._real32
      offset = [1, 1, first]
      dataset_shape = shape(real_r3)

      call h5f%write('datasets', 'real_r3-3', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r3-3', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 3, dim 3 distributed: Read array differs from written array.')


      ! Distribute writing over second rank
      chunk_size = nofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)

      datachunk = real_r3(:, first : last, :)
      datachunk_read = datachunk
      datachunk_read = 0._real32
      offset = [1, first, 1]
      dataset_shape = shape(real_r3)

      call h5f%write('datasets', 'real_r3-2', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r3-2', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 3, dim 2 distributed: Read array differs from written array.')


      ! Distribute writing over first rank
      chunk_size = nofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)

      datachunk = real_r3(first : last, :, :)
      datachunk_read = datachunk
      datachunk_read = 0._real32
      offset = [first, 1, 1]
      dataset_shape = shape(real_r3)

      call h5f%write('datasets', 'real_r3-1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r3-1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 3, dim 1 distributed: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_real_real32_rank_3


   !> Test read and write for real(real32) rank 4 arrays.
   subroutine test_solhdf5_write_and_read_real_real32_rank_4(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l, m

      real(real32) :: real_r4(8, 3, 5, 4)
      real(real32), allocatable :: datachunk(:, :, :, :), datachunk_read(:, :, :, :)
      integer :: offset(4), dataset_shape(4), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)

      do m = 1, size(real_r4, 4)
         do l = 1, size(real_r4, 3)
            do k = 1, size(real_r4, 2)
               do j = 1, size(real_r4, 1)

                  real_r4(j, k, l, m) = real((j - k) * (l - m), real32)
               end do
            end do
         end do
      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      ! Distribute writing over fourth rank
      datachunk = real_r4(:, :, :, first : last)
      datachunk_read = datachunk
      datachunk_read = 0._real32
      offset = [1, 1, 1, first]
      dataset_shape = shape(real_r4)

      call h5f%write('datasets', 'real_r4-4', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r4-4', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 4, dim 4 distributed: Read array differs from written array.')


      ! Distribute writing over third rank
      chunk_size = nofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)

      datachunk = real_r4(:, :, first : last, :)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real32)
      offset = [1, 1, first, 1]
      dataset_shape = shape(real_r4)

      call h5f%write('datasets', 'real_r4-3', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r4-3', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 4, dim 3 distributed: Read array differs from written array.')

      ! Distribute writing over second rank
      chunk_size = nofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)

      datachunk = real_r4(:, first : last, :, :)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real32)
      offset = [1, first, 1, 1]
      dataset_shape = shape(real_r4)

      call h5f%write('datasets', 'real_r4-2', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r4-2', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 4, dim 2 distributed: Read array differs from written array.')


      ! Distribute writing over first rank
      chunk_size = nofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)

      datachunk = real_r4(first : last, :, :, :)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real32)
      offset = [first, 1, 1, 1]
      dataset_shape = shape(real_r4)

      call h5f%write('datasets', 'real_r4-1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r4-1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real32) rank 4, dim 1 distributed: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_real_real32_rank_4


   ! REAL(real64)


   !> Test read and write for real(real64) scalars.
   subroutine test_solhdf5_write_and_read_real_real64_rank_0(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      real(real64) :: real_write, real_read

      real_write = 0.213_real64

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')

      call h5f%write('datasets', 'real_r1', real_write)
      call h5f%read('datasets', 'real_r1', real_read)
      call test_report%assert(all_close(real_write, real_read), &
         'real(real64) rank 0: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_real_real64_rank_0


   !> Test read and write for real(real64) rank 1 arrays.
   subroutine test_solhdf5_write_and_read_real_real64_rank_1(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      real(real64) :: real_r1(8)
      real(real64), allocatable :: datachunk(:), datachunk_read(:)
      integer :: offset(1), dataset_shape(1), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r1), mpiglobal%procs)

      do l = 1, size(real_r1)

         real_r1(l) =  real(l, real64)

      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      datachunk = real_r1(first : last)
      datachunk_read = datachunk
      datachunk_read = 0._real64
      offset = first
      dataset_shape = shape(real_r1)

      call h5f%write('datasets', 'real_r1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 1: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_real_real64_rank_1


   !> Test read and write for real(real64) rank 2 arrays.
   subroutine test_solhdf5_write_and_read_real_real64_rank_2(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      real(real64) :: real_r2(8, 3)
      real(real64), allocatable :: datachunk(:, :), datachunk_read(:, :)
      integer :: offset(2), dataset_shape(2), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r2, 2), mpiglobal%procs)

      do l = 1, size(real_r2, 2)
         do k = 1, size(real_r2, 1)

            real_r2(k, l) =  real((l - k), real64)

         end do
      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      ! Distribute writing over second rank
      datachunk = real_r2(:, first : last)
      datachunk_read = datachunk
      datachunk_read = 0._real64
      offset = [1, first]
      dataset_shape = shape(real_r2)

      call h5f%write('datasets', 'real_r2-2', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r2-2', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 2, dim 2 distributed: Read array differs from written array.')


      ! Distribute writing over first rank
      chunk_size = nofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r2, 1), mpiglobal%procs)

      datachunk = real_r2(first : last, :)
      datachunk_read = datachunk
      datachunk_read = 0._real64
      offset = [first, 1]
      dataset_shape = shape(real_r2)

      call h5f%write('datasets', 'real_r2-1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r2-1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 2, dim 1 distributed: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_real_real64_rank_2


   !> Test read and write for real(real64) rank 3 arrays.
   subroutine test_solhdf5_write_and_read_real_real64_rank_3(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      real(real64) :: real_r3(8, 3, 5)
      real(real64), allocatable :: datachunk(:, :, :), datachunk_read(:, :, :)
      integer :: offset(3), dataset_shape(3), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r3, 3), mpiglobal%procs)

      do l = 1, size(real_r3, 3)
         do k = 1, size(real_r3, 2)
            do j = 1, size(real_r3, 1)

               real_r3(j, k, l) =  real((j - k), real64) / l

            end do
         end do
      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      ! Distribute writing over third rank
      datachunk = real_r3(:, :, first : last)
      datachunk_read = datachunk
      datachunk_read = 0._real64
      offset = [1, 1, first]
      dataset_shape = shape(real_r3)

      call h5f%write('datasets', 'real_r3-3', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r3-3', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 3, dim 3 distributed: Read array differs from written array.')


      ! Distribute writing over second rank
      chunk_size = nofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r3, 2), mpiglobal%procs)

      datachunk = real_r3(:, first : last, :)
      datachunk_read = datachunk
      datachunk_read = 0._real64
      offset = [1, first, 1]
      dataset_shape = shape(real_r3)

      call h5f%write('datasets', 'real_r3-2', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r3-2', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 3, dim 2 distributed: Read array differs from written array.')


      ! Distribute writing over first rank
      chunk_size = nofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r3, 1), mpiglobal%procs)

      datachunk = real_r3(first : last, :, :)
      datachunk_read = datachunk
      datachunk_read = 0._real64
      offset = [first, 1, 1]
      dataset_shape = shape(real_r3)

      call h5f%write('datasets', 'real_r3-1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r3-1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 3, dim 1 distributed: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_real_real64_rank_3


   !> Test read and write for real(real64) rank 4 arrays.
   subroutine test_solhdf5_write_and_read_real_real64_rank_4(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l, m

      real(real64) :: real_r4(8, 3, 5, 4)
      real(real64), allocatable :: datachunk(:, :, :, :), datachunk_read(:, :, :, :)
      integer :: offset(4), dataset_shape(4), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r4, 4), mpiglobal%procs)

      do m = 1, size(real_r4, 4)
         do l = 1, size(real_r4, 3)
            do k = 1, size(real_r4, 2)
               do j = 1, size(real_r4, 1)

                  real_r4(j, k, l, m) = real((j - k) * (l - m), real64)
               end do
            end do
         end do
      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      ! Distribute writing over fourth rank
      datachunk = real_r4(:, :, :, first : last)
      datachunk_read = datachunk
      datachunk_read = 0._real64
      offset = [1, 1, 1, first]
      dataset_shape = shape(real_r4)

      call h5f%write('datasets', 'real_r4-4', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r4-4', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 4, dim 4 distributed: Read array differs from written array.')


      ! Distribute writing over third rank
      chunk_size = nofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r4, 3), mpiglobal%procs)

      datachunk = real_r4(:, :, first : last, :)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real64)
      offset = [1, 1, first, 1]
      dataset_shape = shape(real_r4)

      call h5f%write('datasets', 'real_r4-3', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r4-3', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 4, dim 3 distributed: Read array differs from written array.')

      ! Distribute writing over second rank
      chunk_size = nofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r4, 2), mpiglobal%procs)

      datachunk = real_r4(:, first : last, :, :)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real64)
      offset = [1, first, 1, 1]
      dataset_shape = shape(real_r4)

      call h5f%write('datasets', 'real_r4-2', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r4-2', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 4, dim 2 distributed: Read array differs from written array.')


      ! Distribute writing over first rank
      chunk_size = nofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(real_r4, 1), mpiglobal%procs)

      datachunk = real_r4(first : last, :, :, :)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real64)
      offset = [first, 1, 1, 1]
      dataset_shape = shape(real_r4)

      call h5f%write('datasets', 'real_r4-1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'real_r4-1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 4, dim 1 distributed: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_real_real64_rank_4


   ! COMPLEX(real64)


   !> Test read and write for complex(real64) rank 1 arrays.
   subroutine test_solhdf5_write_and_read_cmplx_real64_rank_1(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      complex(real64) :: cmplx_r1(8)
      complex(real64), allocatable :: datachunk(:), datachunk_read(:)
      integer :: offset(1), dataset_shape(1), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(cmplx_r1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(cmplx_r1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(cmplx_r1), mpiglobal%procs)

      do l = 1, size(cmplx_r1)

         cmplx_r1(l) =  cmplx(l, -l, real64)

      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')

      datachunk = cmplx_r1(first : last)
      datachunk_read = datachunk
      datachunk_read = 0._real64
      offset = first
      dataset_shape = shape(cmplx_r1)

      call h5f%write('datasets', 'cmplx_r1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'cmplx_r1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'complex rank 1: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_cmplx_real64_rank_1


   !> Test read and write for complex(real64) rank 2 arrays.
   subroutine test_solhdf5_write_and_read_cmplx_real64_rank_2(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      complex(real64) :: cmplx_r2(8, 3)
      complex(real64), allocatable :: datachunk(:, :), datachunk_read(:, :)
      integer :: offset(2), dataset_shape(2), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(cmplx_r2, 2), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(cmplx_r2, 2), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(cmplx_r2, 2), mpiglobal%procs)

      do l = 1, size(cmplx_r2, 2)
         do k = 1, size(cmplx_r2, 1)

            cmplx_r2(k, l) =  cmplx((l + k), (l - k), real64)

         end do
      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      ! Distribute writing over second rank
      datachunk = cmplx_r2(:, first : last)
      datachunk_read = datachunk
      datachunk_read = 0._real64
      offset = [1, first]
      dataset_shape = shape(cmplx_r2)

      call h5f%write('datasets', 'cmplx_r2-2', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'cmplx_r2-2', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 2, dim 2 distributed: Read array differs from written array.')


      ! Distribute writing over first rank
      chunk_size = nofset(mpiglobal%rank, size(cmplx_r2, 1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(cmplx_r2, 1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(cmplx_r2, 1), mpiglobal%procs)

      datachunk = cmplx_r2(first : last, :)
      datachunk_read = datachunk
      datachunk_read = 0._real64
      offset = [first, 1]
      dataset_shape = shape(cmplx_r2)

      call h5f%write('datasets', 'cmplx_r2-1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'cmplx_r2-1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'real(real64) rank 2, dim 1 distributed: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_cmplx_real64_rank_2


   !> Test read and write for complex(real64) rank 3 arrays.
   subroutine test_solhdf5_write_and_read_cmplx_real64_rank_3(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l

      complex(real64) :: cmplx_r3(8, 3, 5)
      complex(real64), allocatable :: datachunk(:, :, :), datachunk_read(:, :, :)
      integer :: offset(3), dataset_shape(3), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(cmplx_r3, 3), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(cmplx_r3, 3), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(cmplx_r3, 3), mpiglobal%procs)

      do l = 1, size(cmplx_r3, 3)
         do k = 1, size(cmplx_r3, 2)
            do j = 1, size(cmplx_r3, 1)

               cmplx_r3(j, k, l) = cmplx((j - k) * l, (k - j) * l, real64)

            end do
         end do
      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      ! Distribute writing over third rank
      datachunk = cmplx_r3(:, :, first : last)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real64)
      offset = [1, 1, first]
      dataset_shape = shape(cmplx_r3)

      call h5f%write('datasets', 'cmplx_r3-3', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'cmplx_r3-3', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'complex(real64) rank 3, dim 3 distributed: Read array differs from written array.')


      ! Distribute writing over second rank
      chunk_size = nofset(mpiglobal%rank, size(cmplx_r3, 2), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(cmplx_r3, 2), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(cmplx_r3, 2), mpiglobal%procs)

      datachunk = cmplx_r3(:, first : last, :)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real64)
      offset = [1, first, 1]
      dataset_shape = shape(cmplx_r3)

      call h5f%write('datasets', 'cmplx_r3-2', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'cmplx_r3-2', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'complex(real64) rank 3, dim 2 distributed: Read array differs from written array.')


      ! Distribute writing over first rank
      chunk_size = nofset(mpiglobal%rank, size(cmplx_r3, 1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(cmplx_r3, 1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(cmplx_r3, 1), mpiglobal%procs)

      datachunk = cmplx_r3(first : last, :, :)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real64)
      offset = [first, 1, 1]
      dataset_shape = shape(cmplx_r3)

      call h5f%write('datasets', 'cmplx_r3-1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'cmplx_r3-1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'complex(real64) rank 3, dim 1 distributed: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine


   !> Test read and write for complex(real64) rank 4 arrays.
   subroutine test_solhdf5_write_and_read_cmplx_real64_rank_4(test_report, mpiglobal)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(in) :: mpiglobal

      type(h5file_t) :: h5f
      integer :: ierr, j, k ,l, m

      complex(real64) :: cmplx_r4(8, 3, 5, 4)
      complex(real64), allocatable :: datachunk(:, :, :, :), datachunk_read(:, :, :, :)
      integer :: offset(4), dataset_shape(4), chunk_size, first, last

      chunk_size = nofset(mpiglobal%rank, size(cmplx_r4, 4), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(cmplx_r4, 4), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(cmplx_r4, 4), mpiglobal%procs)

      do m = 1, size(cmplx_r4, 4)
         do l = 1, size(cmplx_r4, 3)
            do k = 1, size(cmplx_r4, 2)
               do j = 1, size(cmplx_r4, 1)

                  cmplx_r4(j, k, l, m) = cmplx((j - k) * (l - m), (k - j) * (m - l), real64)

               end do
            end do
         end do
      end do

      call h5f%init('test_file.h5', mpiglobal)
      call h5f%initialize_group('./', 'datasets')


      ! Distribute writing over fourth rank
      datachunk = cmplx_r4(:, :, :, first : last)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real64)
      offset = [1, 1, 1, first]
      dataset_shape = shape(cmplx_r4)

      call h5f%write('datasets', 'cmplx_r4-4', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'cmplx_r4-4', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'complex(real64) rank 4, dim 4 distributed: Read array differs from written array.')


      ! Distribute writing over third rank
      chunk_size = nofset(mpiglobal%rank, size(cmplx_r4, 3), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(cmplx_r4, 3), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(cmplx_r4, 3), mpiglobal%procs)

      datachunk = cmplx_r4(:, :, first : last, :)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real64)
      offset = [1, 1, first, 1]
      dataset_shape = shape(cmplx_r4)

      call h5f%write('datasets', 'cmplx_r4-3', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'cmplx_r4-3', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'complex(real64) rank 4, dim 3 distributed: Read array differs from written array.')


      ! Distribute writing over second rank
      chunk_size = nofset(mpiglobal%rank, size(cmplx_r4, 2), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(cmplx_r4, 2), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(cmplx_r4, 2), mpiglobal%procs)

      datachunk = cmplx_r4(:, first : last, :, :)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real64)
      offset = [1, first, 1, 1]
      dataset_shape = shape(cmplx_r4)

      call h5f%write('datasets', 'cmplx_r4-2', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'cmplx_r4-2', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'complex(real64) rank 4, dim 2 distributed: Read array differs from written array.')


      ! Distribute writing over first rank
      chunk_size = nofset(mpiglobal%rank, size(cmplx_r4, 1), mpiglobal%procs)
      first = firstofset(mpiglobal%rank, size(cmplx_r4, 1), mpiglobal%procs)
      last = lastofset(mpiglobal%rank, size(cmplx_r4, 1), mpiglobal%procs)

      datachunk = cmplx_r4(first : last, :, :, :)
      datachunk_read = datachunk
      datachunk_read = cmplx(0.0, 0.0, real64)
      offset = [first, 1, 1, 1]
      dataset_shape = shape(cmplx_r4)

      call h5f%write('datasets', 'cmplx_r4-1', datachunk, offset, dataset_shape)
      call h5f%read('datasets', 'cmplx_r4-1', datachunk_read, offset)
      call test_report%assert(all_close(datachunk, datachunk_read), &
         'complex(real64) rank 4, dim 1 distributed: Read array differs from written array.')

      call h5f%finalize()
      if(mpiglobal%is_root) ierr = system_cmd('rm test_file.h5')
   end subroutine test_solhdf5_write_and_read_cmplx_real64_rank_4

end module
