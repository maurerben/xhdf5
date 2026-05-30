! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Benedikt Maurer

!> Utilities for describing and validating HDF5 hyperslab selections.
!>
!> The module provides [`hyperslab_type`](/home/benediktmaurer/Codes/exciting/optixciting/src/xhdf5/utilities/hyperslab.f90)
!> as a compact representation of the file-space and memory-space selections used by
!> the HDF5 read and write wrappers. It also offers validation helpers to ensure that
!> those selections are internally consistent before any HDF5 operation is attempted.
module hyperslab
   use iso_fortran_env, only: int32, real32, real64, output_unit
   use hdf5_globals
   use error_handling
   use mpi_utils


   implicit none


   private
   public :: hyperslab_type, check_hyperslabs


   !> Uninitialized rank
   integer, private, parameter :: uninitialized_rank = -1

   !> Description of one HDF5 hyperslab selection.
   !>
   !> A value of this type stores the selection in both file space and memory space:
   !> extents (`datasize`, `memsize`), offsets (`datastart`, `memstart`), and HDF5
   !> hyperslab traversal parameters (`count`, `stride`, `datablock`).
   !>
   !> Offsets are stored internally in HDF5's 0-based convention. The initializer
   !> accepts Fortran-style 1-based offsets and converts them automatically.
   type hyperslab_type
      !> Rank of the dataset
      integer(int32) :: rank = uninitialized_rank
      !> Shape of the local data chunk handled by the current MPI rank
      integer(hdf5_size), allocatable :: count(:)
      !> Shape of the full dataset (over all MPI ranks)
      integer(hdf5_size), allocatable :: datasize(:)
      !> Shape of the full dataset (over all MPI ranks)
      integer(hdf5_size), allocatable :: memsize(:)
      !> Offset of the local data chunk in context of the dataspace
      integer(hdf5_size), allocatable :: datastart(:)
      !> Offset of the local data chunk in context of the memspace
      integer(hdf5_size), allocatable :: memstart(:)
      !> Array of how many elements to move in each direction
      integer(hdf5_size), allocatable :: stride(:)
      !> Array of element size
      integer(hdf5_size), allocatable :: datablock(:)
      !> Flag to identify hyperslabs for complex datasets.
      logical :: is_complex
   contains
      procedure :: init, delete, check_datasize
   end type


contains


   !> init a hyperslab description from Fortran-facing array metadata.
   !>
   !> Sentinel arrays with negative sums, typically `[-1]`, activate default values:
   !> `datasize` and `memsize` default to `count`, `datastart` and `memstart` default
   !> to zero offsets, and `stride` and `datablock` default to one.
   !>
   !> Scalar data is normalized to rank 1 with `count = [1]`. For complex datasets, an
   !> additional leading dimension of size 2 is inserted so the data can be treated as a
   !> real dataset with one extra rank in the underlying HDF5 representation.
   subroutine init(this, count, datasize, memsize, datastart, memstart, stride, datablock, complex_dataset)
      !> Hyperslab instance to init.
      class(hyperslab_type), intent(inout) :: this
      !> Number of selected elements in each direction.
      integer(int32), intent(in) :: count(:)
      !> Extent of the full dataset in file space.
      !>
      !> Pass `[-1]` to default to `count`.
      integer(int32), intent(in) :: datasize(:)
      !> Extent of the full buffer in memory space.
      !>
      !> Pass `[-1]` to default to `count`.
      integer(int32), intent(in) :: memsize(:)
      !> 1-based offset of the selected chunk in file space.
      !>
      !> Pass `[-1]` to default to the first element in every dimension.
      integer(int32), intent(in) :: datastart(:)
      !> 1-based offset of the selected chunk in memory space.
      !>
      !> Pass `[-1]` to default to the first element in every dimension.
      integer(int32), intent(in) :: memstart(:)
      !> Distance between successive selected blocks in each direction.
      !>
      !> Pass `[-1]` to default to unit stride.
      integer(int32), intent(in) :: stride(:)
      !> Size of each selected block in each direction.
      !>
      !> Pass `[-1]` to default to one element per block.
      integer(int32), intent(in) :: datablock(:)
      !> Whether the dataset stores complex values.
      !>
      !> If true, the hyperslab is expanded with a leading real/imaginary dimension.
      logical, intent(in) :: complex_dataset

      this%count = count
      this%rank = size(this%count)

      ! set rank and shape for scalar input
      if (this%rank == 0) then
         this%count = int([1], hdf5_size)
         this%rank = 1
      end if

      if (sum(datasize) < 0) then
         this%datasize = this%count
      else
         this%datasize = int(datasize, hdf5_size)
      end if

      if (sum(memsize) < 0) then
         this%memsize = this%count
      else
         this%memsize = int(memsize, hdf5_size)
      end if

      if (sum(datastart) < 0) then
         this%datastart = spread(int(0, hdf5_size), 1, this%rank)
      else
         this%datastart = int(datastart - 1, hdf5_size)
      end if

      if (sum(memstart) < 0) then
         this%memstart = spread(int(0, hdf5_size), 1, this%rank)
      else
         this%memstart = int(memstart - 1, hdf5_size)
      end if

      if (sum(stride) < 0) then
         this%stride = spread(int(1, hdf5_size), 1, this%rank)
      else
         this%stride = int(stride, hdf5_size)
      end if

      if (sum(datablock) < 0) then
         this%datablock = spread(int(1, hdf5_size), 1, this%rank)
      else
         this%datablock = int(datablock, hdf5_size)
      end if

      ! complex is treated as real with an extra rank
      this%is_complex = complex_dataset
      if(complex_dataset) then
         this%rank      = this%rank + 1
         this%count     = [[int(2, hdf5_size)], this%count]
         this%datasize  = [[int(2, hdf5_size)], this%datasize]
         this%memsize   = [[int(2, hdf5_size)], this%memsize]
         this%datastart = [[int(0, hdf5_size)], this%datastart]
         this%memstart  = [[int(0, hdf5_size)], this%memstart]
         this%stride    = [[int(1, hdf5_size)], this%stride]
         this%datablock = [[int(1, hdf5_size)], this%datablock]
      end if

   end subroutine init

   !> Validate one or more hyperslab descriptions before an HDF5 operation.
   !>
   !> The routine checks initialization state, rank consistency, positive extents,
   !> non-negative offsets, and that the selected region stays within both file-space
   !> and memory-space bounds. All hyperslabs in the array must share the same rank,
   !> `datasize`, and `memsize`.
   subroutine check_hyperslabs(this, calling_routine, mpi_comm)
      !> Hyperslab array to validate.
      class(hyperslab_type), intent(inout) :: this(:)
      !> Name of the calling routine, used in generated error messages.
      character(*), intent(in) :: calling_routine
      !> MPI communicator used by the assertion helper.
      type(mpi_comm_type), intent(in) :: mpi_comm

      integer :: hsdx
      integer(hdf5_size), allocatable :: file_last_exclusive(:), mem_last_exclusive(:)

      call assert_true(size(this) > 0, &
         'Error(xhdf5: '// calling_routine //'): hyperslab array is empty.')

      do hsdx=1, size(this)

         call assert_true(this(hsdx)%rank /= uninitialized_rank, &
            'Error(xhdf5: '// calling_routine //'): hyperslab is not initialized.')

         call assert_true(allocated(this(hsdx)%count), &
            'Error(xhdf5: '// calling_routine //'): hyperslab count is not allocated.')
         call assert_true(allocated(this(hsdx)%datasize), &
            'Error(xhdf5: '// calling_routine //'): hyperslab datasize is not allocated.')
         call assert_true(allocated(this(hsdx)%memsize), &
            'Error(xhdf5: '// calling_routine //'): hyperslab memsize is not allocated.')
         call assert_true(allocated(this(hsdx)%datastart), &
            'Error(xhdf5: '// calling_routine //'): hyperslab datastart is not allocated.')
         call assert_true(allocated(this(hsdx)%memstart), &
            'Error(xhdf5: '// calling_routine //'): hyperslab memstart is not allocated.')
         call assert_true(allocated(this(hsdx)%stride), &
            'Error(xhdf5: '// calling_routine //'): hyperslab stride is not allocated.')
         call assert_true(allocated(this(hsdx)%datablock), &
            'Error(xhdf5: '// calling_routine //'): hyperslab datablock is not allocated.')

         call assert_true(size(this(hsdx)%count) == this(hsdx)%rank, &
            'Error(xhdf5: '// calling_routine //'): size(count) /= rank.')
         call assert_true(size(this(hsdx)%datasize) == this(hsdx)%rank, &
            'Error(xhdf5: '// calling_routine //'): size(datasize) /= rank.')
         call assert_true(size(this(hsdx)%memsize) == this(hsdx)%rank, &
            'Error(xhdf5: '// calling_routine //'): size(memsize) /= rank.')
         call assert_true(size(this(hsdx)%datastart) == this(hsdx)%rank, &
            'Error(xhdf5: '// calling_routine //'): size(datastart) /= rank.')
         call assert_true(size(this(hsdx)%memstart) == this(hsdx)%rank, &
            'Error(xhdf5: '// calling_routine //'): size(memstart) /= rank.')
         call assert_true(size(this(hsdx)%stride) == this(hsdx)%rank, &
            'Error(xhdf5: '// calling_routine //'): size(stride) /= rank.')
         call assert_true(size(this(hsdx)%datablock) == this(hsdx)%rank, &
            'Error(xhdf5: '// calling_routine //'): size(datablock) /= rank.')

         call assert_true(this(hsdx)%rank == this(1)%rank, &
            'Error(xhdf5: '// calling_routine //'): hyperslabs have inconsistent rank.')
         call assert_true(all(this(hsdx)%datasize == this(1)%datasize), &
            'Error(xhdf5: '// calling_routine //'): hyperslabs have inconsistent datasize.')
         call assert_true(all(this(hsdx)%memsize == this(1)%memsize), &
            'Error(xhdf5: '// calling_routine //'): hyperslabs have inconsistent memsize.')

         call assert_true(all(this(hsdx)%datasize > 0), &
            'Error(xhdf5: '// calling_routine //'): Some elements of datasize <= 0.')
         call assert_true(all(this(hsdx)%memsize >= 0), &
            'Error(xhdf5: '// calling_routine //'): Some elements of memsize < 0.')
         call assert_true(all(this(hsdx)%count >= 0), &
            'Error(xhdf5: '// calling_routine //'): Some elements of count < 0.')
         call assert_true(all(this(hsdx)%stride > 0), &
            'Error(xhdf5: '// calling_routine //'): Some elements of stride <= 0.')
         call assert_true(all(this(hsdx)%datablock > 0), &
            'Error(xhdf5: '// calling_routine //'): Some elements of datablock <= 0.')
         call assert_true(all(this(hsdx)%datastart >= 0), &
            'Error(xhdf5: '// calling_routine //'): Some elements of datastart < 0.')
         call assert_true(all(this(hsdx)%memstart >= 0), &
            'Error(xhdf5: '// calling_routine //'): Some elements of memstart < 0.')

         ! Hyperslab bound: start + (count - 1) * stride + datablock <= extent.
         if(all(this(hsdx)%count > 0)) then
            file_last_exclusive = this(hsdx)%datastart + (this(hsdx)%count - 1) * this(hsdx)%stride + this(hsdx)%datablock
            mem_last_exclusive = this(hsdx)%memstart + (this(hsdx)%count - 1) * this(hsdx)%stride + this(hsdx)%datablock

            call assert_true(all(file_last_exclusive <= this(hsdx)%datasize), &
               'Error(xhdf5: '// calling_routine //'): hyperslab exceeds datasize.')
            call assert_true(all(mem_last_exclusive <= this(hsdx)%memsize), &
               'Error(xhdf5: '// calling_routine //'): hyperslab exceeds memsize.')
         end if
      end do
   end subroutine check_hyperslabs

   !> Reset a hyperslab instance and release its allocated metadata arrays.
   subroutine delete(this)
      !> Hyperslab array to validate.
      class(hyperslab_type), intent(inout) :: this

      this%rank = uninitialized_rank
      if(allocated(this%count)) deallocate(this%count)
      if(allocated(this%datasize)) deallocate(this%datasize)
      if(allocated(this%datastart)) deallocate(this%datastart)
   end subroutine delete

   !> Given an integer array representing the full shape of an array, return
   !> `.true.` if  the attribute `[[datasize]]` matches
   logical function check_datasize(this, size_to_check)
      !> Hyperslab array to validate.
      class(hyperslab_type), intent(in) :: this
      !> Datasize to check
      integer(int32), intent(in) :: size_to_check(:)


      if(this%is_complex) then
         if(size(size_to_check) /= this%rank-1) then
            check_datasize = .false.
         else
            check_datasize = all(size_to_check == int(this%datasize(2:), hdf5_size))
         endif
      else
         if(size(size_to_check) /= this%rank) then
            check_datasize = .false.
         else
            check_datasize = all(size_to_check == int(this%datasize, hdf5_size))
         endif
      endif
   end function


end module hyperslab
