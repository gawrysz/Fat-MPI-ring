module ring

   use constants, only: FP64, INT64

   implicit none

   private
   public :: ring_t

   type :: ring_t
      real(kind=FP64), allocatable, dimension(:) :: sbuf
      real(kind=FP64), allocatable, dimension(:) :: rbuf
      integer, private :: test_type
      integer(kind=INT64), private :: n
      logical :: give_up
   contains
      procedure :: init
      procedure, private :: setup
      procedure :: run
      procedure, private :: check
      procedure :: cleanup
   end type ring_t

contains

   subroutine init(this, n, test_type)

      use composition, only: factorization_t
      use memory,      only: memcheck

      implicit none

      class(ring_t),         intent(inout) :: this       ! object invoking type-bound procedure
      type(factorization_t), intent(in)    :: n          ! number of doubles in the buffer
      integer,               intent(in)    :: test_type  ! test to perform

      this%give_up = .false.
      this%n = n%number
      this%test_type = test_type

      call this%cleanup
      allocate(this%rbuf(this%n), &
           &   this%sbuf(this%n))

      if (.not. memcheck()) call exit(-13)

   end subroutine init

   subroutine setup(this)

      use mpisetup, only: proc

      implicit none

      class(ring_t), intent(inout) :: this  ! object invoking type-bound procedure

      this%sbuf = real(proc)
      this%rbuf = huge(1._FP64)

   end subroutine setup

   subroutine check(this)

      use mpisetup, only: proc, nproc

      implicit none

      class(ring_t), intent(in) :: this  ! object invoking type-bound procedure

      if (any(this%rbuf /= real(mod(proc + 1, nproc)))) &
           write(*,*)"### ", count(this%rbuf /= real(mod(proc + 1, nproc))), "/", size(this%rbuf), " wrong values @", proc

   end subroutine check

   subroutine cleanup(this)

      implicit none

      class(ring_t), intent(inout) :: this  ! object invoking type-bound procedure

      if (allocated(this%rbuf)) deallocate(this%rbuf)
      if (allocated(this%sbuf)) deallocate(this%sbuf)

   end subroutine cleanup

   subroutine run(this, n_chunk)

      use constants, only: T_MPI_SR
      use mpisetup,  only: tag_ub

      implicit none

      class(ring_t),       intent(inout) :: this     ! object invoking type-bound procedure
      integer(kind=INT64), intent(in)    :: n_chunk  ! number of pieces to communicate

      call this%setup

      select case (this%test_type)
         case (T_MPI_SR)
            if (n_chunk < tag_ub) then
               if (.not. sendrecv()) then
                  ! write(*, '(a)')"# Too big memory usage"
                  this%give_up = .true.
                  return
               endif
            else
               write(*, '(a)')"# Reached limit for MPI tags"
               return
            endif
         case default
            write(*,*)" Unknown test type: ", this%test_type
            call exit(-11)
      end select

      call this%check

   contains

      logical function sendrecv()

         use memory,   only: memcheck
         use mpi,      only: MPI_INTEGER_KIND, MPI_DOUBLE_PRECISION, MPI_LOGICAL, MPI_COMM_WORLD, &
              &              MPI_STATUSES_IGNORE, MPI_IN_PLACE, MPI_LAND, MPI_REQUEST_NULL
         use mpisetup, only: proc, nproc, ierr

         implicit none

         integer(kind=MPI_INTEGER_KIND), allocatable, dimension(:) :: req
         integer(kind=INT64) :: i, j

         sendrecv = .true.

         allocate(req(2*n_chunk))
         req = MPI_REQUEST_NULL
         do i = 1, n_chunk
            call MPI_Isend(this%sbuf(1+(i-1)*(this%n/n_chunk)), int(this%n/n_chunk, kind=MPI_INTEGER_KIND), MPI_DOUBLE_PRECISION, &
                 &         mod(nproc + proc - 1, nproc), i, MPI_COMM_WORLD, req(2*i-1), ierr)
            call MPI_Irecv(this%rbuf(1+(i-1)*(this%n/n_chunk)), int(this%n/n_chunk, kind=MPI_INTEGER_KIND), MPI_DOUBLE_PRECISION, &
                 &         mod(        proc + 1, nproc), i, MPI_COMM_WORLD, req(2*i  ), ierr)
            if (mod(i, 2_INT64**18) == 0) then
               sendrecv = memcheck()
               if (.not. sendrecv) exit
            endif
         enddo
         call MPI_Allreduce(MPI_IN_PLACE, sendrecv, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
         if (sendrecv) then
            call MPI_Waitall(size(req, kind=MPI_INTEGER_KIND), req, MPI_STATUSES_IGNORE, ierr)
         else
            do j = 1, 2*i
               call MPI_Cancel(req(j), ierr)
            enddo
         endif
         deallocate(req)

      end function sendrecv

   end subroutine run

end module ring
