! Actual executor of the ring communication.

module ring

   use iso_fortran_env, only: REAL64, INT64

   implicit none

   private
   public :: ring_t

   type :: ring_t
      real(kind=REAL64), allocatable, dimension(:) :: sbuf
      real(kind=REAL64), allocatable, dimension(:) :: rbuf
      integer, private :: test_type
      integer(kind=INT64), private :: n_chunk_max
      integer(kind=INT64), private :: n
      logical :: give_up
   contains
      procedure :: init
      procedure, private :: setup
      procedure :: run
      procedure, private :: check
      procedure :: cleanup
   end type ring_t

   enum, bind(C)
      enumerator :: F_ZERO, F_PROC, F_FANCY
   end enum
   integer, parameter :: fillstyle = F_FANCY

contains

   ! Initialize buffers.

   subroutine init(this, n, test_type, chunk_max)

      use composition, only: factorization_t
      use memory,      only: memcheck

      implicit none

      class(ring_t),         intent(inout) :: this       ! object invoking type-bound procedure
      type(factorization_t), intent(in)    :: n          ! number of doubles in the buffer
      integer,               intent(in)    :: test_type  ! test to perform
      integer(kind=INT64),   intent(in)    :: chunk_max  ! maximum allowed fragmentation

      this%give_up = .false.
      this%n = n%number
      this%test_type = test_type
      this%n_chunk_max = chunk_max

      call this%cleanup
      allocate(this%rbuf(this%n), &
           &   this%sbuf(this%n))

      if (.not. memcheck()) call exit(-13)

   end subroutine init

   ! Put some values into buffers.

   subroutine setup(this)

      use mpisetup, only: proc

      implicit none

      class(ring_t), intent(inout) :: this  ! object invoking type-bound procedure

      integer(kind=INT64) :: i

      select case (fillstyle)
         case (F_ZERO)
            this%sbuf = 0.
         case (F_PROC)
            this%sbuf = real(proc)
         case (F_FANCY)
            do i = lbound(this%sbuf, 1), ubound(this%sbuf, 1)
               this%sbuf(i) = real(mod(proc * i, 17_INT64))
            end do
         case default
            write(*,*)"[ring_t:setup] fillstyle not implemented"
            call exit(-17)
      end select
      this%rbuf = huge(1._REAL64)
   end subroutine setup

   ! Check if the received values are correct.

   subroutine check(this)

      use mpisetup, only: proc, nproc

      implicit none

      class(ring_t), intent(in) :: this  ! object invoking type-bound procedure

      integer(kind=INT64) :: i
      integer :: cnt

      cnt = 0
      select case (fillstyle)
         case (F_ZERO)
            cnt = count(this%rbuf /= 0.)
         case (F_PROC)
            cnt = count(this%rbuf /= real(mod(proc + 1, nproc)))
         case (F_FANCY)
            do i = lbound(this%rbuf, 1), ubound(this%rbuf, 1)
               if (this%rbuf(i) /= real(mod(mod(proc + 1, nproc) * i, 17_INT64))) cnt = cnt + 1
            end do
         case default
            write(*,*)"[ring_t:check] fillstyle not implemented"
            call exit(-19)
      end select
      if (cnt > 0) &
           write(*,*)"### ", cnt, "/", size(this%rbuf), " wrong values @", proc

   end subroutine check

   ! Clean up.

   subroutine cleanup(this)

      implicit none

      class(ring_t), intent(inout) :: this  ! object invoking type-bound procedure

      if (allocated(this%rbuf)) deallocate(this%rbuf)
      if (allocated(this%sbuf)) deallocate(this%sbuf)

   end subroutine cleanup

   ! Do the communication in the specified mode.

   subroutine run(this, n_chunk)

      use constants, only: T_MPI_SR, buflen
      use mpi,       only: MPI_Wtime
      use mpisetup,  only: tag_ub, proc, main_proc

      implicit none

      class(ring_t),       intent(inout) :: this     ! object invoking type-bound procedure
      integer(kind=INT64), intent(in)    :: n_chunk  ! number of pieces to communicate
      enum, bind(C)
         enumerator :: W_START, W_SETUP, W_SR, W_WAITALL, W_CHECK
      end enum
      double precision, dimension(W_START:W_CHECK) :: wtime
      character(len=buflen) :: buf
      integer :: i

      if (proc == main_proc) &
           write(*, '(a5,2a12,5a13)') "#proc", "chunks", "doubles", "T_setup", "T_SendRecv", "T_Waitall", "T_check", "MiB/s"

      wtime(W_START) = MPI_Wtime()

      call this%setup
      wtime(W_SETUP) = MPI_Wtime()

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
      wtime(W_CHECK) = MPI_Wtime()

      write(buf, '(i5,2i12)')proc, n_chunk, this%n/n_chunk
      do i = lbound(wtime, dim=1) + 1, ubound(wtime, dim=1)
         write(buf, '(a," ",f12.6)') trim(buf), wtime(i) - wtime(i-1)
      enddo
      associate (dt => wtime(W_WAITALL) - wtime(W_SETUP))
         if (dt > 0.) write(buf, '(a," ",f12.3)') trim(buf), real(this%n) / 2.**17 / dt
      end associate
      write(*, '(a)') trim(buf)

   contains

      ! Use Isend/Irecv, the default mode.

      logical function sendrecv()

         use memory,   only: memcheck
         use mpi,      only: MPI_INTEGER_KIND, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, MPI_STATUSES_IGNORE
         use mpisetup, only: proc, nproc, ierr

         implicit none

         integer(kind=MPI_INTEGER_KIND), allocatable, dimension(:) :: req
         integer(kind=INT64) :: i

         sendrecv = (n_chunk <= this%n_chunk_max)
         if (.not. sendrecv) return

         allocate(req(2*n_chunk))

         do i = 1, n_chunk
            call MPI_Isend(this%sbuf(1+(i-1)*(this%n/n_chunk)), int(this%n/n_chunk, kind=MPI_INTEGER_KIND), MPI_DOUBLE_PRECISION, &
                 &         mod(nproc + proc - 1, nproc), i, MPI_COMM_WORLD, req(2*i-1), ierr)
            call MPI_Irecv(this%rbuf(1+(i-1)*(this%n/n_chunk)), int(this%n/n_chunk, kind=MPI_INTEGER_KIND), MPI_DOUBLE_PRECISION, &
                 &         mod(        proc + 1, nproc), i, MPI_COMM_WORLD, req(2*i  ), ierr)
         enddo
         wtime(W_SR) = MPI_Wtime()

         call MPI_Waitall(size(req, kind=MPI_INTEGER_KIND), req, MPI_STATUSES_IGNORE, ierr)
         wtime(W_WAITALL) = MPI_Wtime()

         deallocate(req)

      end function sendrecv

   end subroutine run

end module ring
