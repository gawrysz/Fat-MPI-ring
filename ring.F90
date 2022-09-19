! Actual executor of the ring communication.

module ring

   use constants,       only: lablen
   use iso_fortran_env, only: REAL64, INT64

   implicit none

   private
   public :: ring_t

   type :: tmr_t  ! Keeps timers with labels
      character(len=lablen) :: label
      double precision :: wtime
   end type tmr_t

   type :: ring_t
      private
      real(kind=REAL64), allocatable, dimension(:) :: sbuf
      real(kind=REAL64), allocatable, dimension(:) :: rbuf
      integer :: test_type
      integer(kind=INT64) :: n_chunk_max
      integer(kind=INT64) :: n
      logical, public :: give_up
      type(tmr_t) :: w_start, w_setup, w_check
      type(tmr_t), dimension(:), allocatable :: w_tst
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

      use constants, only: T_MPI_SR, T_MPI_GET1, T_MPI_PUT1, T_MPI_GETN, T_MPI_PUTN
      use memory,    only: memcheck

      implicit none

      class(ring_t),       intent(inout) :: this       ! object invoking type-bound procedure
      integer(kind=INT64), intent(in)    :: n          ! number of doubles in the buffer
      integer,             intent(in)    :: test_type  ! test to perform
      integer(kind=INT64), intent(in)    :: chunk_max  ! maximum allowed fragmentation

      this%give_up = .false.
      this%n = n
      this%test_type = test_type
      this%n_chunk_max = chunk_max

      call this%cleanup
      allocate(this%rbuf(this%n), &
           &   this%sbuf(this%n))

      if (.not. memcheck()) call exit(-13)

      this%w_start%label = "T_start"  ! unused
      this%w_setup%label = "T_setup"
      this%w_check%label = "T_check"

      select case (this%test_type)
         case (T_MPI_SR)
            allocate(this%w_tst(2))
            this%w_tst(1)%label = "T_SendRecv"
            this%w_tst(2)%label = "T_Waitall"
         case (T_MPI_GET1, T_MPI_PUT1, T_MPI_GETN, T_MPI_PUTN)
            allocate(this%w_tst(5))
            this%w_tst(1)%label = "T_WinCreate"
            this%w_tst(2)%label = "T_WinFence1"
            this%w_tst(3)%label = merge("T_Get", "T_Put", any(this%test_type == [T_MPI_GET1, T_MPI_GETN]))
            this%w_tst(4)%label = "T_WinFence2"
            this%w_tst(5)%label = "T_WinFree"
         case default
            write(*,*)"[ring:init] Unknown test type: ", this%test_type
            call exit(-47)
      end select

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
            write(*,*)"[ring:setup] fillstyle not implemented"
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
            write(*,*)"[ring:check] fillstyle not implemented"
            call exit(-19)
      end select
      if (cnt > 0) then
         write(*,*)"### ", cnt, "/", size(this%rbuf), " wrong values @", proc
         call exit(-5)
      end if

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

      use constants, only: T_MPI_SR, T_MPI_GET1, T_MPI_PUT1, T_MPI_GETN, T_MPI_PUTN
      use mpi,       only: MPI_Wtime
      use mpisetup,  only: tag_ub

      implicit none

      class(ring_t),       intent(inout) :: this     ! object invoking type-bound procedure
      integer(kind=INT64), intent(in)    :: n_chunk  ! number of pieces to communicate

      this%w_start%wtime = MPI_Wtime()

      call this%setup
      this%w_setup%wtime = MPI_Wtime()

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
         case (T_MPI_GET1, T_MPI_PUT1)
            call getput1
         case (T_MPI_GETN, T_MPI_PUTN)
            call getputn
         case default
            write(*,*)"[ring:run] Unknown test type: ", this%test_type
            call exit(-11)
      end select

      call this%check
      this%w_check%wtime = MPI_Wtime()

      call print_wtime([this%w_start, this%w_setup, this%w_tst, this%w_check])

   contains

      ! Collect wtime on main_proc, print it in order and also print extrema

      subroutine print_wtime(wt)

         use constants, only: buflen
         use mpi,       only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, MPI_STATUS_IGNORE
         use mpisetup,  only: main_proc, nproc, proc, ierr

         implicit none

         type(tmr_t), dimension(:), intent(in) :: wt

         integer :: p
         character(len=buflen) :: buf
         double precision, dimension(size(wt)) :: p_wtime, dt_max, dt_min, dt
         integer :: i

         if (proc == main_proc) then
            write(buf, '(a5,2a12)') "#proc", "chunks", "doubles"
            do i = lbound(wt, dim=1) + 1, ubound(wt, dim=1)  ! skip T_start
               write(buf(len_trim(buf) + 1:), '(a14)') trim(wt(i)%label)
            end do
            write(buf(len_trim(buf) + 1:), '(a14)') "MiB/s"
            write(*, '(a)') trim(buf)

            dt_min(:) = huge(1.)
            dt_max(:) = -huge(1.)
            do p = main_proc, nproc-1
               if (p == main_proc) then
                  p_wtime(:) = wt(:)%wtime
               else
                  call MPI_Recv(p_wtime, size(p_wtime), MPI_DOUBLE_PRECISION, p, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
               end if

               dt(lbound(p_wtime, dim=1)) = p_wtime(ubound(wt, dim=1) - 1) - p_wtime(lbound(wt, dim=1) + 1)
               do i = lbound(p_wtime, dim=1) + 1, ubound(p_wtime, dim=1)
                  dt(i) = p_wtime(i) - p_wtime(i-1)
               enddo

               dt_min(:) = min(dt_min(:), dt(:))
               dt_max(:) = max(dt_max(:), dt(:))

               write(buf, '(i5,2i12)')p, n_chunk, this%n/n_chunk
               call pr_wt(dt, buf)
            end do
            write(buf, '(a5,2i12)')"min", n_chunk, this%n/n_chunk
            call pr_wt(dt_min, buf)
            write(buf, '(a5,2i12)')"max", n_chunk, this%n/n_chunk
            call pr_wt(dt_max, buf)
         else
            p_wtime(:) = wt(:)%wtime
            call MPI_Send(p_wtime, size(p_wtime), MPI_DOUBLE_PRECISION, main_proc, proc, MPI_COMM_WORLD, ierr)
         end if

      end subroutine print_wtime

      ! Extracted for convenience

      subroutine pr_wt(dt, buf)

         implicit none

         double precision, dimension(:), intent(in)    :: dt
         character(len=*),               intent(inout) :: buf

         integer :: i

         do i = lbound(dt, dim=1) + 1, ubound(dt, dim=1)
            write(buf, '(a," ",f13.7)') trim(buf), dt(i)
         enddo
         if (dt(lbound(dt, dim=1)) > 0.) write(buf, '(a," ",f13.3)') trim(buf), real(this%n) / 2.**17 / dt(lbound(dt, dim=1))
         write(*, '(a)') trim(buf)

      end subroutine pr_wt

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
         this%w_tst(1)%wtime = MPI_Wtime()

         call MPI_Waitall(size(req, kind=MPI_INTEGER_KIND), req, MPI_STATUSES_IGNORE, ierr)
         this%w_tst(2)%wtime = MPI_Wtime()

         deallocate(req)

      end function sendrecv

      subroutine getput1

         use mpi,      only: MPI_DOUBLE_PRECISION, MPI_INTEGER_KIND, MPI_ADDRESS_KIND, MPI_INFO_NULL, MPI_COMM_WORLD
         use mpisetup, only: proc, nproc, ierr

         implicit none

         integer(kind=INT64) :: i
         integer(kind=MPI_INTEGER_KIND) :: lb, d_ext, ww

         call MPI_Type_get_extent(MPI_DOUBLE_PRECISION, lb, d_ext, ierr)

         select case (this%test_type)
            case (T_MPI_GET1)
               call MPI_Win_create(this%sbuf, size(this%sbuf, kind=MPI_ADDRESS_KIND) * d_ext, d_ext, MPI_INFO_NULL, MPI_COMM_WORLD, ww, ierr)
            case (T_MPI_PUT1)
               call MPI_Win_create(this%rbuf, size(this%rbuf, kind=MPI_ADDRESS_KIND) * d_ext, d_ext, MPI_INFO_NULL, MPI_COMM_WORLD, ww, ierr)
            case default
               write(*,*)"[ring:run:getput1] test type not implemented (win create)"
               call exit(-59)
         end select
         this%w_tst(1)%wtime = MPI_Wtime()

         call MPI_Win_fence(0, ww, ierr)
         this%w_tst(2)%wtime = MPI_Wtime()

         do i = 1, n_chunk
            select case (this%test_type)
               ! Note that despite of indexing from 1 on local array, we have to provide target_disp starting from 0
               case (T_MPI_GET1)
                  call MPI_Get(this%rbuf(1+(i-1)*(this%n/n_chunk)), int(this%n/n_chunk, kind=MPI_ADDRESS_KIND), MPI_DOUBLE_PRECISION, &
                       &       mod(proc + 1, nproc), int((i-1)*(this%n/n_chunk), kind=MPI_ADDRESS_KIND), int(this%n/n_chunk, kind=MPI_ADDRESS_KIND), &
                       &       MPI_DOUBLE_PRECISION, ww, ierr)
               case (T_MPI_PUT1)
                  call MPI_Put(this%sbuf(1+(i-1)*(this%n/n_chunk)), int(this%n/n_chunk, kind=MPI_ADDRESS_KIND), MPI_DOUBLE_PRECISION, &
                       &       mod(nproc + proc - 1, nproc), int((i-1)*(this%n/n_chunk), kind=MPI_ADDRESS_KIND), int(this%n/n_chunk, kind=MPI_ADDRESS_KIND), &
                       &       MPI_DOUBLE_PRECISION, ww, ierr)
               case default
                  write(*,*)"[ring:run:getput1] test type not implemented (get/put)"
                  call exit(-53)
            end select
         end do
         this%w_tst(3)%wtime = MPI_Wtime()

         call MPI_Win_fence(0, ww, ierr)
         this%w_tst(4)%wtime = MPI_Wtime()

         call MPI_Win_free(ww, ierr)
         this%w_tst(5)%wtime = MPI_Wtime()

      end subroutine getput1

      subroutine getputn

         use mpi,      only: MPI_DOUBLE_PRECISION, MPI_INTEGER_KIND, MPI_ADDRESS_KIND, MPI_INFO_NULL, MPI_COMM_WORLD
         use mpisetup, only: proc, nproc, ierr

         implicit none

         integer(kind=INT64) :: i
         integer(kind=MPI_INTEGER_KIND) :: lb, d_ext
         integer(kind=MPI_INTEGER_KIND), dimension(:), allocatable :: ww

         call MPI_Type_get_extent(MPI_DOUBLE_PRECISION, lb, d_ext, ierr)
         allocate(ww(n_chunk))

         do i = 1, n_chunk
            select case (this%test_type)
               case (T_MPI_GETN)
                  call MPI_Win_create(this%sbuf(1+(i-1)*(this%n/n_chunk)), int(this%n/n_chunk * d_ext, kind=MPI_ADDRESS_KIND), d_ext, MPI_INFO_NULL, MPI_COMM_WORLD, ww(i), ierr)
               case (T_MPI_PUTN)
                  call MPI_Win_create(this%rbuf(1+(i-1)*(this%n/n_chunk)), int(this%n/n_chunk * d_ext, kind=MPI_ADDRESS_KIND), d_ext, MPI_INFO_NULL, MPI_COMM_WORLD, ww(i), ierr)
               case default
                  write(*,*)"[ring:run:getputn] test type not implemented (win create)"
                  call exit(-59)
            end select
         end do
         this%w_tst(1)%wtime = MPI_Wtime()

         do i = 1, n_chunk
            call MPI_Win_fence(0, ww(i), ierr)
         end do
         this%w_tst(2)%wtime = MPI_Wtime()

         do i = 1, n_chunk
            select case (this%test_type)
               ! Note that despite of indexing from 1 on local array, we have to provide target_disp starting from 0
               case (T_MPI_GETN)
                  call MPI_Get(this%rbuf(1+(i-1)*(this%n/n_chunk)), int(this%n/n_chunk, kind=MPI_ADDRESS_KIND), MPI_DOUBLE_PRECISION, &
                       &       mod(proc + 1, nproc), int(0, kind=MPI_ADDRESS_KIND), int(this%n/n_chunk, kind=MPI_ADDRESS_KIND), &
                       &       MPI_DOUBLE_PRECISION, ww(i), ierr)
               case (T_MPI_PUTN)
                  call MPI_Put(this%sbuf(1+(i-1)*(this%n/n_chunk)), int(this%n/n_chunk, kind=MPI_ADDRESS_KIND), MPI_DOUBLE_PRECISION, &
                       &       mod(nproc + proc - 1, nproc), int(0, kind=MPI_ADDRESS_KIND), int(this%n/n_chunk, kind=MPI_ADDRESS_KIND), &
                       &       MPI_DOUBLE_PRECISION, ww(i), ierr)
               case default
                  write(*,*)"[ring:run:getputn] test type not implemented (get/put)"
                  call exit(-53)
            end select
         end do
         this%w_tst(3)%wtime = MPI_Wtime()

         do i = 1, n_chunk
            call MPI_Win_fence(0, ww(i), ierr)
         end do
         this%w_tst(4)%wtime = MPI_Wtime()

         do i = 1, n_chunk
            call MPI_Win_free(ww(i), ierr)
         end do
         this%w_tst(5)%wtime = MPI_Wtime()

         deallocate(ww)

      end subroutine getputn

   end subroutine run

end module ring
