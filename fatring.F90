! Unlike most MPI ring demonstrations here we focus on transferring a huge
! amount of data. The transfer employs various MPI mechanisms and a range of
! message sizes is used. Each stage of the communication is profiled so the
! performance of the machine(s) and MPI implementation can be analyzed in detail.

program fat_ring

   use composition,     only: factorization_t
   use constants,       only: buflen, T_MPI_SR  !, V_SPEED, V_STATS
   use divisor,         only: factored_divisor
   use iso_fortran_env, only: output_unit, INT64, REAL64
   use mpisetup,        only: parallel_init, parallel_finalize, main_proc, proc

   implicit none

   ! constants
   real(kind=REAL64) :: two = 2_REAL64

   ! defaults for main parameters
   integer(kind=INT64), save :: n_doubles = 10810800  ! The amount of data to operate on, HCN(47), also HCN-3 (approx. 82.5 MB per process max.)
   integer(kind=INT64), save :: n_chunk_max = 2**15   ! Skip too big tests due to excessive allocations of memory inside MPI and generally slow execution
   real(kind=REAL64), save :: minskip = 2.**(1./3.)   ! Pick the triples that are no closer to each other than this ratio

   integer, save :: test_type = T_MPI_SR  !, verbosity = V_SPEED

   ! local variables
   integer(kind=INT64) :: i, p1, p2, d3
   real(kind=REAL64) :: a
   character(len=buflen) :: arg, buf
   type(factorization_t) :: n1, n2
   type(factored_divisor) :: d1, d2
   character(len=*), parameter :: usage_str = "Usage: fatring [number_of_doubles [maximum_number_of_chunks [minimum_sample_distance]]"

   call parallel_init

   ! ToDo parse arguments
   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg)
      read(arg, '(i30)') i ! expect overflow on some 19-digit numbers but who cares?
      if (i <= 0) then
         if (proc == main_proc) then
            write(*,*) usage_str
            write(*,*)"Invalid input for size (non-negative integer expected)"
         end if
         call exit(-23)
      end if
      n_doubles = i
   endif

   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg)
      read(arg, '(i30)') i ! expect overflow on some 19-digit numbers but who cares?
      if (i <= 0) then
         if (proc == main_proc) then
            write(*,*) usage_str
            write(*,*)"Invalid input for chunks (non-negative integer expected)"
         end if
         call exit(-29)
      end if
      n_chunk_max = i
   endif

   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg)
      read(arg, '(f30.20)') a ! expect overflow on some 19-digit numbers but who cares?
      if (a <= 1.) then
         if (proc == main_proc) then
            write(*,*) usage_str
            write(*,*)"Invalid ratio (>=1. real expected, read: ", a, ")"
         end if
         call exit(-31)
      end if
      minskip = a
   endif

   write(arg, *) n_doubles
   select case (n_doubles*8_INT64)  ! Warning: 8-byte sizeof(double) assumed
      case (2**30:2_INT64**40-1)
         write(buf, '(f6.1,a)') n_doubles / two**27, " GiB"
      case (2**20:2**30-1)
         write(buf, '(f6.1,a)') n_doubles / two**17, " MiB"
      case (2**10:2**20-1)
         write(buf, '(f6.1,a)') n_doubles / two**7, " kiB"
      case default
         write(buf, '(i20,a)') int(n_doubles * two**3, kind=INT64), " B"
   end select

   !  Scan triangular region of chunk numer and sizes within the limits, try to cover the accessible area as evenly as possible
   call n1%factorize(n_doubles)
   ! Need something better
   !flush(output_unit)
   !call MPI_Barrier(MPI_COMM_WORLD, ierr)
   if (proc == main_proc) then
      write(*, '(/,7a)')"# Starting fat MPI ring test with ", trim(adjustl(arg)), " == ", trim(n1%factor_str), " doubles (", trim(adjustl(buf)), ")"
      write(*, '(a,g10.3)')"# Minimum sample spacing ratio = ", minskip
   end if
   i = 1
   p1 = 0
   call d1%reset(n1)
   do while (d1%is_valid())
      ! Let's scan through all triplets that make n_doubles and its divisors but skip some of them to limit density of the points
      if (real(d1%total()) >= minskip * p1) then
         call n2%factorize(n1%number/d1%total())
         p2 = 0
         call d2%reset(n2)
         do while (d2%is_valid())
            if (real(d2%total()) >= minskip * p2) then
               d3 = n1%number/(d1%total()*d2%total())
               if (d3 >= d1%total() .and.  d2%total() >= d1%total()) then
                  call run_test(d1%total(), d2%total()) !, d3)
                  if (d2%total() >= d1%total() * sqrt(minskip) .and. d3 >= d1%total() * sqrt(minskip)) &
                       call run_test(d2%total(), d1%total()) !, d3)
                  if (d3 >= d1%total() * sqrt(minskip)) &
                       call run_test(d3, d2%total()) !, d1%total())
               end if
               p2 = d2%total()
            endif
            call d2%next_div
         enddo
         call d2%clear
         call n2%erase
         p1 = d1%total()
      endif
      call d1%next_div
   enddo
   call d1%clear
   call n1%erase

   call parallel_finalize
   if (proc == main_proc) write(*, '(a)')"# End."

contains

   subroutine run_test(ch_size, ch_num) !, c)

      use mpi,      only: MPI_Barrier, MPI_COMM_WORLD
      use mpisetup, only: ierr
      use ring,     only: ring_t

      implicit none

      integer(kind=INT64), intent(in) :: ch_size  ! Chunk size
      integer(kind=INT64), intent(in) :: ch_num   ! Number of chunks
      !, c  ! c * ch_size * ch_num == n_doubles

      integer, save :: cnt = 1
      type(ring_t) :: r

      if (ch_num < n_chunk_max) then

         call r%init(ch_size * ch_num, test_type, n_chunk_max)
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         if (proc == main_proc) &
              write(*, '(/,a,i4,a,2(i10,a))')"# Test ", cnt, ": " , ch_num, " chunk" // trim(merge("s", " ", ch_num>1)) // " of ", ch_size, " doubles"
         call r%run(ch_num)
         cnt = cnt + 1
         if (r%give_up .and. (proc == main_proc)) write(*, '(a)')"# Some tests have been skipped"
         call r%cleanup
      end if

   end subroutine run_test

end program fat_ring
