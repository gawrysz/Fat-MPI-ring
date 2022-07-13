! Unlike most MPI ring demonstrations here we focus on transferring a huge
! amount of data. The transfer employs various MPI mechanisms and a range of
! message sizes is used. Each stage of the communication is profiled so the
! performance of the machine(s) and MPI implementation can be analyzed in detail.

program fat_ring

   use composition,     only: factorization_t
   use constants,       only: buflen, T_MPI_SR  !, V_SPEED, V_STATS
   use divisor,         only: factored_divisor
   use iso_fortran_env, only: output_unit, INT64, REAL64
   use mpi,             only: MPI_Barrier, MPI_COMM_WORLD
   use mpisetup,        only: parallel_init, parallel_finalize, main_proc, proc, ierr
   use ring,            only: ring_t

   implicit none

   ! constants
   real(kind=REAL64) :: two = 2_REAL64
   integer :: n_chunk_max = 2**15  ! Skip too big tests due to excessive allocations of memory inside MPI and generally slow execution

   ! defaults for main parameters
!   integer(kind=INT64), save :: n_doubles = 5**2 * 2**19  ! the amount of data to operate on (approx. 100 MB per process)
   integer(kind=INT64), save :: n_doubles = 5**2 * 2**15  ! the amount of data to operate on (approx. 6 MB per process)
   integer, save :: test_type = T_MPI_SR  !, verbosity = V_SPEED

   ! local variables
   integer(kind=INT64) :: i
   character(len=buflen) :: arg, buf
   type(factorization_t) :: n
   type(factored_divisor) :: d
   type(ring_t) :: r

   call parallel_init

   ! ToDo parse arguments
   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg)
      read(arg, '(i30)') i ! expect overflow on some 19-digit numbers but who cares?
      if (i > 0) n_doubles = i
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

   call n%factorize(n_doubles)

   if (proc == main_proc) &
        write(*, '(7a)')"# Starting fat MPI ring test with ", trim(adjustl(arg)), " == ", trim(n%factor_str), " doubles (", trim(adjustl(buf)), ")"

   call d%reset(n)
   i = 1
   call r%init(n, test_type, n_chunk_max)
   do while (d%is_valid())
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (.not. r%give_up) then
         if (proc == main_proc) &
              write(*, '(/,a,i4,a,2(i10,a))')"# Test ", i, ": " , d%total(), " chunks of ", n%number/d%total(), " doubles"
         call r%run(d%total())
      endif
      flush(output_unit)
      i = i + 1
      call d%next_div
   enddo
   if (r%give_up .and. (proc == main_proc)) write(*, '(a)')"# Some tests have been skipped"
   call r%cleanup
   call d%clear
   call n%erase

!   call parallel_finalize
   if (proc == main_proc) write(*, '(a)')"# End."

end program fat_ring
