! Unlike most MPI ring demonstrations here we focus on transferring a huge
! amount of data. The transfer employs various MPI mechanisms and a range of
! message sizes is used. Each stage of the communication is profiled so the
! performance of the machine(s) and MPI implementation can be analyzed in detail.

program fat_ring

   use composition, only: factorization_t
   use constants,   only: INT64, FP64, buflen, V_SPEED, V_STATS, T_MPI
   use divisor,     only: factored_divisor
   use mpi,         only: MPI_Barrier, MPI_COMM_WORLD
   use mpisetup,    only: parallel_init, parallel_finalize, main_proc, proc, ierr

   implicit none

   ! constants
   real(kind=FP64) :: two = 2_FP64

   ! defaults for main parameters
   integer(kind=INT64), save :: n_doubles = 5**2 * 2**21  ! the amount of data to operate on
   integer, save :: verbosity = V_SPEED, test_mask = T_MPI

   ! local variables
   integer(kind=INT64) :: i
   character(len=buflen) :: arg, buf
   type(factorization_t) :: n
   type(factored_divisor) :: d

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
   do while (d%is_valid())
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if ((proc == main_proc) .or. verbosity >= V_STATS) &
         write(*,'(2(a,i4),a,2(i10,a))')"# Test ", i, " @", proc, ": " , d%total(), " chunks of ", n%number/d%total(), " doubles"
      i = i + 1
      call d%next_div
   enddo
   call d%clear
   call n%erase

   if (.false.) i = test_mask  ! temporarily suppress -Wunused-variable

   call parallel_finalize
   if (proc == main_proc) write(*, '(a)')"# End."

end program fat_ring
