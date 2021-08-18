! Unlike most MPI ring demonstrations here we focus on transferring a huge
! amount of data. The transfer employs various MPI mechanisms and a range of
! message sizes is used. Each stage of the communication is profiled so the
! performance of the machine(s) and MPI implementation can be analyzed in detail.

program fat_ring

   use composition, only: factorization_t
   use constants,   only: INT64, FP64, buflen

   implicit none

   ! constants
   enum, bind(C)
      enumerator :: V_NONE = 0, V_SPEED, V_STATS, V_DETAILED  ! verbosity levels
   end enum
   enum, bind(C)
      enumerator :: T_MPI = 1 !, T_COARRAY = 2 * T_MPI, T_OMP = 2 * T_COARRAY  ! masks for test types to perform
   end enum
   real(kind=FP64) :: two = 2_FP64

   ! main parameters
   integer(kind=INT64) :: n_doubles = 3 * 2**24  ! the amount of data to operate on
   integer, save :: verbosity = V_SPEED, test_mask = T_MPI

   ! local variables
   integer(kind=INT64) :: i
   character(len=buflen) :: arg, buf
   type(factorization_t) :: n

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
   write(*,*) "Starting fat MPI ring test with ", trim(adjustl(arg)), " doubles (", trim(adjustl(buf)), ")"

   call n%factorize(n_doubles)

   call n%erase

   if (.false.) i = verbosity * test_mask  ! temporarily suppress -Wunused-variable

   write(*,*)"End."

end program fat_ring
