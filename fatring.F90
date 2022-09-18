! Unlike most MPI ring demonstrations here we focus on transferring a huge
! amount of data. The transfer employs various MPI mechanisms and a range of
! message sizes is used. Each stage of the communication is profiled so the
! performance of the machine(s) and MPI implementation can be analyzed in detail.

program fat_ring

   use composition,     only: factorization_t
   use constants,       only: buflen, T_MPI_SR, T_MPI_GET1, T_MPI_PUT1, T_MPI_GETN, T_MPI_PUTN, T_MPI_SHMEM  !, V_SPEED, V_STATS
   use divisor,         only: factored_divisor
   use iso_fortran_env, only: output_unit, INT64, REAL64
   use mpisetup,        only: parallel_init, parallel_finalize, main_proc, proc
   use primes_utils,    only: primes_t

   implicit none

   ! constants
   real(kind=REAL64) :: two = 2_REAL64

   ! defaults for main parameters
   integer(kind=INT64), parameter :: n_doubles_def = 10810800  ! The amount of data to operate on, HCN(47), also HCN-3 (approx. 82.5 MB per process max.)
   integer(kind=INT64), parameter :: n_chunk_max_def = 2**15   ! Skip too big tests due to excessive allocations of memory inside MPI and generally slow execution
   real(kind=REAL64), parameter :: minskip_def = 2.**(1./3.)   ! Pick the triples that are no closer to each other than this ratio

   integer(kind=INT64), save :: n_doubles = n_doubles_def
   integer(kind=INT64), save :: n_chunk_max = n_chunk_max_def
   real(kind=REAL64), save :: minskip = minskip_def
   integer, save :: test_type = T_MPI_SR  !, verbosity = V_SPEED

   ! local variables
   integer(kind=INT64) :: i, p1, p2, d3
   real(kind=REAL64) :: a
   character(len=buflen) :: arg, buf
   type(factorization_t) :: n1, n2
   type(factored_divisor) :: d1, d2
   type(primes_t) :: primes
   integer(kind=INT64), parameter :: max_allowed_noncomposite = 1000000  ! prime sizes of the buffer up to 999983 are allowed
   integer :: ia, in
   enum, bind(C)
      enumerator :: A_ND, A_NC, A_MS
   end enum

   call parallel_init

   in = A_ND
   ia = 1
   do while (ia <= command_argument_count())
      call get_command_argument(ia, arg)
      select case (trim(arg))
         case ("-h")
            call usage("")
            call exit(-43)
         case ("-t")
            ia = ia + 1
            if (ia <= command_argument_count()) then
               call get_command_argument(ia, buf)
               select case (trim(buf))
                  case ("SR", "sr")
                     test_type = T_MPI_SR
                  case ("P1", "p1", "Put1")
                     test_type = T_MPI_PUT1
                     call usage("TEST_TYPE '" // trim(buf) // "' not implemented yet")
                     call exit(-41)
                  case ("G1", "g1", "Get1")
                     test_type = T_MPI_GET1
                  case ("PN", "pn", "PutN")
                     test_type = T_MPI_PUTN
                     call usage("TEST_TYPE '" // trim(buf) // "' not implemented yet")
                     call exit(-41)
                  case ("GN", "gn", "GetN")
                     test_type = T_MPI_GETN
                     call usage("TEST_TYPE '" // trim(buf) // "' not implemented yet")
                     call exit(-41)
                  case ("Sh", "SH", "sh", "shmem")
                     test_type = T_MPI_SHMEM
                     call usage("TEST_TYPE '" // trim(buf) // "' not implemented yet")
                     call exit(-41)
                  case default
                     call usage("Unrecognized TEST_TYPE provided: '" // trim(buf) // "'")
                     call exit(-41)
               end select
            else
               call usage("No TEST_TYPE provided")
               call exit(-37)
            end if
         case default
            select case (in)
               case (A_ND)
                  read(arg, '(i30)') i ! expect overflow on some 19-digit numbers but who cares?
                  if (i <= 0) then
                     call usage("Invalid input for size (non-negative integer expected)")
                     call exit(-23)
                  end if
                  n_doubles = i
               case (A_NC)
                  read(arg, '(i30)') i ! expect overflow on some 19-digit numbers but who cares?
                  if (i <= 0) then
                     call usage("Invalid input for chunks (non-negative integer expected)")
                     call exit(-29)
                  end if
                  n_chunk_max = i
               case (A_MS)
                  read(arg, '(f30.20)') a ! expect overflow on some 19-digit numbers but who cares?
                  if (a <= 1.) then
                     write(buf,'(a,f6.3,a)')"Invalid ratio (>=1. real expected, read: ", a, ")"
                     call usage(buf)
                     call exit(-31)
                  end if
                  minskip = a
               case default
                  write(*,'(a)')"Too many numbers provided. Ignored: '" // trim(buf) // "'"
            end select
            in = in + 1
      end select
      ia = ia + 1
   end do

   if (proc == main_proc) write(*,'(2a,i9,i7,f7.3)')"Called fatring -t ", trim(tt(test_type)), n_doubles, n_chunk_max, minskip

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

   ! Restrict the test to nicely divisible numbers to avoid computing excessively large list of primes.
   ! Calling primes%sieve(n_doubles) will be perfectlyalways correct but would take long time for GB-sized buffer.
   ! call primes%sieve(n_doubles)
   ! This will cause failure in factorize() for prime-sized buffer larger than max_allowed_noncomposite
   call primes%sieve(max(int(sqrt(real(n_doubles)), kind=INT64), min(n_doubles, max_allowed_noncomposite)))
   ! if (proc == main_proc) call primes%print

   !  Scan triangular region of chunk numbers and sizes within the limits, try to cover the accessible area evenly (in a relatively cheap way)
   call n1%factorize(n_doubles, primes)
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
         call n2%factorize(n1%number/d1%total(), primes)
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
   call primes%erase

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

   subroutine usage(extra_msg)

      use mpisetup, only: main_proc, proc

      implicit none

      character(len=*), intent(in) :: extra_msg

      if (proc /= main_proc) return

      write(*,'(a)')"Usage: fatring [-h] [-t SR|P1|G1|Pn|Gn|Sh] [number_of_doubles [maximum_number_of_chunks [minimum_sample_distance]]"
      write(*,'(a)')"       -h: help"
      write(*,'(a)')"       -t TEST_TYPE: choose test type"
      write(*,'(a)')"       TEST_TYPE: one of"
      write(*,'(a)')"           SR: Send/Recv"
      write(*,'(a)')"           P1/G1: Put/Get in single window"
      write(*,'(a)')"           Pn/Gn: Put/Get in many windows"
      write(*,'(a)')"           Sh: Shared memory"
      write(*,'(a)')""
      write(*,'(a,i9,i7,f7.3)')"defaults: -t SR ", n_doubles_def, n_chunk_max_def, minskip_def

      if (len_trim(extra_msg) > 0) then
         write(*,'(a)')""
         write(*,'(a)')trim(extra_msg)
      end if

   end subroutine usage

   function tt(type) result(t_str)

      implicit none

      integer, intent(in) :: type

      character(len=2) :: t_str

      t_str = "!!"
      select case(type)
         case (T_MPI_SR)
            t_str = "SR"
         case (T_MPI_GET1)
            t_str = "G1"
         case (T_MPI_PUT1)
            t_str = "P1"
         case (T_MPI_GETN)
            t_str = "Gn"
         case (T_MPI_PUTN)
            t_str = "Pn"
         case (T_MPI_SHMEM)
            t_str = "Sh"
         case default
            t_str = "??"
      end select

   end function tt

end program fat_ring
