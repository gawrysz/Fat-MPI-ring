! This module serves an array of prime numbers.
! Adopted from Piernik:
!     https://github.com/piernik-dev/piernik/blob/master/src/base/primes.F90

module primes_utils

   use iso_fortran_env, only: INT64

   implicit none

   private
   public :: primes_t

   ! The type that contain table of primes and initializer
   type :: primes_t
      integer(kind=INT64), allocatable, dimension(:) :: tab  ! table of prime numbers
      integer(kind=INT64), private :: max                    ! max number to which the search was performed
   contains
      procedure :: sieve  ! routine used to initialize (and extendif necessary) the table of prime numbers
      procedure :: erase  ! restore initial state
      procedure :: print  ! print what was found
   end type primes_t

contains

! Eratosthenes sieve.
!
! On first call it creates a table of primes up to a given number (does an initialization).
! On subsequent calls it may extend the table with primes, if required.

   subroutine sieve(this, n)

      implicit none

      class(primes_t),     intent(inout) :: this  ! object invoking type-bound procedure
      integer(kind=INT64), intent(in)    :: n     ! max value of prime number

      integer(kind=INT64), dimension(n) :: numb   ! here we could use a 1-byte integers but then we would
                                                  ! need to write something more elaborate than pack()
      integer(kind=INT64) :: i

      if (allocated(this%tab)) then ! we may expect this%max was already set
         if (n <= this%max) return
         deallocate(this%tab)
      endif

      if (n < 1) then
         write(*,*) "Prime numbers are usually greater than 1 ..."
         call exit(-2)
      endif

      ! The sieve
      numb = [0_INT64, (i, i = 2_INT64, n)]
      do i = 2_INT64, n
         if (numb(i) /= 0) numb( 2_INT64*i : n : i ) = 0
      enddo

      ! Gather the primes into a table
      allocate(this%tab(count(numb /= 0)))
      this%tab(:) = pack(numb, numb /= 0)
      this%max = n

   end subroutine sieve

! Print the table of found primes.

   subroutine print(this)

      implicit none

      class(primes_t), intent(inout) :: this  ! object invoking type-bound procedure

      write(*, '(a,i5,a,i9,a)') "There are ", size(this%tab), " prime numbers smaller than", this%max,":"
      print "(15i11)", this%tab

   end subroutine print

! Restore initial state.

   subroutine erase(this)

      use constants, only: INVALID

      implicit none

      class(primes_t), intent(inout) :: this  ! object invoking type-bound procedure

      deallocate(this%tab)
      this%max = INVALID

   end subroutine erase

end module primes_utils
