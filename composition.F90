! Decompose given integer into prime factors

! To achieve good coverage of the message size spectrum it is advised to use
! a buffer length that is nicely divisible by many different factors.
! Some candidates are:
!     2**n  - for fast, coarse scans (~3 point/decade)
!     3*2**n  - less coarse scan (6-7 points per decade, default)
!     3*5*7*2**n  - fine scan (~13 points per decade)
!     43243200, 61261200 or other highly composite number for best resolution (and slowest scan)

! For a list of Highly Composite Numbers, see e.g.:
! http://wwwhomes.uni-bielefeld.de/achim/highly.txt

module composition

   use constants, only: INT64

   implicit none

   private
   public :: factorization_t

   type component
      integer(kind=INT64) :: prime
      integer :: power
   end type component

   type factorization_t
      integer(kind=INT64) :: number
      type(component), allocatable, dimension(:) :: factors
   contains
      procedure :: factorize
      procedure :: erase
   end type factorization_t

contains

   subroutine factorize(this, n)

      use constants,    only: buflen
      use primes_utils, only: primes_t

      implicit none

      class(factorization_t), intent(inout) :: this  ! object invoking type-bound procedure
      integer(kind=INT64),    intent(in)    :: n     ! the number to factorize

      type(primes_t) :: primes
      integer :: i
      integer(kind=INT64) :: auxn
      character(len=buflen) :: buf, n1, n2

      this%number = n
      call primes%sieve(int(sqrt(real(this%number)), kind=INT64))
      !call primes%print

      allocate(this%factors(0))
      auxn = this%number
      do i = lbound(primes%tab, dim=1), ubound(primes%tab, dim=1)
         if (mod(auxn, primes%tab(i)) == 0) then
            this%factors = [ this%factors, component(primes%tab(i), 0) ]
            do while (mod(auxn, primes%tab(i)) == 0)
               auxn = auxn / primes%tab(i)
               this%factors(ubound(this%factors, dim=1))%power = this%factors(ubound(this%factors, dim=1))%power + 1
            enddo
         endif
         if (auxn <= 1) exit
      enddo
      if (auxn > 1) then
         write(*,*)"Decomposition aborted due to prime factors exceeding sqrt(number of doubles). Please try with more composite number"
         call exit(-7)
      endif

      auxn = 1
      write(buf, *) this%number, " = "
      do i = lbound(this%factors, dim=1), ubound(this%factors, dim=1)
         write(n1, *) this%factors(i)%prime
         write(n2, *) this%factors(i)%power
         write(buf, *) trim(adjustl(buf)), " ", trim(adjustl(n1)), "**", trim(adjustl(n2)), merge("  ", " *", i == ubound(this%factors, dim=1))
         auxn = auxn * this%factors(i)%prime**this%factors(i)%power
      enddo
      if (auxn == this%number) then
         write(*,*) trim(adjustl(buf))
      else
         write(*,*)"Decomposition failed: ", auxn, " /= ", this%number
         call exit(-3)
      endif

      call primes%erase

   end subroutine factorize

   subroutine erase(this)

      use constants, only: INVALID

      implicit none

      class(factorization_t), intent(inout) :: this  ! object invoking type-bound procedure

      this%number = INVALID
      if (allocated(this%factors)) deallocate(this%factors)

   end subroutine erase

end module composition
