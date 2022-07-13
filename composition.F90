! Decompose given integer into prime factors

! To achieve good coverage of the message size spectrum it is advised to use
! a buffer length that is nicely divisible by many different factors.
! Some candidates are:
!     2**n, 3**n              for fast, coarse scans (~3 and ~2 point/decade, respectively)
!     3*2**n                  less coarse scan (6-7 points per decade, 2nd favourite)
!     11*2**n, 23*2**n        similar to above but a bit better for logarithmic scales
!     13**2*3**n              similar to above, free of powers of 2, very even on logarithmic scale
!     5**2*2**n               ~10 points per decade, very even, default
!     3*7*2**n                fine scan (~13 points per decade)
!     11*19*2**n, 23*19*2**n  similar to above but a bit better for logarithmic scales
!     43243200, 61261200      (or other highly composite number) for best resolution, especially in the middle (and slowest scan)

! For a list of Highly Composite Numbers, see e.g.:
! http://wwwhomes.uni-bielefeld.de/achim/highly.txt

module singleprime

   use iso_fortran_env, only: INT64

   implicit none

   private
   public :: component

   type :: component
      integer(kind=INT64) :: prime
      integer :: power
   contains
      procedure :: val
   end type component

contains

   integer(kind=INT64) pure function val(this)

      use constants, only: INVALID

      implicit none

      class(component), intent(in) :: this  ! object invoking type-bound procedure

      if (this%power > INVALID) then
         val = this%prime ** this%power
      else
         val = INVALID
      endif

   end function val

end module singleprime

module setofprimes

   use iso_fortran_env, only: INT64
   use singleprime,     only: component

   implicit none

   private
   public :: component_set

   type :: component_set
      type(component), allocatable, dimension(:) :: factors
   contains
      procedure :: total
   end type component_set

contains

   integer(kind=INT64) pure function total(this)

      use constants, only: INVALID

      implicit none

      class(component_set), intent(in) :: this  ! object invoking type-bound procedure

      integer :: i

      if (allocated(this%factors)) then
         total = 1
         do i = lbound(this%factors, dim=1), ubound(this%factors, dim=1)
            if (this%factors(i)%val() > INVALID) then
               total = total * this%factors(i)%val()
            else
               total = INVALID
               exit
            endif
         enddo
      else
         total = INVALID
      endif

   end function total

end module setofprimes

module divisor

   use setofprimes, only: component_set

   implicit none

   private
   public :: factored_divisor

   type, extends(component_set) :: factored_divisor
      type(component_set) :: reference
   contains
      procedure :: reset
      procedure, private :: next
      procedure :: next_div
      procedure :: clear
      procedure :: is_valid
   end type factored_divisor

contains

   subroutine reset(this, f)

      implicit none

      class(factored_divisor), intent(inout) :: this  ! object invoking type-bound procedure
      class(component_set),    intent(in)    :: f     ! reference number, should be already factored

      call this%clear
      allocate(this%factors(size(f%factors)))
      allocate(this%reference%factors(size(f%factors)))

      this%factors(:)%prime = f%factors(:)%prime
      this%factors(:)%power = 0

      this%reference%factors = f%factors

   end subroutine reset

   subroutine next(this)

      use constants, only: INVALID

      implicit none

      class(factored_divisor), intent(inout) :: this  ! object invoking type-bound procedure

      integer :: i
      logical :: done

      done = .false.

      do i = lbound(this%factors, dim=1), ubound(this%factors, dim=1)
         associate (p => this%factors(i)%power)
            if (p < this%reference%factors(i)%power) then
               p = p + 1
               done = .true.
            else
               p = 0
            endif
         end associate
         if (done) exit
      enddo

      if (.not. done) this%factors(:)%power = INVALID

   end subroutine next

   subroutine next_div(this)

      use iso_fortran_env, only: INT64

      implicit none

      class(factored_divisor), intent(inout) :: this  ! object invoking type-bound procedure

      integer(kind=INT64) :: cur, nxt

      cur = this%total()
      nxt = huge(1_INT64)
      this%factors(:)%power = 0
      do while (this%is_valid())  ! for sure we can use more efficient searching or create a list and sort it but here it is fast enough for now
         call this%next
         if (this%is_valid()) then
            if (this%total() > cur .and. this%total() < nxt) nxt = this%total()
         endif
      enddo
      this%factors(:)%power = 0
      do while (this%is_valid())
!         write(*,*), cur, nxt, this%is_valid(), this%total()
         call this%next
         if (this%total() == nxt) exit
      enddo

   end subroutine next_div

   subroutine clear(this)

      implicit none

      class(factored_divisor), intent(inout) :: this  ! object invoking type-bound procedure

      if (allocated(this%factors)) deallocate(this%factors)
      if (allocated(this%reference%factors)) deallocate(this%reference%factors)

   end subroutine clear

   logical function is_valid(this)

      use constants, only: INVALID

      implicit none

      class(factored_divisor), intent(in) :: this  ! object invoking type-bound procedure

      is_valid = (allocated(this%factors) .and. allocated(this%reference%factors))
      if (is_valid) is_valid = all(this%factors(:)%power > INVALID)

   end function is_valid

end module divisor

module composition

   use constants,       only: buflen
   use iso_fortran_env, only: INT64
   use setofprimes,     only: component_set
   use singleprime,     only: component

   implicit none

   private
   public :: factorization_t

   type, extends(component_set) :: factorization_t
      integer(kind=INT64) :: number
      character(len=buflen) :: factor_str
   contains
      procedure :: factorize
      procedure :: erase
   end type factorization_t

contains

   subroutine factorize(this, n)

      use primes_utils, only: primes_t

      implicit none

      class(factorization_t), intent(inout) :: this  ! object invoking type-bound procedure
      integer(kind=INT64),    intent(in)    :: n     ! the number to factorize

      type(primes_t) :: primes
      integer :: i
      integer(kind=INT64) :: auxn
      character(len=buflen) :: n1, n2

      this%number = n
      ! Restrict the test to nicely divisible numbers to avoid computing excessively large list of primes.
      ! Calling primes%sieve(this%number) will be perfectly correct but would take long time for GB-sized buffer.
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

      this%factor_str = ""
      do i = lbound(this%factors, dim=1), ubound(this%factors, dim=1)
         write(n1, *) this%factors(i)%prime
         write(n2, *) this%factors(i)%power
         this%factor_str = adjustl(trim(this%factor_str) // " " // trim(adjustl(n1)) // "**" // trim(adjustl(n2)) // merge("  ", " *", i == ubound(this%factors, dim=1)))
      enddo

      if (this%total() /= this%number) then
         write(n1, *) this%number / this%total()
         this%factor_str = trim(this%factor_str) // " * " // trim(adjustl(n1)) // " (incomplete?)"
         write(*,*)"Decomposition failed: ", this%total(), " /= ", this%number
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
