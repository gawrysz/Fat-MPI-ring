! various constants

module constants

   implicit none

   public

   integer, parameter :: INT64 = selected_int_kind(16)
   integer, parameter :: FP64 = selected_real_kind(12)

   integer, parameter :: INVALID = -1

   integer, parameter :: buflen = 512

end module constants
