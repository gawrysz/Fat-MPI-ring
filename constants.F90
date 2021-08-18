! various constants

module constants

   implicit none

   public

   integer, parameter :: INT64 = selected_int_kind(16)
   integer, parameter :: FP64 = selected_real_kind(12)

   integer, parameter :: INVALID = -1

   integer, parameter :: buflen = 512

   ! code options
   enum, bind(C)
      enumerator :: V_NONE = 0, V_SPEED, V_STATS, V_DETAILED  ! verbosity levels
   end enum
   enum, bind(C)
      enumerator :: T_MPI = 1, T_COARRAY = 2 * T_MPI, T_OMP = 2 * T_COARRAY  ! masks for test types to perform
   end enum

end module constants
