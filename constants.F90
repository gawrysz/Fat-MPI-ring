! various constants

module constants

   implicit none

   public

   ! Not needed, provided by iso_fortran_env
   ! integer, parameter :: I64 = selected_int_kind(16)
   ! integer, parameter :: FP64 = selected_real_kind(12)

   integer, parameter :: INVALID = -1

   integer, parameter :: buflen = 512

   ! code options
   enum, bind(C)
      enumerator :: V_NONE = 0, V_SPEED, V_STATS, V_DETAILED  ! verbosity levels
   end enum
   enum, bind(C)
      enumerator :: T_MPI_SR, T_MPI_G, T_MPI_P, T_MPI_SHMEM, T_COARRAY, T_OMP  ! test types to perform
   end enum

end module constants
