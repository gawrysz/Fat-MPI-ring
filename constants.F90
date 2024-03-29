! various constants

module constants

   implicit none

   public

   integer, parameter :: INVALID = -1

   integer, parameter :: buflen = 512  ! length of buffers
   integer, parameter :: lablen = 32   ! length of labels

   ! code options
   enum, bind(C)
      enumerator :: V_NONE = 0, V_SPEED, V_STATS, V_DETAILED  ! verbosity levels
   end enum
   enum, bind(C)
      enumerator :: T_MPI_SR, T_MPI_GET1, T_MPI_PUT1, T_MPI_GETN, T_MPI_PUTN, T_MPI_SHMEM, T_COARRAY, T_OMP  ! test types to perform
   end enum

end module constants
