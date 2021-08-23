module mpisetup

   use constants, only: INVALID
   use mpi,       only: MPI_INTEGER_KIND

   implicit none

   private
   public :: parallel_init, parallel_finalize, &
        &    main_proc, nproc, proc, ierr

   integer(kind=MPI_INTEGER_KIND), parameter :: main_proc = 0
   integer(kind=MPI_INTEGER_KIND) :: nproc = INVALID, proc =INVALID, ierr

contains

   subroutine parallel_init

      use mpi, only: MPI_Init, MPI_Comm_size, MPI_Comm_rank, MPI_COMM_WORLD

      implicit none

      call MPI_Init(ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, proc, ierr)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

   end subroutine parallel_init

   subroutine parallel_finalize

      use mpi, only: MPI_Finalize, MPI_Barrier, MPI_COMM_WORLD

      implicit none

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_Finalize(ierr)

   end subroutine parallel_finalize

end module mpisetup
