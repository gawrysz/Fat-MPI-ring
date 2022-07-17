! Take care of the parallel environment.

module mpisetup

   use constants, only: INVALID
   use mpi,       only: MPI_INTEGER_KIND, MPI_ADDRESS_KIND, MPI_MAX_PROCESSOR_NAME

   implicit none

   private
   public :: parallel_init, parallel_finalize, &
        &    main_proc, nproc, proc, local_nproc, myname, ierr, tag_ub

   integer(kind=MPI_INTEGER_KIND), parameter :: main_proc = 0
   integer(kind=MPI_INTEGER_KIND), protected :: nproc = INVALID, proc = INVALID, local_nproc = INVALID, local_proc = INVALID
   integer(kind=MPI_INTEGER_KIND), protected :: local_comm
   integer(kind=MPI_INTEGER_KIND)            :: ierr
   integer(kind=MPI_ADDRESS_KIND), protected :: tag_ub = INVALID
   character(len=MPI_MAX_PROCESSOR_NAME), protected :: myname

contains

   ! Initialize the MPI,
   ! find out ranks and nodes,
   ! say "Hello".

   subroutine parallel_init

      use mpi, only: MPI_Init, MPI_Comm_size, MPI_Comm_rank, &
           &         MPI_COMM_WORLD, MPI_TAG_UB, MPI_COMM_TYPE_SHARED, MPI_INFO_NULL

      implicit none

      logical :: flag
      integer(kind=MPI_INTEGER_KIND) :: key, mynamelen

      call MPI_Init(ierr)

      call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, proc, ierr)
      call MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, tag_ub, flag, ierr)

      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, key, MPI_INFO_NULL, local_comm, ierr)
      call MPI_Comm_size(local_comm, local_nproc, ierr)
      call MPI_Comm_rank(local_comm, local_proc, ierr)
      call MPI_Get_processor_name(myname, mynamelen, ierr)

      if (.not. flag) write(*,*)"Problem with MPI_TAG_UB: ", tag_ub, " ???"

      write(*,'(2(a,i4),2a,2(a,i4))')"Hello from ", proc, "/", nproc, " @", trim(myname), ":", local_proc, "/", local_nproc

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if (proc == main_proc) write(*,'()')

   end subroutine parallel_init

   ! Clean up.

   subroutine parallel_finalize

      use mpi, only: MPI_Finalize, MPI_Barrier, MPI_COMM_WORLD

      implicit none

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_Finalize(ierr)

   end subroutine parallel_finalize

end module mpisetup
