module memory

   use constants, only: INVALID, INT64

   implicit none

   private
   public :: memcheck, totmem

   integer(kind=INT64), protected :: totmem = INVALID

contains

   logical function memcheck()

      use constants, only: buflen
      use mpisetup,  only: local_nproc

      implicit none

      integer, parameter :: pidlen = 8
      character(len=buflen) :: filename, line
      character(len=pidlen) :: pid_char
      character(len=7) :: vm = "VmData:"
      integer :: stat_lun, io, system_mem_usage
      real, parameter :: warnlevel = 0.5

      if (totmem <= INVALID) call find_totmem
      memcheck = (totmem > INVALID)

      write(pid_char, '(i8)') getpid()
      filename = '/proc/' // trim(adjustl(pid_char)) // '/status'

      if (.not. file_exists(filename)) then
         memcheck = .false.
         return
      endif

      open(newunit=stat_lun, file=filename, status='old')
      do
         read (stat_lun,'(a)', iostat=io) line
         if (io /= 0) exit
         if (line(1:7) == vm) then
            read (line(8:), *) system_mem_usage
            exit
         endif
      enddo
      close(stat_lun)

      if (real(system_mem_usage) * local_nproc > warnlevel * real(totmem)) then
         memcheck = .false.
         write(*, '(a,i10,a,f5.1,a)')"# " // vm // " = ", system_mem_usage, ", ", 100. * real(system_mem_usage)/real(totmem), "% of RAM. Too much."
      endif

   end function memcheck

   subroutine find_totmem

      use constants, only: buflen

      implicit none

      character(len=buflen) :: filename, line
      integer :: stat_lun, io

      ! grep MemTotal: /proc/meminfo | awk '{print $2}'
      filename = '/proc/meminfo'

      if (.not. file_exists(filename)) return

      open(newunit=stat_lun, file=filename, status='old')
      do
         read (stat_lun,'(a)', iostat=io) line
         if (io /= 0) exit
         if (line(1:9) == 'MemTotal:') then
            read (line(10:), *) totmem
            exit
         endif
      enddo
      close(stat_lun)

   end subroutine find_totmem

   logical function file_exists(fname)

      implicit none

      character(len=*), intent(in) :: fname

      inquire (file=fname, exist=file_exists)
      if (.not. file_exists) write(*, '(3a)') "# Cannot find '", trim(fname), "'"

   end function file_exists

end module memory

