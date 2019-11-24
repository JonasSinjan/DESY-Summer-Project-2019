! NB NB NB NB NB
! Note on how FFTW stores the $k$ values in the arrays
!  first dimension: 0, 1, 2, ... n/2  (take integer division into account on last term)
!  second dimension - n even: 0, 1, 2, ... n/2, -n/2 + 1, - n/2 + 2, ... -1
!  second dimension - n odd:  0, 1, 2, ... n/2, -n/2, n/2 + 1, ... -1

program main

  !use spectrum, only: power_kk, dft_filter
  use, intrinsic :: iso_c_binding
  use omp_lib
  implicit none

  ! ------------------------------------------------------------------------
  ! include FFTW
  ! ------------------------------------------------------------------------
  !include 'fftw3.f03'

  ! ------------------------------------------------------------------------
  ! define fixed parameters
  ! ------------------------------------------------------------------------
#ifdef DP
  integer, parameter :: sp = kind(1.d0)
#else
  integer, parameter :: sp = kind(1.0)
#endif
  !double preicision, real number (with decimal) not integer
  integer, parameter :: dp = kind(1.d0)

  real(dp), parameter :: pi = acos(-1.d0)
  real(dp), parameter :: twopi = 2.d0*pi
  real(dp), allocatable :: wtime

  ! ------------------------------------------------------------------------
  ! define and initialize problem parameters
  ! ------------------------------------------------------------------------
  integer :: ngrids = 9
  real(sp) :: lx = twopi ! this the box space in fourier space?
  real(sp) :: ly = twopi
  real(sp) :: lz = twopi
  real(sp) :: bx0 = 1.
  real(sp) :: by0 = 0.
  real(sp) :: bz0 = 0. !3d

  ! ------------------------------------------------------------------------
  ! define variables
  ! ------------------------------------------------------------------------
  integer :: n, num_seed
  integer, allocatable :: seed(:)

  real(sp), dimension(:), allocatable :: x, y, z, amp_list
  real(sp), dimension(:,:,:), allocatable :: bx, by, bz !3d
  real(sp) :: time, dx, dy, dz


  ! ------------------------------------------------------------------------
  ! define dummy variables
  ! ------------------------------------------------------------------------
  integer :: i, j, k
  integer :: lun

  character(len=400) :: data_dir
  character(len=1024) :: file_out
  character(len=400) :: cmd
  
  time = 0.

  wtime = omp_get_wtime()
  ! ------------------------------------------------------------------------
  ! specify folder for output data
  ! ------------------------------------------------------------------------
  data_dir = './Runs/512_B_amp1/'

  cmd = 'mkdir -p ' // trim(data_dir)
  call system(cmd)

  ! ------------------------------------------------------------------------
  ! calculate grid parameters
  ! ------------------------------------------------------------------------
  n = 2**(ngrids) + 1

  call random_seed(size = num_seed)
  allocate(seed(num_seed))

  ! creating random seed
  do i = 1, num_seed
    seed(i) = i*4251
  enddo

  call random_seed(put=seed)

! ------------------------------------------------------------------------
  ! allocate arrays
  ! ------------------------------------------------------------------------

  allocate(x(n))
  allocate(y(n))
  allocate(z(n))!3d

  ! ------------------------------------------------------------------------
  ! set space grid and time of output files
  ! ------------------------------------------------------------------------
  dx = lx/real(n)
  dy = ly/real(n)
  dz = lz/real(n) !3d

  do i = 1, n
    x(i) = (real(i) - 0.5)*dx
  enddo

  do j = 1, n
    y(j) = (real(j) - 0.5)*dy
  enddo

  do k = 1, n !3d n instead of 1,1
    z(k) = (real(k) - 0.5)*dz
  enddo


  ! ------------------------------------------------------------------------
  ! build bx, by and bz and write to file
  ! ------------------------------------------------------------------------
  allocate (bx(n,n,n)) !n,n,n for all?
  allocate (by(n,n,n))
  allocate (bz(n,n,n)) !3d
  allocate (amp_list(24))
  !-------------------------------------------------------------------
  !3 arrays allocated
  !-------------------------------------------------------------------

  bx(:,:,:) = bx0 ! :,:,:? for 3d
  by(:,:,:) = by0
  bz(:,:,:) = bz0 !3d

  !do I need to vary in 3rd direction now too? k?
  amp_list = (/1.0,2.3,1.8,0.2,0.5,1.6,1.7,1.1,0.1,1.2,2.1,1.5,0.3,0.9,1.4,0.7,1.2,1.3,2.3,0.4,0.8,1.8,1.9,0.6/)
  ! M_A = 4.9251199812242934
  !amp_list = amp_list*0.5
  ! M_A = 2.4625599906121467
  do k = 1, n
    do j = 1, n
      do i = 1, n

        by(i,j,k) = by(i,j,k) + amp_list(1)*sin(3.0*x(i) - 8.0*z(k))
        by(i,j,k) = by(i,j,k) + amp_list(2)*sin(8.0*x(i) + 1.6)
        by(i,j,k) = by(i,j,k) + amp_list(3)*sin(7.0*x(i) - 7.0*z(k))
        by(i,j,k) = by(i,j,k) + amp_list(4)*sin(4.0*x(i) + 1.6)
        by(i,j,k) = by(i,j,k) + amp_list(5)*sin(2.0*x(i) - 4.0*z(k))
        by(i,j,k) = by(i,j,k) + amp_list(6)*sin(1.0*x(i) + 1.6)
        by(i,j,k) = by(i,j,k) + amp_list(7)*sin(2.0*x(i) - 3.0*z(k))
        by(i,j,k) = by(i,j,k) + amp_list(8)*sin(3.0*x(i) + 1.6)

        bx(i,j,k) = bx(i,j,k) + amp_list(9)*cos(6.0*y(j) + 6.0*z(k))
        bx(i,j,k) = bx(i,j,k) + amp_list(10)*cos(7.0*y(j)+ 1.6 + 5.0*z(i))
        bx(i,j,k) = bx(i,j,k) + amp_list(11)*cos(3.0*y(j) + 6*z(k))
        bx(i,j,k) = bx(i,j,k) + amp_list(12)*cos(4.0*y(j)+ 1.6 + 3.0*z(i))
        bx(i,j,k) = bx(i,j,k) + amp_list(13)*cos(5.0*y(j) + 2.0*z(k))
        bx(i,j,k) = bx(i,j,k) + amp_list(14)*cos(2.0*y(j) + 4.0*z(i))
        bx(i,j,k) = bx(i,j,k) + amp_list(15)*cos(9.0*y(j) + 2.0*z(k))
        bx(i,j,k) = bx(i,j,k) + amp_list(16)*cos(y(j)+ 1.6 + 1.0*z(i))
   
        bz(i,j,k) = bz(i,j,k) + amp_list(17)*cos(3.0*x(i) + 3.0*y(j))
        bz(i,j,k) = bz(i,j,k) + amp_list(18)*sin(2.0*x(i)- 5.0*y(j))
        bz(i,j,k) = bz(i,j,k) + amp_list(19)*cos(3.0*x(i) + 1.6 + 3.0*y(j))
        bz(i,j,k) = bz(i,j,k) + amp_list(20)*sin(2.0*x(i)- 2.0*y(j))
        bz(i,j,k) = bz(i,j,k) + amp_list(21)*cos(3.0*x(i) + 7.0*y(j))
        bz(i,j,k) = bz(i,j,k) + amp_list(22)*sin(2.0*x(i) + 1.6 - 10.0*y(j))
        bz(i,j,k) = bz(i,j,k) + amp_list(23)*cos(3.0*x(i) + 3.0*y(j))
        bz(i,j,k) = bz(i,j,k) + amp_list(24)*sin(2.0*x(i)- 9.0*y(j))

      enddo
    enddo
  enddo

  lun = 701
  file_out = trim(data_dir) // '/' // 'BX.BIN'
  ! bx(:,:,:)?
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) bx(:,:,:)
  close(lun)

  lun = 701
  file_out = trim(data_dir) // '/' // 'BY.BIN'
  ! by(:,:,:)?
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) by(:,:,:)
  close(lun)

  !3d
  lun = 701
  file_out = trim(data_dir) // '/' // 'BZ.BIN'
  ! bz(:,:,:)?
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) bz(:,:,:)
  close(lun)
 
  wtime = omp_get_wtime() - wtime

  print *, "Time taken = ", wtime

end program main
