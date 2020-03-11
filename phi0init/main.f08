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
  include 'fftw3.f03'

  integer :: start_time, stop_time !timing fortran program
  integer :: count_rate, count_max

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
  real(dp), allocatable :: wtime, tot_time
  real(dp), allocatable :: time_init, time_final, elapsed_time
  integer, allocatable :: threadno
  integer :: thread_id

  ! ------------------------------------------------------------------------
  ! define and initialize problem parameters
  ! ------------------------------------------------------------------------
  integer :: ngrids = 9
  real(sp) :: anis = 1.
  real(sp) :: lx = twopi ! this the box space in fourier space?
  real(sp) :: ly = twopi
  real(sp) :: lz = twopi

  ! ------------------------------------------------------------------------
  ! define variables
  ! ------------------------------------------------------------------------
  integer :: n, num_seed
  integer, dimension(1) :: rand_seed
  integer, allocatable :: seed(:)
  real(sp) :: num
  real(sp) :: h
  real(sp) :: amp
  real(sp) :: kx, ky, kz

  real(sp), dimension(:), allocatable :: x, y, z

  real(sp), dimension(:,:,:), allocatable :: phi0
  complex(sp), dimension(:,:,:), allocatable :: phi0k,fk

  integer :: nk, nkb, nkt

  logical :: lsum_power

  ! ------------------------------------------------------------------------
  ! define auxiliary variables
  ! ------------------------------------------------------------------------
  type(C_PTR) :: plan_phi0
  type(C_PTR) :: plan1, plan2
  type(C_PTR) :: plan
  integer*8 :: dftplan
  real(sp), dimension(:,:,:), allocatable :: f

  real(sp) :: kmax, kmod, mag
  real(sp) :: test_var
  real(sp) :: k_para, k_perp, E_coeff, ph


  ! ------------------------------------------------------------------------
  ! define dummy variables
  ! ------------------------------------------------------------------------
  integer :: l,q
  integer :: m
  integer :: i, j, k
  integer :: ii, jj, k_k
  integer :: ki, kj, kk !3d
  integer :: lun

  character(len=400) :: data_dir
  character(len=1024) :: file_out
  character(len=400) :: cmd


  ! ------------------------------------------------------------------------
  ! specify folder for output data
  ! ------------------------------------------------------------------------
  data_dir = './Runs/512_15_kpara/'

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
  ! build scalar field phi0 and write to file
  ! ------------------------------------------------------------------------
  m = n - 1

  allocate (phi0(n,n,n)) 
  !aux arrays
  allocate (phi0k((m/2 + 1), m, m))
  allocate (fk((m/2 + 1), m, m))
  allocate (f(m, m, m))

  

  phi0(:,:,:) = 0 ! initialise all entries to zero
  
  print*, omp_get_max_threads(), "calculating phi0k"
  wtime = omp_get_wtime()

  call dfftw_plan_dft_c2r_3d(dftplan, m,m,m, fk, f, FFTW_ESTIMATE)

  phi0k(:,:,:) = 0

  kmax = real(m)/2 !? was /2 for 2D - at high wavenumber not adding turbulent modes

  ! SKIP NYQUIST FREQUENCY
  do kk = min((-m/2 + 1), 0), m/2-1
    ! do kj = 0, 0
    if (kk >= 0) then
      k = kk + 1
    else
      k = m + kk + 1
    endif  
    
    !print*, kk

    ! SKIP NYQUIST FREQUENCY
    do kj = min((-m/2 + 1), 0), m/2-1
      ! do kj = 0, 0
      if (kj >= 0) then
        j = kj + 1
      else
        j = m + kj + 1
      endif

      ! SKIP NYQUIST FREQUENCY
      do ki = 0, m/2-1
        i = ki + 1
        
        if ((ki == 0) .and. (kj == 0) .and. (kk == 0)) cycle

        kmod = sqrt(real(ki)**2 + real(kj)**2 + real(kk)**2)

        k_para = abs(ki)
        k_perp = sqrt(max((kmod**2 - k_para**2), 0.))

        ! GS95
        if (k_perp > 0.) then
          E_coeff = k_perp**(-10./3.)*exp(-k_para/k_perp**(2./3.))  ! 3D
        else
          E_coeff = 15*k_para**(-2) !different amplitude?
          !E_coeff = 0
        endif

        !sort random phase
        call random_number(ph)
        ph = ph*twopi
        ! if (kmod > kmax) then from michael's method
        !   phi0k(i,j,k) = (0.d0, 0.d0)
        ! else
        !   phi0k(i,j,k) = sqrt(E_coeff)*(cos(ph) + (0., 1.)*sin(ph))
        ! endif
        !/home/jonas/Documents/VSCode/DESY/phi0init/Runs/512_test
        ! reality condition of FFT - conjugate
        ! sort random phase
        call random_number(ph)
        ph = ph*twopi
        if (ki == 0) then
          if (kj > 0) then
            phi0k(i,j,k) = sqrt(E_coeff)*(cos(ph) + (0., 1.)*sin(ph))
            if (k/=1) then
              phi0k(i,m-j+2,m-k+2) = sqrt(E_coeff)*(cos(ph) - (0., 1.)*sin(ph))
            else 
              phi0k(i,m-j+2,k) = sqrt(E_coeff)*(cos(ph) - (0., 1.)*sin(ph))
            endif
          else if (kj < 0) then
            cycle
          else
            if (kk > 0) then
              phi0k(i,j,k) = sqrt(E_coeff)*(cos(ph) + (0., 1.)*sin(ph))
              phi0k(i,j,m-k+2) = sqrt(E_coeff)*(cos(ph) - (0., 1.)*sin(ph))
            else
              cycle
            endif
          endif  
        else
          phi0k(i,j,k) = sqrt(E_coeff)*(cos(ph) + (0., 1.)*sin(ph))
        endif
      enddo ! ki
    enddo ! kj
  enddo !kk
  
  ! execute inverse DFT
  ! attention with the normalization of the DFT
  fk(:,:,:) = phi0k(:,:,:)

  call dfftw_execute_dft_c2r(dftplan, fk, f)
  phi0(1:m,1:m,1:m) = f(:,:,:) !inner cube

  phi0(m+1,1:m,1:m) = f(1,:,:) !outer faces
  phi0(1:m,m+1,1:m) = f(:,1,:)
  phi0(1:m,1:m,m+1) = f(:,:,1)

  phi0(1:m,m+1,m+1) = f(:,1,1) !outer edges
  phi0(m+1,m+1,1:m) = f(1,1,:)
  phi0(m+1,1:m,m+1) = f(1,:,1)

  phi0(m+1,m+1,m+1) = f(1,1,1) !last point

  call dfftw_destroy_plan(dftplan)

  deallocate (phi0k)
  deallocate (fk)
  deallocate (f)
  
  print*, 'The loop has successfully completed'

  wtime = omp_get_wtime() - wtime

  print *, wtime, "Loop Time Phi0 init" 

  print*, '* Writing file: phi_0'

  file_out = 'PHI0.BIN'
  lun = 701
  file_out = trim(data_dir) // '/' // file_out
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) phi0(:,:,:)
  close(lun)

  
  deallocate (phi0)

  stop

end program main
