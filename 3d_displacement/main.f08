! NB NB NB NB NB
! Note on how FFTW stores the $k$ values in the arrays
!  first dimension: 0, 1, 2, ... n/2  (take integer division into account on last term)
!  second dimension - n even: 0, 1, 2, ... n/2, -n/2 + 1, - n/2 + 2, ... -1
!  second dimension - n odd:  0, 1, 2, ... n/2, -n/2, n/2 + 1, ... -1

program main

  use spectrum, only: power_kk, dft_filter
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
  real(dp), allocatable :: wtime
  real(dp), allocatable :: time_init, time_final, elapsed_time
  integer, allocatable :: threadno
  integer :: thread_id
  ! ------------------------------------------------------------------------
  ! define types
  ! ------------------------------------------------------------------------
  type :: sgrid
    real(sp), dimension(:,:,:), allocatable :: bx, by, bz !3d IS BY JUST THE PERPENDICULAR OR DO WE NEED BX
    real(sp), dimension(:,:,:), allocatable :: dbx, dby, dbz
    real(sp), dimension(:,:,:), allocatable :: etz_x,etz_y, etz_z !3d
    real(sp), dimension(:,:,:), allocatable :: drx, dry, drz
  end type sgrid


  ! ------------------------------------------------------------------------
  ! define and initialize problem parameters
  ! ------------------------------------------------------------------------
  integer :: ngrids = 6
  real(sp) :: bx0 = 1.
  real(sp) :: by0 = 0.
  real(sp) :: bz0 = 0. !3d
  real(sp) :: anis = 1.
  real(sp) :: lx = twopi ! this the box space in fourier space?
  real(sp) :: ly = twopi
  real(sp) :: lz = twopi

  ! ------------------------------------------------------------------------
  ! define variables
  ! ------------------------------------------------------------------------
  integer :: n
  integer, dimension(1) :: rand_seed
  real(sp) :: num
  real(sp) :: h
  real(sp) :: amp
  real(sp) :: kx, ky, kz
  real(sp), dimension(:,:,:), allocatable :: bx, by, bz !3d
  real(sp), dimension(:,:,:), allocatable :: rx0, ry0, rz0 !3d
  real(sp), dimension(:,:,:), allocatable :: drx, dry, drz !3d
  real(sp), dimension(:,:,:), allocatable :: x_arr, y_arr, z_arr, input_1
  type(sgrid), dimension(:), allocatable :: mgrid

  real(sp), dimension(:), allocatable :: x, y, z
  real(sp) :: time, dx, dy, dz

  real(sp), dimension(:,:,:), allocatable :: phi0
  real(sp), dimension(:,:,:), allocatable :: phi
  complex(sp), dimension(:,:,:), allocatable :: phi0k

  integer :: nk, nkb, nkt
  real(dp), dimension(:), allocatable :: ps_k !not sure which one this is used for
  real(dp), dimension(:), allocatable :: ps_kb !parallel, only 1D along x
  real(dp), dimension(:,:), allocatable :: ps_kt !perpendicular
  real(dp), dimension(:,:,:), allocatable :: ps_kk !all components
  real(dp), dimension(:,:,:), allocatable :: tmp3d

  logical :: lsum_power

  ! ------------------------------------------------------------------------
  ! define auxiliary variables
  ! ------------------------------------------------------------------------
  type(C_PTR) :: plan_phi0
  type(C_PTR) :: plan1, plan2
  real(sp), dimension(:,:,:), allocatable :: f
  complex(sp), dimension(:,:,:), allocatable :: bxk, byk, bzk !3d
  complex(sp), dimension(:,:,:), allocatable :: dbxk, dbyk, dbzk !3d
  complex(sp), dimension(:,:,:), allocatable :: etzk_x, etzk_y, etzk_z

  real(sp) :: b(3), b2 !3d
  real(sp) :: b1(3), ph1, ph2, ph3, aux1, aux2, aux3 !3d
  real(sp) :: rxp, ryp, drxp, dryp, rzp, drzp !3d
  real(sp) :: wx0, wx1, wy0, wy1, wz0, wz1
  real(sp) :: kmax, kmod, mag

  real(sp) :: k_para, k_perp, E_coeff, ph


  ! ------------------------------------------------------------------------
  ! define dummy variables
  ! ------------------------------------------------------------------------
  integer :: l
  integer :: m
  integer :: i, j, k
  integer :: ii, jj, k_k
  integer :: iip1, jjp1, kkp1
  integer :: ki, kj, kk !3d
  integer :: lun
  real(sp) :: tmp, tmp2
  real(sp) :: c_00,c_01,c_10,c_11,c_0,c_1,c !for the trilinear interpolations

  character(len=400) :: data_dir
  character(len=1024) :: file_out
  character(len=400) :: cmd


  ! ------------------------------------------------------------------------
  ! specify folder for output data
  ! ------------------------------------------------------------------------
  data_dir = './64run3D_OpenMP/'

  cmd = 'mkdir -p ' // trim(data_dir)
  call system(cmd)

  cmd = 'mkdir -p ' // trim(data_dir) // 'power_spectra/'
  call system(cmd)

  cmd = 'mkdir -p ' // trim(data_dir) // 'structure_function/'
  call system(cmd)

  ! ------------------------------------------------------------------------
  ! calculate grid parameters
  ! ------------------------------------------------------------------------
  n = 2**(ngrids) + 1
  h = twopi/real(n - 1) !box length?


  nk = n/2 + 1
  nkb = n/2 + 1
  nkt = n/2 + 1


  ! ------------------------------------------------------------------------
  ! allocate arrays
  ! ------------------------------------------------------------------------
  allocate (tmp3d(n,n,n))

  allocate (ps_k(nk))
  allocate (ps_kb(nkb))
  allocate (ps_kt(nkt,nkt))
  allocate (ps_kk(nkb,nkt,nkt))

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

  time = 0.


  ! ------------------------------------------------------------------------
  ! build bx, by and bz and write to file
  ! ------------------------------------------------------------------------
  allocate (bx(n,n,n)) !n,n,n for all?
  allocate (by(n,n,n))
  allocate (bz(n,n,n)) !3d

  bx(:,:,:) = bx0 ! :,:,:? for 3d
  by(:,:,:) = by0
  bz(:,:,:) = bz0 !3d

  !do I need to vary in 3rd direction now too? k?
  do k = 1, n
    do j = 1, n
      do i = 1, n
        by(i,j,k) = by(i,j,k) + 0.5*sin(2.0*x(i))
        by(i,j,k) = by(i,j,k) + 0.5*sin(4.0*x(i)+1.6)
      enddo
    enddo
  enddo

  file_out = trim(data_dir) // 'MAGNETIC.DAT' 
  open(unit=400, file=trim(file_out), form='formatted', status='replace', action='write')
  do i = 1, n
    do j = 1, n
      do k = 1, n !3d
        write(400,'(4(es24.16, 1x))') real(i - 1)*h, real(j - 1)*h, bx(i,j,k), by(i,j,k)
      enddo
    enddo
    write(400,*)
  enddo
  close(400)

  lun = 701
  file_out = trim(data_dir) // '/' // 'X.BIN'
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) x(:)
  close(lun)

  lun = 701
  file_out = trim(data_dir) // '/' // 'Y.BIN'
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) y(:)
  close(lun)

  lun = 701 !3D
  file_out = trim(data_dir) // '/' // 'Z.BIN'
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) z(:)
  close(lun)

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
  ! bZ(:,:,:)?
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) bZ(:,:,:)
  close(lun)


  ! ------------------------------------------------------------------------
  ! generate grids hierarchy and allocate memory for each grid level
  ! ------------------------------------------------------------------------
  allocate (mgrid(ngrids))

  do l = 1, ngrids
    allocate (mgrid(l)%bx(n,n,n)) ! n,n,n ?
    allocate (mgrid(l)%by(n,n,n))
    allocate (mgrid(l)%bz(n,n,n)) !3d
    allocate (mgrid(l)%dbx(n,n,n))
    allocate (mgrid(l)%dby(n,n,n))
    allocate (mgrid(l)%dbz(n,n,n)) !3d
    allocate (mgrid(l)%etz_x(n,n,n))
    allocate (mgrid(l)%etz_y(n,n,n))
    allocate (mgrid(l)%etz_z(n,n,n))
    allocate (mgrid(l)%drx(n,n,n))
    allocate (mgrid(l)%dry(n,n,n))
    allocate (mgrid(l)%drz(n,n,n)) !3d
  enddo


  ! ------------------------------------------------------------------------
  ! build bx and by in each grid level, using spectral filtering
  ! ------------------------------------------------------------------------
  mgrid(1)%bx(:,:,:) = bx(:,:,:) !:,:,: 3d
  mgrid(1)%by(:,:,:) = by(:,:,:)
  mgrid(1)%bz(:,:,:) = bz(:,:,:) !3d

  do l = 2, ngrids

    kmax = sqrt(2.d0)*(2.d0**(ngrids - l))

    m = n - 1

    ! allocate auxiliary matrices
    allocate (f(m, m, m)) !3d
    allocate (bxk((m/2 + 1), m, m)) !3d added extra m dimension - split first dimension in half
    allocate (byk((m/2 + 1), m, m))
    allocate (bzk((m/2 + 1), m, m)) 

    ! prepare plans for the dft (plan1) and dft inverse (plan2)
#ifdef DP
    plan1 = fftw_plan_dft_r2c_3d(m, m, m, f, bxk, FFTW_ESTIMATE) !3D
    plan2 = fftw_plan_dft_c2r_3d(m, m, m, bxk, f, FFTW_ESTIMATE)
#else
    plan1 = fftwf_plan_dft_r2c_3d(m, m, m, f, bxk, FFTW_ESTIMATE)
    plan2 = fftwf_plan_dft_c2r_3d(m, m, m, bxk, f, FFTW_ESTIMATE)
#endif

    ! transform mgrid(l-1)%bx and mgrid(l-1)%by to fourier space, destroy plan1

#ifdef DP
    f(:,:,:) = mgrid(l-1)%bx(1:m,1:m,1:m) !1:m for all dimensions?
    call fftw_execute_dft_r2c(plan1, f, bxk)
    !r2c = real to complex
    f(:,:,:) = mgrid(l-1)%by(1:m,1:m,1:m)
    call fftw_execute_dft_r2c(plan1, f, byk)

    f(:,:,:) = mgrid(l-1)%bz(1:m,1:m,1:m) !3d again do i need bzk?
    call fftw_execute_dft_r2c(plan1, f, bzk)

    call fftw_destroy_plan(plan1)
#else
    f(:,:,:) = mgrid(l-1)%bx(1:m,1:m,1:m)
    call fftwf_execute_dft_r2c(plan1, f, bxk)

    f(:,:,:) = mgrid(l-1)%by(1:m,1:m,1:m)
    call fftwf_execute_dft_r2c(plan1, f, byk)

    f(:,:,:) = mgrid(l-1)%bz(1:m,1:m,1:m)
    call fftwf_execute_dft_r2c(plan1, f, bzk)

    call fftwf_destroy_plan(plan1)
#endif

    ! normalize bxk and byk & bzk???
    bxk(:,:,:) = bxk(:,:,:) / real(m*m*m) !3d
    byk(:,:,:) = byk(:,:,:) / real(m*m*m)
    bzk(:,:,:) = bzk(:,:,:) / real(m*m*m)

    ! loop over all the modes
    ! (i, j) are the indices of the matrices generated by the dft
    ! corresponding to the wave numbers (ki, kj)
    do kk = min((-m/2 + 1), 0), m/2
      if (kk >= 0) then
        k = kk + 1
      else
        k = m + kk + 1
      endif  
      do kj = min((-m/2 + 1), 0), m/2
        if (kj >= 0) then
          j = kj + 1
        else
          j = m + kj + 1
        endif

        do ki = 0, m/2
          i = ki + 1

          kmod = sqrt(real(ki)**2 + real(kj)**2 + real(kk)**2)

          if (kmod > kmax) then
            bxk(i,j,k) = (0., 0.)
            byk(i,j,k) = (0., 0.)
            bzk(i,j,k) = (0., 0.)
          endif

        enddo
      enddo
    enddo

    ! transform bxk, byk, bzk to real space, apply periodic boundaries, destroy plan2
#ifdef DP
    call fftw_execute_dft_c2r(plan2, bxk, f)
    !c2r = complex to real
    !need to add the slices for periodic boundaries
    mgrid(l)%bx(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%bx(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%bx(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%bx(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%bx(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%bx(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%bx(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%bx(m+1,m+1,m+1) = f(1,1,1) !last point

    call fftw_execute_dft_c2r(plan2, byk, f)
    mgrid(l)%by(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%by(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%by(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%by(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%by(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%by(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%by(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%by(m+1,m+1,m+1) = f(1,1,1) !last point

    call fftw_execute_dft_c2r(plan2, bzk, f) !3d
    mgrid(l)%bz(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%bz(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%bz(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%bz(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%bz(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%bz(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%bz(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%bz(m+1,m+1,m+1) = f(1,1,1) !last point

    call fftw_destroy_plan(plan2)
#else
    call fftw_execute_dft_c2r(plan2, bxk, f)
    !c2r = complex to real
    !need to add the slices for periodic boundaries
    mgrid(l)%bx(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%bx(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%bx(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%bx(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%bx(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%bx(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%bx(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%bx(m+1,m+1,m+1) = f(1,1,1) !last point

    call fftw_execute_dft_c2r(plan2, byk, f)
    mgrid(l)%by(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%by(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%by(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%by(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%by(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%by(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%by(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%by(m+1,m+1,m+1) = f(1,1,1) !last point

    call fftw_execute_dft_c2r(plan2, bzk, f) !3d
    mgrid(l)%bz(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%bz(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%bz(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%bz(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%bz(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%bz(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%bz(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%bz(m+1,m+1,m+1) = f(1,1,1) !last point
    
    call fftw_destroy_plan(plan2)

#endif

    ! deallocate auxiliary matrices
    deallocate (f)
    deallocate (bxk)
    deallocate (byk)
    deallocate (bzk)

    ! mgrid(l)%bx(:,:) = bx0
    ! mgrid(l)%by(:,:) = by0

  enddo


  ! ------------------------------------------------------------------------
  ! write bx and by in each grid to a different file
  ! ------------------------------------------------------------------------
  do l = 1, ngrids

    write (file_out, "('BX_GRID', i0, '.BIN')") l
    lun = 701
    file_out = trim(data_dir) // '/' // file_out
    open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
      write(lun) mgrid(l)%bx(:,:,:)
    close(lun)

    write (file_out, "('BY_GRID', i0, '.BIN')") l
    !why 702 all of a sudden??
    lun = 702
    file_out = trim(data_dir) // '/' // file_out
    open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
      write(lun) mgrid(l)%by(:,:,:)
    close(lun)

    write (file_out, "('BZ_GRID', i0, '.BIN')") l
    !why 702 all of a sudden??
    lun = 702
    file_out = trim(data_dir) // '/' // file_out
    open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
      write(lun) mgrid(l)%bz(:,:,:) !3d
    close(lun)

  enddo


  ! ------------------------------------------------------------------------
  ! calculate vector field db in each grid level
  ! ------------------------------------------------------------------------
  do l = 1, ngrids
    do k = 1, n !3d
      do j = 1, n
        do i = 1, n

          if (l < ngrids) then
            b(1) = mgrid(l+1)%bx(i,j,k)
            b(2) = mgrid(l+1)%by(i,j,k)
            b(3) = mgrid(l+1)%bz(i,j,k)
          else
            ! use average uniform field
            b(1) = bx0
            b(2) = by0
            b(3) = bz0
          endif

          mgrid(l)%dbx(i,j,k) = b(1) - mgrid(l)%bx(i,j,k)
          mgrid(l)%dby(i,j,k) = b(2) - mgrid(l)%by(i,j,k)
          mgrid(l)%dbz(i,j,k) = b(3) - mgrid(l)%bz(i,j,k) !3d

        enddo
      enddo
    enddo
  enddo


  ! ------------------------------------------------------------------------
  ! calculate etz (only non-null component) from equation 
  ! curl(et) = db in each grid level
  ! ------------------------------------------------------------------------
  do l = 1, ngrids

    m = n - 1

    ! allocate auxiliary matrices
    allocate (f(m, m, m)) !3d
    allocate (dbxk((m/2 + 1), m, m))
    allocate (dbyk((m/2 + 1), m, m))
    allocate (dbzk((m/2 + 1), m, m))
    allocate (etzk_x((m/2 + 1), m, m))
    allocate (etzk_y((m/2 + 1), m, m))
    allocate (etzk_z((m/2 + 1), m, m))

    ! prepare plans for the dft (plan1) and dft inverse (plan2)
#ifdef DP
    plan1 = fftw_plan_dft_r2c_3d(m, m, m, f, dbxk, FFTW_ESTIMATE) !3d
    plan2 = fftw_plan_dft_c2r_3d(m, m, m, etzk_x, f, FFTW_ESTIMATE) !3d do i need plans for each component?
#else
    plan1 = fftwf_plan_dft_r2c_3d(m, m, m, f, dbxk, FFTW_ESTIMATE)
    plan2 = fftwf_plan_dft_c2r_3d(m, m, m, etzk_x, f, FFTW_ESTIMATE) !3d do i need plans for each component?
#endif

    ! transform dbx and dby to fourier space, destroy plan1
#ifdef DP
    f(:,:,:) = mgrid(l)%dbx(1:m,1:m,1:m)
    call fftw_execute_dft_r2c(plan1, f, dbxk)

    f(:,:,:) = mgrid(l)%dby(1:m,1:m,1:m)
    call fftw_execute_dft_r2c(plan1, f, dbyk)

    f(:,:,:) = mgrid(l)%dbz(1:m,1:m,1:m) !3d
    call fftw_execute_dft_r2c(plan1, f, dbzk)

    call fftw_destroy_plan(plan1)
#else
    f(:,:,:) = mgrid(l)%dbx(1:m,1:m,1:m)
    call fftwf_execute_dft_r2c(plan1, f, dbxk)

    f(:,:,:) = mgrid(l)%dby(1:m,1:m,1:m)
    call fftwf_execute_dft_r2c(plan1, f, dbyk)

    f(:,:,:) = mgrid(l)%dbz(1:m,1:m,1:m) !3d
    call fftw_execute_dft_r2c(plan1, f, dbzk)

    call fftwf_destroy_plan(plan1)
#endif

    ! normalize dbxk and dbyk
    dbxk(:,:,:) = dbxk(:,:,:) / real(m*m*m)
    dbyk(:,:,:) = dbyk(:,:,:) / real(m*m*m)
    dbzk(:,:,:) = dbzk(:,:,:) / real(m*m*m)

    ! loop over all the modes
    ! (i, j) are the indices of the matrices generated by the dft
    ! corresponding to the wave numbers (ki, kj)
    
    !this unsure about
    
    do kk = min((-m/2 + 1), 0), m/2
      if (kk >= 0) then
        k = kk + 1
      else
        k = m + kk + 1
      endif    
      do kj = min((-m/2 + 1), 0), m/2
        if (kj >= 0) then
          j = kj + 1
        else
          j = m + kj + 1
        endif

        do ki = 0, m/2
          i = ki + 1

          if ((ki==0) .and. (kj==0) .and. (kk==0)) then
            etzk_x(i,j,k) = (0., 0.)
            etzk_y(i,j,k) = (0., 0.)
            etzk_z(i,j,k) = (0., 0.)
            cycle
          endif

          ! etk = - i (dbk x k) / k^2
          !in 3d etzk will still be a vector
          mag = (real(ki)**2 + real(kj)**2+real(kk)**2) !k**2 magnitude

          etzk_x(i,j,k) = - (0., 1.)*(dbyk(i,j,k)*real(kk) - dbzk(i,j,k)*real(kj))/mag

          etzk_y(i,j,k) = - (0., 1.)*(dbzk(i,j,k)*real(ki) - dbxk(i,j,k)*real(kk))/mag

          etzk_z(i,j,k) = - (0., 1.)*(dbxk(i,j,k)*real(kj) - dbyk(i,j,k)*real(ki))/mag

        enddo
      enddo
    enddo

    ! transform each component of etzk to real space, destroy plan2
#ifdef DP
    call fftw_execute_dft_c2r(plan2, etzk_x, f)
    mgrid(l)%etz_x(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%etz_x(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%etz_x(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%etz_x(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%etz_x(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%etz_x(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%etz_x(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%etz_x(m+1,m+1,m+1) = f(1,1,1) !last point

    call fftw_execute_dft_c2r(plan2, etzk_y, f)
    mgrid(l)%etz_y(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%etz_y(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%etz_y(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%etz_y(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%etz_y(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%etz_y(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%etz_y(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%etz_y(m+1,m+1,m+1) = f(1,1,1) !last point

    call fftw_execute_dft_c2r(plan2, etzk_z, f)
    mgrid(l)%etz_z(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%etz_z(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%etz_z(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%etz_z(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%etz_z(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%etz_z(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%etz_z(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%etz_z(m+1,m+1,m+1) = f(1,1,1) !last point
    call fftw_destroy_plan(plan2)
#else
    call fftw_execute_dft_c2r(plan2, etzk_x, f)
    mgrid(l)%etz_x(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%etz_x(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%etz_x(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%etz_x(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%etz_x(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%etz_x(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%etz_x(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%etz_x(m+1,m+1,m+1) = f(1,1,1) !last point
   
    call fftw_execute_dft_c2r(plan2, etzk_y, f)
    mgrid(l)%etz_y(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%etz_y(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%etz_y(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%etz_y(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%etz_y(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%etz_y(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%etz_y(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%etz_y(m+1,m+1,m+1) = f(1,1,1) !last point

    call fftw_execute_dft_c2r(plan2, etzk_z, f)
    mgrid(l)%etz_z(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(l)%etz_z(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(l)%etz_z(1:m,m+1,1:m) = f(:,1,:)
    mgrid(l)%etz_z(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(l)%etz_z(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(l)%etz_z(m+1,m+1,1:m) = f(1,1,:)
    mgrid(l)%etz_z(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(l)%etz_z(m+1,m+1,m+1) = f(1,1,1) !last point
    
    call fftwf_destroy_plan(plan2)
#endif

    ! deallocate auxiliary matrices
    deallocate (f)
    deallocate (dbxk)
    deallocate (dbyk)
    deallocate (dbzk)
    deallocate (etzk_x)
    deallocate (etzk_y)
    deallocate (etzk_z)

  enddo


  ! ------------------------------------------------------------------------
  ! calculate displacement field dr from equation 
  ! et = dr x b in each grid level
  ! (need to assume dr . b = 0)
  ! 3d SO cross product
  ! ------------------------------------------------------------------------
  do l = 1, ngrids

    do k = 1, n
      do j = 1, n
        do i = 1, n

          b(1) = mgrid(l)%bx(i,j,k) + mgrid(l)%dbx(i,j,k)
          b(2) = mgrid(l)%by(i,j,k) + mgrid(l)%dby(i,j,k)
          b(3) = mgrid(l)%bz(i,j,k) + mgrid(l)%dbz(i,j,k)

          b2 = b(1)**2 + b(2)**2+b(3)**2 !magnitude of 3d B field

          mgrid(l)%drx(i,j,k) = (b(2)*mgrid(l)%etz_z(i,j,k) - b(3)*mgrid(l)%etz_y(i,j,k))/b2
          mgrid(l)%dry(i,j,k) = (b(3)*mgrid(l)%etz_x(i,j,k) - b(1)*mgrid(l)%etz_z(i,j,k))/b2
          mgrid(l)%drz(i,j,k) = (b(1)*mgrid(l)%etz_y(i,j,k) - b(2)*mgrid(l)%etz_x(i,j,k))/b2

        enddo
      enddo
    enddo
  enddo


  ! ------------------------------------------------------------------------
  ! write drx and dry in each grid to file - UNSURE FOR 3D
  ! ------------------------------------------------------------------------
  
  !NOT SURE ABOUT THIS WRITING IN 3D - COMMENTED OUT FOR NOW -MAINLY THE BIT AFTER 400 IN WRITE IN NESTED LOOP

  ! do l = 1, ngrids
  !
  !   write (file_out, "('DISPLACEMENT_GRID', i0, '.DAT')") l
  !   file_out = trim(data_dir) // trim(file_out)
  !   open(unit=400, file=trim(file_out), form='formatted', status='replace', action='write')
  !   do k = 1, n  
  !     do i = 1, n
  !       do j = 1, n
  !         write(400,'(4(es24.16, 1x))') real(i - 1)*h, real(j - 1)*h,real(j - 1)*h, mgrid(l)%drx(i,j,k), mgrid(l)%dry(i,j,k), mgrid(l)%drz(i,j,k)
  !       enddo
  !       write(400,*)
  !     enddo
  !     write(400,*)
  !   enddo
  !   close(400)
  ! enddo


  ! ------------------------------------------------------------------------
  ! calculate origin points rx0, ry0 and rz0
  ! ------------------------------------------------------------------------
  allocate (rx0(n,n,n))
  allocate (ry0(n,n,n))
  allocate (rz0(n,n,n))

  do k = 1, n
    do j = 1, n
      do i = 1, n

        rxp = real(i - 1)*h
        ryp = real(j - 1)*h
        rzp = real(k - 1)*h
  
        rx0(i,j,k) = rxp
        ry0(i,j,k) = ryp
        ry0(i,j,k) = rzp

        do l = 1, ngrids

          ! gives account of the periodicity
          rxp = mod(rxp, twopi)
          ryp = mod(ryp, twopi)
          rzp = mod(rzp, twopi)

          if (rxp < 0.) rxp = rxp + twopi
          if (ryp < 0.) ryp = ryp + twopi
          if (rzp < 0.) rzp = rzp + twopi

          ! calculate the left lower corner grid point indices
          ii = floor(rxp/h) + 1
          jj = floor(ryp/h) + 1
          k_k = floor(rzp/h) + 1

          ! calculate the right upper corner grid point indices
          iip1 = ii + 1
          jjp1 = jj + 1
          kkp1 = k_k + 1

          if (iip1 > n) iip1 = 2
          if (jjp1 > n) jjp1 = 2
          if (kkp1 > n) kkp1 = 2

          ! calculate linear weigths for the interpolation
          wx1 = mod(rxp, h)/h !x_d
          wy1 = mod(ryp, h)/h !y_d
          wz1 = mod(rzp, h)/h !z_d

          wx0 = 1.d0 - wx1
          wy0 = 1.d0 - wy1
          wz0 = 1.d0 - wz1

          ! perform three sets of trilinear interpolation

          !rx - linear in first index (following steps in trilinear wikipedia page)
          c_00 = wx0*mgrid(l)%drx(ii,jj,k_k) + wx1*mgrid(l)%drx(iip1,jj,k_k)
          c_01 = wx0*mgrid(l)%drx(ii,jj,kkp1) + wx1*mgrid(l)%drx(iip1,jj,kkp1)
          c_10 = wx0*mgrid(l)%drx(ii,jjp1,k_k) + wx1*mgrid(l)%drx(iip1,jjp1,k_k)
          c_11 = wx0*mgrid(l)%drx(ii,jjp1,kkp1) + wx1*mgrid(l)%drx(iip1,jjp1,kkp1)
          !linear in second index
          c_0 = c_00*wy0 + c_10*wy1
          c_1 = c_01*wy0 + c_11*wy1
          !final linear in last indez
          c = c_0*wz0 + c_1*wz1
          drxp = c

          !ry - linear interp in first index
          c_00 = wx0*mgrid(l)%dry(ii,jj,k_k) + wx1*mgrid(l)%dry(iip1,jj,k_k)
          c_01 = wx0*mgrid(l)%dry(ii,jj,kkp1) + wx1*mgrid(l)%dry(iip1,jj,kkp1)
          c_10 = wx0*mgrid(l)%dry(ii,jjp1,k_k) + wx1*mgrid(l)%dry(iip1,jjp1,k_k)
          c_11 = wx0*mgrid(l)%dry(ii,jjp1,kkp1) + wx1*mgrid(l)%dry(iip1,jjp1,kkp1)
          !linear in second index
          c_0 = c_00*wy0 + c_10*wy1
          c_1 = c_01*wy0 + c_11*wy1
          !final linear in last index
          c = c_0*wz0 + c_1*wz1
          dryp = c

          !rz - linear interp in first index
          c_00 = wx0*mgrid(l)%drz(ii,jj,k_k) + wx1*mgrid(l)%drz(iip1,jj,k_k)
          c_01 = wx0*mgrid(l)%drz(ii,jj,kkp1) + wx1*mgrid(l)%drz(iip1,jj,kkp1)
          c_10 = wx0*mgrid(l)%drz(ii,jjp1,k_k) + wx1*mgrid(l)%drz(iip1,jjp1,k_k)
          c_11 = wx0*mgrid(l)%drz(ii,jjp1,kkp1) + wx1*mgrid(l)%drz(iip1,jjp1,kkp1)
          !linear in second index
          c_0 = c_00*wy0 + c_10*wy1
          c_1 = c_01*wy0 + c_11*wy1
          !final linear in last index
          c = c_0*wz0 + c_1*wz1
          drzp = c

          ! new normalised positions (periodic)
          rxp = rxp + drxp
          ryp = ryp + dryp
          rzp = rzp + drzp

          !new positions
          rx0(i,j,k) = rx0(i,j,k) + drxp
          ry0(i,j,k) = ry0(i,j,k) + dryp
          rz0(i,j,k) = rz0(i,j,k) + drzp

        enddo
        
      enddo
    enddo
  enddo

  ! ------------------------------------------------------------------------
  ! calculate total displacement fields drx and dry and write to file - UNSURE FOR 3D
  ! ------------------------------------------------------------------------
  allocate (drx(n,n,n))
  allocate (dry(n,n,n))
  allocate (drz(n,n,n))

  do k = 1, n
    do j = 1, n
      do i = 1, n
        drx(i,j,k) = rx0(i,j,k) - real(i - 1)*h
        dry(i,j,k) = ry0(i,j,k) - real(j - 1)*h
        drz(i,j,k) = ry0(i,j,k) - real(k - 1)*h
      enddo
    enddo
  enddo
  ! file_out = trim(data_dir) // 'DISPLACEMENT.DAT'
  ! open(unit=400, file=trim(file_out), form='formatted', status='replace', action='write')
  ! do i = 1, n
  !   do j = 1, n
  !     write(400,'(4(es24.16, 1x))') real(i - 1)*h, real(j - 1)*h, drx(i,j), dry(i,j)
  !   enddo
  !   write(400,*)
  ! enddo
  ! close(400)


  ! ------------------------------------------------------------------------
  ! build scalar field phi0 and write to file
  ! ------------------------------------------------------------------------
  m = n !- 1

  allocate (phi0(n,n,n))
  allocate (x_arr(n,n,n)) 
  allocate (y_arr(n,n,n)) 
  allocate (z_arr(n,n,n))
  allocate (input_1(n,n,n))   
  ! allocate (phi0k((m/2 + 1), m))

  do i = 1, n
    x_arr(i,:,:) = i*twopi/n
  enddo
  do j = 1, n
    y_arr(:,j,:) = j*twopi/n
  enddo
  do k = 1, n
    z_arr(:,:,k) = k*twopi/n
  enddo

  phi0(:,:,:) = 0 ! initialise all entries to zero

  print*, omp_get_max_threads()
  wtime = omp_get_wtime()
  
  !not thread safe - phi0 magnitudes greater when using OpenMP - distributed memory also not good for extending into much larger scales

  call omp_set_num_threads(30)
  !$OMP PARALLEL
  do ki = 0, n-3 
    kx = (-(n-1)/2 + 1) + ki
    print*, kx
    if (ki == 2) then
      threadno = omp_get_num_threads() !check to see if openmp working
      print*, "Total running threads", threadno
    endif
    thread_id = omp_get_thread_num()
    if (thread_id == 2) then
      print*, "ki value", ki
    endif
    
    do kj = 0, n-3 ! up to nyquist frequency
      ky = (-(n-1)/2 + 1) + kj

      do kk = 0, n-3 !3d
        kz = (-(n-1)/2 + 1) + kk

        if ((abs (ky) < 1.0D-5 .and. abs(kz) < 1.0D-5)) then !cant root 0 - now for 3d
            continue

        else

          call random_number(num)
          tmp = abs(ky**2+kz**2)**(-10.0d0/6.0d0) !3D
          tmp2 = exp(-(twopi)**(1.0d0/3.0d0)*abs(kx)/(abs(ky**2+kz**2)**(2.0d0/6.0d0)))
          amp = sqrt(tmp*tmp2) !amplitude

          input_1 = kx*x_arr + ky*y_arr + kz*z_arr + num*twopi

          phi0(:,:,:) = phi0(:,:,:) + amp*cos(input_1) !3d
          
        endif
      
      enddo  
    enddo
  enddo
  !$OMP END PARALLEL
  print*, 'The loop has successfully completed'

  wtime = omp_get_wtime() - wtime

  print *, wtime 

  print*, '* Writing file: phi_0'

  ! file_out = trim(data_dir) // 'PHI0.DAT'
  ! open(unit=400, file=trim(file_out), form='formatted', status='replace', action='write')
  ! do i = 1, n
  !   do j = 1, n
  !     write(400,'(3(es24.16, 1x))') real(i - 1)*h, real(j - 1)*h, phi0(i,j)
  !   enddo
  !   write(400,*)
  ! enddo
  ! close(400)

  file_out = 'PHI0.BIN'
  lun = 701
  file_out = trim(data_dir) // '/' // file_out
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) phi0(:,:,:)
  close(lun)


  !----------------------------------------------------------------
  ! calculate power spectrum of phi0 and write to file
  !----------------------------------------------------------------
  lsum_power = .false.

  tmp3d(:,:,:) = phi0(:,:,:)
  !call power_kk(tmp3d, ps_k, ps_kb, ps_kt, ps_kk, lsum_power, lx, ly, lz)
  !call write_power_spectra('PHI0.DAT', 400) - no DAT file created yet


  ! ------------------------------------------------------------------------
  ! remap scalar field phi0 into phi and write to file
  ! ------------------------------------------------------------------------
  allocate (phi(n,n,n))

  do k = 1, n
    do j = 1, n
      do i = 1, n

        rxp = real(i - 1)*h
        ryp = real(j - 1)*h
        rzp = real(k - 1)*h
  
        rx0(i,j,k) = rxp
        ry0(i,j,k) = ryp
        ry0(i,j,k) = rzp

        ! gives account of the periodicity
        rxp = mod(rxp, twopi)
        ryp = mod(ryp, twopi)
        rzp = mod(rzp, twopi)

        if (rxp < 0.) rxp = rxp + twopi
        if (ryp < 0.) ryp = ryp + twopi
        if (rzp < 0.) rzp = rzp + twopi

        ! calculate the left lower corner grid point indices
        ii = floor(rxp/h) + 1
        jj = floor(ryp/h) + 1
        k_k = floor(rzp/h) + 1

        ! calculate the right upper corner grid point indices
        iip1 = ii + 1
        jjp1 = jj + 1
        kkp1 = k_k + 1

        if (iip1 > n) iip1 = 2
        if (jjp1 > n) jjp1 = 2
        if (kkp1 > n) kkp1 = 2

        ! calculate linear weigths for the interpolation
        wx1 = mod(rxp, h)/h !x_d
        wy1 = mod(ryp, h)/h !y_d
        wz1 = mod(rzp, h)/h !z_d

        wx0 = 1.d0 - wx1
        wy0 = 1.d0 - wy1
        wz0 = 1.d0 - wz1

        !trilinear interpolation for phi

        !linear interp in first index
        c_00 = wx0*phi0(ii,jj,k_k) + wx1*phi0(iip1,jj,k_k)
        c_01 = wx0*phi0(ii,jj,kkp1) + wx1*phi0(iip1,jj,kkp1)
        c_10 = wx0*phi0(ii,jjp1,k_k) + wx1*phi0(iip1,jjp1,k_k)
        c_11 = wx0*phi0(ii,jjp1,kkp1) + wx1*phi0(iip1,jjp1,kkp1)
        !linear in second index
        c_0 = c_00*wy0 + c_10*wy1
        c_1 = c_01*wy0 + c_11*wy1
        !final linear in last index
        c = c_0*wz0 + c_1*wz1
        
        phi(i,j,k) = c
        
      enddo
    enddo
  enddo

  !print*, phi(23,:)

  ! file_out = trim(data_dir) // 'PHI.DAT'
  ! open(unit=400, file=trim(file_out), form='formatted', status='replace', action='write')
  ! do i = 1, n
  !   do j = 1, n
  !     write(400,'(3(es24.16, 1x))') real(i - 1)*h, real(j - 1)*h, phi(i,j)
  !   enddo
  !   write(400,*)
  ! enddo
  ! close(400)

  file_out = 'PHI.BIN'
  lun = 701
  file_out = trim(data_dir) // '/' // file_out
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) phi(:,:,:)
  close(lun)

  !----------------------------------------------------------------
  ! calculate power spectrum of phi and write to file
  !----------------------------------------------------------------
  lsum_power = .false.

  tmp3d(:,:,:) = phi(:,:,:)
  !call power_kk(tmp3d, ps_k, ps_kb, ps_kt, ps_kk, lsum_power, lx, ly, lz)
  !call write_power_spectra('PHI.DAT', 400) - no DAT files yet


  ! ------------------------------------------------------------------------
  ! data for Reinaldo's structure function code 
  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------
  ! ATTENTION: SAVE FIELD IN DN (TEST-PURPOSE ONLY!)
  ! ------------------------------------------------------------------------

  !3d
  lun = 700
  file_out = trim(data_dir) // 'DN0.BIN'
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write')
  write(lun) n, n, n

  tmp3d(:,:,:) = phi(:,:,:)

  write(lun) tmp3d(:,:,:)
  write(lun) time, x, y, z, dx, dy, dz
  close(lun)


  lun = 701
  file_out = trim(data_dir) // 'BB0.BIN'
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write')
  write(lun) n, n, n

  tmp3d(:,:,:) = bx(:,:,:)
  write(lun) tmp3d(:,:,:)

  tmp3d(:,:,:) = by(:,:,:)
  write(lun) tmp3d(:,:,:)
  
  tmp3d(:,:,:) = bz(:,:,:)
  write(lun) tmp3d(:,:,:)

  write(lun) time, x, y, z, dx, dy, dz
  close(lun)


  lun = 702
  file_out = trim(data_dir) // 'VV0.BIN'
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write')
  write(lun) n, n, n
  write(lun) tmp3d(:,:,:)
  write(lun) tmp3d(:,:,:)
  write(lun) tmp3d(:,:,:)
  write(lun) time, x, y, z, dx, dy, dz
  close(lun)


  ! ------------------------------------------------------------------------
  ! deallocate memory
  ! ------------------------------------------------------------------------
  do l = 1, ngrids
    deallocate (mgrid(l)%bx)
    deallocate (mgrid(l)%by)
    deallocate (mgrid(l)%bz)
    deallocate (mgrid(l)%dbx)
    deallocate (mgrid(l)%dby)
    deallocate (mgrid(l)%dbz)
    deallocate (mgrid(l)%etz_x)
    deallocate (mgrid(l)%etz_y)
    deallocate (mgrid(l)%etz_z)
    deallocate (mgrid(l)%drx)
    deallocate (mgrid(l)%dry)
    deallocate (mgrid(l)%drz)
  enddo
  deallocate (mgrid)

  deallocate (bx)
  deallocate (by)
  deallocate (bz)
  deallocate (rx0)
  deallocate (ry0)
  deallocate (rz0)
  deallocate (drx)
  deallocate (dry)
  deallocate (drz)

  deallocate(x)
  deallocate(y)
  deallocate(z)

  deallocate (tmp3d)

  deallocate (ps_k)
  deallocate (ps_kb)
  deallocate (ps_kt)
  deallocate (ps_kk)

  deallocate (phi0)
  deallocate (phi)
  deallocate (x_arr)
  deallocate (y_arr)
  deallocate (z_arr)
  deallocate (input_1)

  stop

  ! contains

  !   !----------------------------------------------------------------
  !   ! calculate power spectrum of phi and write to file
  !   !----------------------------------------------------------------
  !   subroutine write_power_spectra(partial_file_name, lun)

  !     implicit none

  !     character(len=*) :: partial_file_name
  !     integer :: lun

  !     ! local dummy variables
  !     integer :: k, kb, kt
  !     character(len=400) :: file_name


  !     file_name = trim(data_dir) // 'power_spectra/PS_K_' // trim(partial_file_name)
  !     open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
  !     write(lun, "(a1, a4, 1x, 1(a24, 1x))") '#', 'k', &
  !       '|phi(k)|^2'
  !     do k = 1, nk
  !       write(lun, '(i5, 1x, 1(es24.16, 1x))') k-1, &
  !         ps_k(k)
  !     enddo
  !     close(lun)

  !     file_name = trim(data_dir) // 'power_spectra/PS_KB_' // trim(partial_file_name)
  !     open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
  !     write(lun, "(a1, a4, 1x, 1(a24, 1x))") '#', 'kx', &
  !       '|phi(kx)|^2'
  !     do k = 1, nkb
  !       write(lun, '(i5, 1x, 1(es24.16, 1x))') k-1, &
  !         ps_kb(k)
  !     enddo
  !     close(lun)

  !     file_name = trim(data_dir) // 'power_spectra/PS_KT_' // trim(partial_file_name)
  !     open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
  !     write(lun, "(a1, a4, 1x, 1(a24, 1x))") '#', 'kyz', &
  !       '|phi(kyz)|^2'
  !     do k = 1, nkt
  !       write(lun, '(i5, 1x, 1(es24.16, 1x))') k-1, &
  !         ps_kt(k)
  !     enddo
  !     close(lun)

  !     file_name = trim(data_dir) // 'power_spectra/PS_KK_' // trim(partial_file_name)
  !     open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
  !     write(lun, "(a1, a4, 1x, a5, 1x, 1(a24, 1x))") '#', 'kx', 'kyz', &
  !       '|phi(kx, kyz)|^2'
  !     do kt = 1, nkt
  !       do kb = 1, nkb
  !         write(lun, '(2(i5, 1x), 1(es24.16, 1x))') kb-1, kt-1, &
  !           ps_kk(kb,kt)
  !       enddo
  !       write(lun,*)
  !     enddo
  !     close(lun)

  !     return
  !   end subroutine write_power_spectra

end program main
