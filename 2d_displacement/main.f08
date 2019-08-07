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
    real(sp), dimension(:,:), allocatable :: bx
    real(sp), dimension(:,:), allocatable :: by
    real(sp), dimension(:,:), allocatable :: dbx
    real(sp), dimension(:,:), allocatable :: dby
    real(sp), dimension(:,:), allocatable :: etz
    real(sp), dimension(:,:), allocatable :: drx
    real(sp), dimension(:,:), allocatable :: dry
  end type sgrid


  ! ------------------------------------------------------------------------
  ! define and initialize problem parameters
  ! ------------------------------------------------------------------------
  integer :: ngrids = 9
  real(sp) :: bx0 = 1.
  real(sp) :: by0 = 0.
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
  real(sp) :: kx
  real(sp) :: ky
  real(sp), dimension(:,:), allocatable :: bx, by
  real(sp), dimension(:,:), allocatable :: rx0, ry0
  real(sp), dimension(:,:), allocatable :: drx, dry
  type(sgrid), dimension(:), allocatable :: mgrid

  real(sp), dimension(:), allocatable :: x, y, z
  real(sp) :: time, dx, dy, dz

  real(sp), dimension(:,:), allocatable :: phi0
  real(sp), dimension(:,:), allocatable :: phi
  complex(sp), dimension(:,:), allocatable :: phi0k

  integer :: nk, nkb, nkt
  real(dp), dimension(:), allocatable :: ps_k
  real(dp), dimension(:), allocatable :: ps_kb
  real(dp), dimension(:), allocatable :: ps_kt
  real(dp), dimension(:,:), allocatable :: ps_kk
  real(dp), dimension(:,:,:), allocatable :: tmp3d

  logical :: lsum_power

  ! ------------------------------------------------------------------------
  ! define auxiliary variables
  ! ------------------------------------------------------------------------
  type(C_PTR) :: plan_phi0
  type(C_PTR) :: plan1, plan2
  real(sp), dimension(:,:), allocatable :: f
  complex(sp), dimension(:,:), allocatable :: bxk, byk
  complex(sp), dimension(:,:), allocatable :: dbxk, dbyk
  complex(sp), dimension(:,:), allocatable :: etzk

  real(sp) :: b(2), b2
  real(sp) :: b1(2), ph1, ph2, aux1, aux2
  real(sp) :: rxp, ryp, drxp, dryp
  real(sp) :: wx0, wx1, wy0, wy1
  real(sp) :: kmax, kmod

  real(sp) :: k_para, k_perp, E_coeff, ph


  ! ------------------------------------------------------------------------
  ! define dummy variables
  ! ------------------------------------------------------------------------
  integer :: l
  integer :: m
  integer :: i, j, k
  integer :: ii, jj
  integer :: iip1, jjp1
  integer :: ki, kj
  integer :: lun
  real(sp) :: tmp, tmp2

  character(len=400) :: data_dir
  character(len=1024) :: file_out
  character(len=400) :: cmd


  ! ------------------------------------------------------------------------
  ! specify folder for output data
  ! ------------------------------------------------------------------------
  data_dir = './512run2D_73_paratest/'
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
  allocate (tmp3d(n,n,1))

  allocate (ps_k(nk))
  allocate (ps_kb(nkb))
  allocate (ps_kt(nkt))
  allocate (ps_kk(nkb,nkt))

  allocate(x(n))
  allocate(y(n))
  allocate(z(1))

  ! ------------------------------------------------------------------------
  ! set space grid and time of output files
  ! ------------------------------------------------------------------------
  dx = lx/real(n)
  dy = ly/real(n)
  dz = lz/real(1)

  do i = 1, n
    x(i) = (real(i) - 0.5)*dx
  enddo

  do j = 1, n
    y(j) = (real(j) - 0.5)*dy
  enddo

  do k = 1, 1
    z(k) = (real(k) - 0.5)*dz
  enddo

  time = 0.


  ! ------------------------------------------------------------------------
  ! build bx and by and write to file
  ! ------------------------------------------------------------------------
  allocate (bx(n,n))
  allocate (by(n,n))

  bx(:,:) = bx0
  by(:,:) = by0

  do j = 1, n
    do i = 1, n
      by(i,j) = by(i,j) + 0.5*sin(2.0*x(i))
      by(i,j) = by(i,j) + 0.5*sin(4.0*x(i)+1.6)
    enddo
  enddo

  ! ! build a random magnetic field
  ! do kj = 0, 2
  !   do ki = 0, 2

  !     kmod = sqrt(real(ki)**2 + real(kj)**2)

  !     ! only continue if kmod is inside interval [1,2]
  !     if ((kmod < 1) .or. (kmod > 2)) cycle

  !     ! sort a random direction perpendicular to k
  !     call random_number(b(:))
  !     b(1) = b(1) - (b(1)*real(ki) + b(2)*real(kj))*real(ki)/(real(ki)**2 + real(kj)**2)
  !     b(2) = b(2) - (b(1)*real(ki) + b(2)*real(kj))*real(kj)/(real(ki)**2 + real(kj)**2)
  !     b(:) = b(:) / sqrt( sum(b(:)**2) )

  !     call random_number(b1(:))
  !     b1(1) = b1(1) - (b1(1)*real(ki) + b1(2)*real(kj))*real(ki)/(real(ki)**2 + real(kj)**2)
  !     b1(2) = b1(2) - (b1(1)*real(ki) + b1(2)*real(kj))*real(kj)/(real(ki)**2 + real(kj)**2)
  !     b1(:) = b1(:) / sqrt( sum(b1(:)**2) )

  !     ! sort random phases
  !     call random_number(ph1)
  !     ph1 = ph1*twopi

  !     call random_number(ph2)
  !     ph2 = ph2*twopi

  !     do j = 1, n
  !       do i = 1, n

  !         aux1 = b(1)*cos(real(ki)*x(i) + real(kj)*y(j) + ph1)
  !         aux2 = b(2)*cos(real(ki)*x(i) + real(kj)*y(j) + ph1)

  !         if ((ki > 0) .and. (kj > 0)) then
  !           aux1 = aux1 + b1(1)*cos(real(ki)*x(i) - real(kj)*y(j) + ph2)
  !           aux2 = aux2 + b1(2)*cos(real(ki)*x(i) - real(kj)*y(j) + ph2)
  !         endif

  !         if (kmod > 0.) then
  !           aux1 = 2.*aux1
  !           aux2 = 2.*aux2
  !         endif

  !         bx(i,j) = bx(i,j) + aux1
  !         by(i,j) = by(i,j) + aux2

  !       enddo
  !     enddo

  !   enddo ! ki
  ! enddo ! kj

  !print*, by(:,1)  ! printing by in x direction - should see variation

  file_out = trim(data_dir) // 'MAGNETIC.DAT' 
  open(unit=400, file=trim(file_out), form='formatted', status='replace', action='write')
  do i = 1, n
    do j = 1, n
      write(400,'(4(es24.16, 1x))') real(i - 1)*h, real(j - 1)*h, bx(i,j), by(i,j)
    enddo
    write(400,*)
  enddo
  close(400)

  print*, '* Writing file: magnetic field'
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

  lun = 701
  file_out = trim(data_dir) // '/' // 'BX.BIN'
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) bx(:,:)
  close(lun)

  lun = 701
  file_out = trim(data_dir) // '/' // 'BY.BIN'
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) by(:,:)
  close(lun)


  ! ------------------------------------------------------------------------
  ! generate grids hierarchy and allocate memory for each grid level
  ! ------------------------------------------------------------------------
  allocate (mgrid(ngrids))

  do l = 1, ngrids
    allocate (mgrid(l)%bx(n,n))
    allocate (mgrid(l)%by(n,n))
    allocate (mgrid(l)%dbx(n,n))
    allocate (mgrid(l)%dby(n,n))
    allocate (mgrid(l)%etz(n,n))
    allocate (mgrid(l)%drx(n,n))
    allocate (mgrid(l)%dry(n,n))
  enddo


  ! ------------------------------------------------------------------------
  ! build bx and by in each grid level, using spectral filtering
  ! ------------------------------------------------------------------------
  mgrid(1)%bx(:,:) = bx(:,:)
  mgrid(1)%by(:,:) = by(:,:)

  do l = 2, ngrids

    kmax = sqrt(2.d0)*(2.d0**(ngrids - l))

    m = n - 1

    ! allocate auxiliary matrices
    allocate (f(m, m))
    allocate (bxk((m/2 + 1), m))
    allocate (byk((m/2 + 1), m))

    ! prepare plans for the dft (plan1) and dft inverse (plan2)
#ifdef DP
    plan1 = fftw_plan_dft_r2c_2d(m, m, f, bxk, FFTW_ESTIMATE)
    plan2 = fftw_plan_dft_c2r_2d(m, m, bxk, f, FFTW_ESTIMATE)
#else
    plan1 = fftwf_plan_dft_r2c_2d(m, m, f, bxk, FFTW_ESTIMATE)
    plan2 = fftwf_plan_dft_c2r_2d(m, m, bxk, f, FFTW_ESTIMATE)
#endif

    ! transform mgrid(l-1)%bx and mgrid(l-1)%by to fourier space, destroy plan1
#ifdef DP
    f(:,:) = mgrid(l-1)%bx(1:m,1:m)
    call fftw_execute_dft_r2c(plan1, f, bxk)
    !r2c = real to complex
    f(:,:) = mgrid(l-1)%by(1:m,1:m)
    call fftw_execute_dft_r2c(plan1, f, byk)

    call fftw_destroy_plan(plan1)
#else
    f(:,:) = mgrid(l-1)%bx(1:m,1:m)
    call fftwf_execute_dft_r2c(plan1, f, bxk)

    f(:,:) = mgrid(l-1)%by(1:m,1:m)
    call fftwf_execute_dft_r2c(plan1, f, byk)

    call fftwf_destroy_plan(plan1)
#endif

    ! normalize bxk and byk
    bxk(:,:) = bxk(:,:) / real(m*m)
    byk(:,:) = byk(:,:) / real(m*m)

    ! loop over all the modes
    ! (i, j) are the indices of the matrices generated by the dft
    ! corresponding to the wave numbers (ki, kj)
    do kj = min((-m/2 + 1), 0), m/2
      if (kj >= 0) then
        j = kj + 1
      else
        j = m + kj + 1
      endif

      do ki = 0, m/2
        i = ki + 1

        kmod = sqrt(real(ki)**2 + real(kj)**2)

        if (kmod > kmax) then
          bxk(i,j) = (0.d0, 0.d0)
          byk(i,j) = (0.d0, 0.d0)
        endif

      enddo
    enddo

    ! transform bxk and byk to real space, apply periodic boundaries, destroy plan2
#ifdef DP
    call fftw_execute_dft_c2r(plan2, bxk, f)
    !c2r = complex to real
    mgrid(l)%bx(1:m,1:m) = f(:,:)
    mgrid(l)%bx(m+1,1:m) = f(1,:)
    mgrid(l)%bx(1:m,m+1) = f(:,1)
    mgrid(l)%bx(m+1,m+1) = f(1,1)

    call fftw_execute_dft_c2r(plan2, byk, f)
    mgrid(l)%by(1:m,1:m) = f(:,:)
    mgrid(l)%by(m+1,1:m) = f(1,:)
    mgrid(l)%by(1:m,m+1) = f(:,1)
    mgrid(l)%by(m+1,m+1) = f(1,1)

    call fftw_destroy_plan(plan2)
#else
    call fftwf_execute_dft_c2r(plan2, bxk, f)
    mgrid(l)%bx(1:m,1:m) = f(:,:)
    mgrid(l)%bx(m+1,1:m) = f(1,:)
    mgrid(l)%bx(1:m,m+1) = f(:,1)
    mgrid(l)%bx(m+1,m+1) = f(1,1)

    call fftwf_execute_dft_c2r(plan2, byk, f)
    mgrid(l)%by(1:m,1:m) = f(:,:)
    mgrid(l)%by(m+1,1:m) = f(1,:)
    mgrid(l)%by(1:m,m+1) = f(:,1)
    mgrid(l)%by(m+1,m+1) = f(1,1)

    call fftwf_destroy_plan(plan2)
#endif

    ! deallocate auxiliary matrices
    deallocate (f)
    deallocate (bxk)
    deallocate (byk)

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
      write(lun) mgrid(l)%bx(:,:)
    close(lun)

    write (file_out, "('BY_GRID', i0, '.BIN')") l
    !why 702 all of a sudden??
    lun = 702
    file_out = trim(data_dir) // '/' // file_out
    open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
      write(lun) mgrid(l)%by(:,:)
    close(lun)

  enddo


  ! ------------------------------------------------------------------------
  ! calculate vector field db in each grid level
  ! ------------------------------------------------------------------------
  do l = 1, ngrids

    do j = 1, n
      do i = 1, n

        if (l < ngrids) then
          b(1) = mgrid(l+1)%bx(i,j)
          b(2) = mgrid(l+1)%by(i,j)
        else
          ! use average uniform field
          b(1) = bx0
          b(2) = by0
        endif

        mgrid(l)%dbx(i,j) = b(1) - mgrid(l)%bx(i,j)
        mgrid(l)%dby(i,j) = b(2) - mgrid(l)%by(i,j)

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
    allocate (f(m, m))
    allocate (dbxk((m/2 + 1), m))
    allocate (dbyk((m/2 + 1), m))
    allocate (etzk((m/2 + 1), m))

    ! prepare plans for the dft (plan1) and dft inverse (plan2)
#ifdef DP
    plan1 = fftw_plan_dft_r2c_2d(m, m, f, dbxk, FFTW_ESTIMATE)
    plan2 = fftw_plan_dft_c2r_2d(m, m, etzk, f, FFTW_ESTIMATE)
#else
    plan1 = fftwf_plan_dft_r2c_2d(m, m, f, dbxk, FFTW_ESTIMATE)
    plan2 = fftwf_plan_dft_c2r_2d(m, m, etzk, f, FFTW_ESTIMATE)
#endif

    ! transform dbx and dby to fourier space, destroy plan1
#ifdef DP
    f(:,:) = mgrid(l)%dbx(1:m,1:m)
    call fftw_execute_dft_r2c(plan1, f, dbxk)

    f(:,:) = mgrid(l)%dby(1:m,1:m)
    call fftw_execute_dft_r2c(plan1, f, dbyk)

    call fftw_destroy_plan(plan1)
#else
    f(:,:) = mgrid(l)%dbx(1:m,1:m)
    call fftwf_execute_dft_r2c(plan1, f, dbxk)

    f(:,:) = mgrid(l)%dby(1:m,1:m)
    call fftwf_execute_dft_r2c(plan1, f, dbyk)

    call fftwf_destroy_plan(plan1)
#endif

    ! normalize dbxk and dbyk
    dbxk(:,:) = dbxk(:,:) / real(m*m)
    dbyk(:,:) = dbyk(:,:) / real(m*m)

    ! loop over all the modes
    ! (i, j) are the indices of the matrices generated by the dft
    ! corresponding to the wave numbers (ki, kj)
    do kj = min((-m/2 + 1), 0), m/2
      if (kj >= 0) then
        j = kj + 1
      else
        j = m + kj + 1
      endif

      do ki = 0, m/2
        i = ki + 1

        if ((ki==0) .and. (kj==0)) then
          etzk(i,j) = (0., 0.)
          cycle
        endif

        ! etk = - i (dbk x k) / k^2
        etzk(i,j) = - (0., 1.)*(dbxk(i,j)*real(kj) - dbyk(i,j)*real(ki))/(real(ki)**2 + real(kj)**2)

      enddo
    enddo

    ! transform etzk to real space, destroy plan2
#ifdef DP
    call fftw_execute_dft_c2r(plan2, etzk, f)
    call fftw_destroy_plan(plan2)
#else
    call fftwf_execute_dft_c2r(plan2, etzk, f)
    call fftwf_destroy_plan(plan2)
#endif

    mgrid(l)%etz(1:m,1:m) = f(:,:)

    ! apply periodic boundaries
    mgrid(l)%etz(m+1,1:m) = f(1,:)
    mgrid(l)%etz(1:m,m+1) = f(:,1)
    mgrid(l)%etz(m+1,m+1) = f(1,1)

    ! deallocate auxiliary matrices
    deallocate (f)
    deallocate (dbxk)
    deallocate (dbyk)
    deallocate (etzk)

  enddo


  ! ------------------------------------------------------------------------
  ! calculate displacement field dr from equation 
  ! et = dr x b in each grid level
  ! (need to assume dr . b = 0)
  ! ------------------------------------------------------------------------
  do l = 1, ngrids

    do j = 1, n
      do i = 1, n

        b(1) = mgrid(l)%bx(i,j) + mgrid(l)%dbx(i,j)
        b(2) = mgrid(l)%by(i,j) + mgrid(l)%dby(i,j)
        b2 = b(1)**2 + b(2)**2

        mgrid(l)%drx(i,j) =   b(2)*mgrid(l)%etz(i,j)/b2
        mgrid(l)%dry(i,j) = - b(1)*mgrid(l)%etz(i,j)/b2

      enddo
    enddo

  enddo


  ! ------------------------------------------------------------------------
  ! write drx and dry in each grid to file
  ! ------------------------------------------------------------------------
  do l = 1, ngrids

    write (file_out, "('DISPLACEMENT_GRID', i0, '.DAT')") l
    file_out = trim(data_dir) // trim(file_out)
    open(unit=400, file=trim(file_out), form='formatted', status='replace', action='write')
    do i = 1, n
      do j = 1, n
        write(400,'(4(es24.16, 1x))') real(i - 1)*h, real(j - 1)*h, mgrid(l)%drx(i,j), mgrid(l)%dry(i,j)
      enddo
      write(400,*)
    enddo
    close(400)

  enddo


  ! ------------------------------------------------------------------------
  ! calculate origin points rx0 and ry0
  ! ------------------------------------------------------------------------
  allocate (rx0(n,n))
  allocate (ry0(n,n))

  do j = 1, n
    do i = 1, n

      ! start point for regression
      rxp = real(i - 1)*h
      ryp = real(j - 1)*h

      rx0(i,j) = rxp
      ry0(i,j) = ryp

      ! add displacement contribution from each grid level
      do l = 1, ngrids

        ! interpolate drx and drz at position (rxp, ryp)

        ! gives account of the periodicity
        rxp = mod(rxp, twopi)
        ryp = mod(ryp, twopi)
        if (rxp < 0.) rxp = rxp + twopi
        if (ryp < 0.) ryp = ryp + twopi

        ! calculate the left lower corner grid point indices
        ii = floor(rxp/h) + 1
        jj = floor(ryp/h) + 1

        ! calculate the right upper corner grid point indices
        iip1 = ii + 1
        jjp1 = jj + 1

        if (iip1 > n) iip1 = 2
        if (jjp1 > n) jjp1 = 2

        ! calculate linear weigths for the interpolation
        wx1 = mod(rxp, h)/h
        wy1 = mod(ryp, h)/h

        wx0 = 1.d0 - wx1
        wy0 = 1.d0 - wy1

        ! perform bilinear interpolation
        drxp =   wx0*wy0*mgrid(l)%drx(ii  , jj  ) &
               + wx1*wy0*mgrid(l)%drx(iip1, jj  ) &
               + wx0*wy1*mgrid(l)%drx(ii  , jjp1) &
               + wx1*wy1*mgrid(l)%drx(iip1, jjp1)

        dryp =   wx0*wy0*mgrid(l)%dry(ii  , jj  ) &
               + wx1*wy0*mgrid(l)%dry(iip1, jj  ) &
               + wx0*wy1*mgrid(l)%dry(ii  , jjp1) &
               + wx1*wy1*mgrid(l)%dry(iip1, jjp1)

        rxp = rxp + drxp
        ryp = ryp + dryp

        rx0(i,j) = rx0(i,j) + drxp
        ry0(i,j) = ry0(i,j) + dryp

      enddo

    enddo
  enddo


  ! ------------------------------------------------------------------------
  ! calculate total displacement fields drx and dry and write to file
  ! ------------------------------------------------------------------------
  allocate (drx(n,n))
  allocate (dry(n,n))

  do j = 1, n
    do i = 1, n
      drx(i,j) = rx0(i,j) - real(i - 1)*h
      dry(i,j) = ry0(i,j) - real(j - 1)*h
    enddo
  enddo

  file_out = trim(data_dir) // 'DISPLACEMENT.DAT'
  open(unit=400, file=trim(file_out), form='formatted', status='replace', action='write')
  do i = 1, n
    do j = 1, n
      write(400,'(4(es24.16, 1x))') real(i - 1)*h, real(j - 1)*h, drx(i,j), dry(i,j)
    enddo
    write(400,*)
  enddo
  close(400)


  ! ------------------------------------------------------------------------
  ! build scalar field phi0 and write to file
  ! ------------------------------------------------------------------------
  m = n !- 1

  allocate (phi0(n,n)) 
  ! allocate (phi0k((m/2 + 1), m))

  phi0(:,:) = 0 ! initialise all entries to zero
  
 
  !rand_seed =(/75421/)
  !call random_seed(put=rand_seed)


  ! generate GS95 Spectrum for Strong Alvenic Turbulence (Goldreich-Sridhar 1995)
  ! print *, n

  !---ALTERNATIVE TIMING---

  !call system_clock(start_time, count_rate, count_max)
  !time_init = start_time*1.0/count_rate

  print*, omp_get_max_threads()
  !SHARED(tmp, tmp2, amp, phi0, i, j, kj)
  wtime = omp_get_wtime()

  call omp_set_num_threads(40)

  !$OMP PARALLEL
  !$OMP DO 
  do ki = 0, n-3 
    kx = (-(n-1)/2 + 1) + ki
    if (ki == 2) then
      threadno = omp_get_num_threads() !check to see if openmp working
      print*, threadno
    endif
    thread_id = omp_get_thread_num()
    if (thread_id == 3) then
      print*, ki
    endif
    do kj = 0, n-3 ! up to nyquist frequency
      ky = (-(n-1)/2 + 1) + kj
      call random_number(num)
      if ((kx == 0) .OR. (ky == 0)) then !cant root 0
            continue
      else
        tmp = abs(ky)**(-7/3) !2D
        tmp2 = exp(-abs(kx)/(abs(ky)**(2/3)))
        amp = sqrt(tmp*tmp2) !amplitude
        do i = 1, n
                do j = 1, n
                phi0(i,j) = phi0(i,j) + amp*cos(kx*i*twopi/n + ky*j*twopi/n + num*twopi)
                enddo
        enddo
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  print*, 'The loop has successfully completed'

  wtime = omp_get_wtime() - wtime

  print *, wtime 

  !---ALTERNATIVE TIMING---

  !call system_clock(stop_time, count_rate, count_max)
  !time_final = stop_time*1.0/count_rate
  !elapsed_time = time_final - time_init
  !print *, elapsed_time


  !print*, phi0(23,67), phi0(13,45), phi0(103,31)


  ! generate chess pattern (eight strips)
  ! do j = 1, n
  !   jj = (j - 1)/(n/8)

  !   do i = 1, n
  !     ii = (i - 1)/(n/8)

  !     if (mod(ii,2) == 0) then
  !       if (mod(jj,2) == 0) then
  !         phi0(i,j) = 1.
  !       else
  !         phi0(i,j) = 2.
  !       endif
  !     else
  !       if (mod(jj,2) == 0) then
  !         phi0(i,j) = 2.
  !       else
  !         phi0(i,j) = 1.
  !       endif
  !     endif

  !   enddo
  ! enddo


!     ! create plan
! #ifdef DP
!     plan_phi0 = fftw_plan_dft_c2r_2d (m, m, phi0k, phi0, FFTW_ESTIMATE)
! #else
!     plan_phi0 = fftwf_plan_dft_c2r_2d (m, m, phi0k, phi0, FFTW_ESTIMATE)
! #endif

!     ! ------------------------------------------------------------------------
!     ! loop over wavenumbers to build scalar field phi0k
!     ! ------------------------------------------------------------------------
!     phi0k(:,:) = 0.

!     kmax = real(m)/2.

!     ! SKIP NYQUIST FREQUENCY
!     do kj = min((-m/2 + 1), 0), m/2-1
!       if (kj >= 0) then
!         j = kj + 1
!       else
!         j = m + kj + 1
!       endif

!       ! SKIP NYQUIST FREQUENCY
!       do ki = 0, m/2-1
!         i = ki + 1

!         if ((ki == 0) .and. (kj == 0)) cycle

!         k_para = real(ki)
!         k_perp = real(kj)
!         kmod = sqrt(k_para**2 + k_perp**2)

!         ! GS95
!         if (k_perp > 0.) then
!          E_coeff = k_perp**(-7./3.)*exp(-k_para/k_perp**(2./3.))  ! 2D
!          !E_coeff = k_perp**(-10./3.)*exp(-k_para/k_perp**(2./3.))  ! 3D
!         else
!          E_coeff = 0.
!         endif

!         ! Kolmogorov
!         !E_coeff = kmod**(-8./3.) ! 2D
!         !E_coeff = kmod**(-11./3.) ! 3D

!         ! Kolmogorov + anisotropy
!         ! E_coeff = (k_para**2 + (k_perp/(1.d0 + anis))**2)**(-8./6.) ! 2D
!         !E_coeff = (k_para**2 + (k_perp/(1.d0 + anis))**2)**(-11./6.) ! 3D

!         if (kmod > kmax) then
!           phi0k(i,j) = (0.d0, 0.d0)
!         else
!         ! sort random phase
!           call random_number(ph)
!           ph = ph*twopi
!           phi0k(i,j) = sqrt(E_coeff)*(cos(ph) + (0., 1.)*sin(ph))
!         endif

!       enddo ! ki
!     enddo ! kj

!     ! execute inverse DFT
! #ifdef DP
!     call fftw_execute_dft_c2r(plan_phi0, phi0k, phi0)
! #else
!     call fftwf_execute_dft_c2r(plan_phi0, phi0k, phi0)
! #endif

!     ! normalise FFT
!     ! phi0(:,:) = phi0(:,:)/real(m*m)

!     ! destroy plan
! #ifdef DP
!     call fftw_destroy_plan(plan_phi0)
! #else
!     call fftwf_destroy_plan(plan_phi0)
! #endif

!   deallocate (phi0k)

  print*, '* Writing file: phi_0'

  file_out = trim(data_dir) // 'PHI0.DAT'
  open(unit=400, file=trim(file_out), form='formatted', status='replace', action='write')
  do i = 1, n
    do j = 1, n
      write(400,'(3(es24.16, 1x))') real(i - 1)*h, real(j - 1)*h, phi0(i,j)
    enddo
    write(400,*)
  enddo
  close(400)


  file_out = 'PHI0.BIN'
  lun = 701
  file_out = trim(data_dir) // '/' // file_out
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) phi0(:,:)
  close(lun)


  !----------------------------------------------------------------
  ! calculate power spectrum of phi0 and write to file
  !----------------------------------------------------------------
  lsum_power = .false.

  tmp3d(:,:,1) = phi0(:,:)
  call power_kk(tmp3d, ps_k, ps_kb, ps_kt, ps_kk, lsum_power, lx, ly, lz)
  call write_power_spectra('PHI0.DAT', 400)


  ! ------------------------------------------------------------------------
  ! remap scalar field phi0 into phi and write to file
  ! ------------------------------------------------------------------------
  allocate (phi(n,n))

  do j = 1, n
    do i = 1, n

      rxp = rx0(i,j)
      ryp = ry0(i,j)

      ! gives account of the periodicity
      rxp = mod(rxp, twopi)
      ryp = mod(ryp, twopi)
      if (rxp < 0.) rxp = rxp + twopi
      if (ryp < 0.) ryp = ryp + twopi

      ! calculate the left lower corner grid point indices
      ii = floor(rxp/h) + 1
      jj = floor(ryp/h) + 1

      ! calculate the right upper corner grid point indices
      iip1 = ii + 1
      jjp1 = jj + 1

      if (iip1 > n) iip1 = 2
      if (jjp1 > n) jjp1 = 2

      ! calculate linear weigths for the interpolation
      wx1 = mod(rxp, h)/h
      wy1 = mod(ryp, h)/h

      wx0 = 1.d0 - wx1
      wy0 = 1.d0 - wy1

      ! perform bilinear interpolation
      phi(i,j) =   wx0*wy0*phi0(ii  , jj  ) &
                 + wx1*wy0*phi0(iip1, jj  ) &
                 + wx0*wy1*phi0(ii  , jjp1) &
                 + wx1*wy1*phi0(iip1, jjp1)

    enddo
  enddo

  !print*, phi(23,:)

  file_out = trim(data_dir) // 'PHI.DAT'
  open(unit=400, file=trim(file_out), form='formatted', status='replace', action='write')
  do i = 1, n
    do j = 1, n
      write(400,'(3(es24.16, 1x))') real(i - 1)*h, real(j - 1)*h, phi(i,j)
    enddo
    write(400,*)
  enddo
  close(400)

  file_out = 'PHI.BIN'
  lun = 701
  file_out = trim(data_dir) // '/' // file_out
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write', access='stream')
    write(lun) phi(:,:)
  close(lun)

  !----------------------------------------------------------------
  ! calculate power spectrum of phi and write to file
  !----------------------------------------------------------------
  lsum_power = .false.

  tmp3d(:,:,1) = phi(:,:)
  call power_kk(tmp3d, ps_k, ps_kb, ps_kt, ps_kk, lsum_power, lx, ly, lz)
  call write_power_spectra('PHI.DAT', 400)


  ! ------------------------------------------------------------------------
  ! data for Reinaldo's structure function code 
  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------
  ! ATTENTION: SAVE FIELD IN DN (TEST-PURPOSE ONLY!)
  ! ------------------------------------------------------------------------

  lun = 700
  file_out = trim(data_dir) // 'DN0.BIN'
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write')
  write(lun) n, n, 1

  tmp3d(:,:,1) = phi(:,:)

  write(lun) tmp3d(:,:,:)
  write(lun) time, x, y, z, dx, dy, dz
  close(lun)


  lun = 701
  file_out = trim(data_dir) // 'BB0.BIN'
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write')
  write(lun) n, n, 1

  tmp3d(:,:,1) = bx(:,:)
  write(lun) tmp3d(:,:,:)

  tmp3d(:,:,1) = by(:,:)
  write(lun) tmp3d(:,:,:)
  ! only 2D atm, hence bz = 0
  tmp3d(:,:,1) = 0.
  write(lun) tmp3d(:,:,:)

  write(lun) time, x, y, z, dx, dy, dz
  close(lun)


  lun = 702
  file_out = trim(data_dir) // 'VV0.BIN'
  open(unit=lun, file=trim(file_out), form='unformatted', status='replace', action='write')
  write(lun) n, n, 1
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
    deallocate (mgrid(l)%dbx)
    deallocate (mgrid(l)%dby)
    deallocate (mgrid(l)%etz)
    deallocate (mgrid(l)%drx)
    deallocate (mgrid(l)%dry)
  enddo
  deallocate (mgrid)

  deallocate (bx)
  deallocate (by)
  deallocate (rx0)
  deallocate (ry0)
  deallocate (drx)
  deallocate (dry)

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

  ! stop

  contains

    !----------------------------------------------------------------
    ! calculate power spectrum of phi and write to file
    !----------------------------------------------------------------
    subroutine write_power_spectra(partial_file_name, lun)

      implicit none

      character(len=*) :: partial_file_name
      integer :: lun

      ! local dummy variables
      integer :: k, kb, kt
      character(len=400) :: file_name


      file_name = trim(data_dir) // 'power_spectra/PS_K_' // trim(partial_file_name)
      open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
      write(lun, "(a1, a4, 1x, 1(a24, 1x))") '#', 'k', &
        '|phi(k)|^2'
      do k = 1, nk
        write(lun, '(i5, 1x, 1(es24.16, 1x))') k-1, &
          ps_k(k)
      enddo
      close(lun)

      file_name = trim(data_dir) // 'power_spectra/PS_KB_' // trim(partial_file_name)
      open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
      write(lun, "(a1, a4, 1x, 1(a24, 1x))") '#', 'kx', &
        '|phi(kx)|^2'
      do k = 1, nkb
        write(lun, '(i5, 1x, 1(es24.16, 1x))') k-1, &
          ps_kb(k)
      enddo
      close(lun)

      file_name = trim(data_dir) // 'power_spectra/PS_KT_' // trim(partial_file_name)
      open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
      write(lun, "(a1, a4, 1x, 1(a24, 1x))") '#', 'kyz', &
        '|phi(kyz)|^2'
      do k = 1, nkt
        write(lun, '(i5, 1x, 1(es24.16, 1x))') k-1, &
          ps_kt(k)
      enddo
      close(lun)

      file_name = trim(data_dir) // 'power_spectra/PS_KK_' // trim(partial_file_name)
      open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
      write(lun, "(a1, a4, 1x, a5, 1x, 1(a24, 1x))") '#', 'kx', 'kyz', &
        '|phi(kx, kyz)|^2'
      do kt = 1, nkt
        do kb = 1, nkb
          write(lun, '(2(i5, 1x), 1(es24.16, 1x))') kb-1, kt-1, &
            ps_kk(kb,kt)
        enddo
        write(lun,*)
      enddo
      close(lun)

      return
    end subroutine write_power_spectra

end program main
