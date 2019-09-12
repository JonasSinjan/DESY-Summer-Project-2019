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
  integer :: ngrids = 7
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
  integer :: n, num_seed
  integer, dimension(1) :: rand_seed
  integer, allocatable :: seed(:)
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
  complex(sp), dimension(:,:,:), allocatable :: phi0k,fk

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
  type(C_PTR) :: plan
  integer*8 :: dftplan
  real(sp), dimension(:,:,:), allocatable :: f
  complex(sp), dimension(:,:,:), allocatable :: bxk, byk, bzk !3d
  complex(sp), dimension(:,:,:), allocatable :: dbxk, dbyk, dbzk !3d
  complex(sp), dimension(:,:,:), allocatable :: etzk_x, etzk_y, etzk_z
  
  real(sp) :: b(3), b2 !3d
  real(sp) :: b1(3), ph1, ph2, ph3, aux1, aux2, aux3 !3d
  real(sp) :: rxp, ryp, drxp, dryp, rzp, drzp !3d
  real(sp) :: wx0, wx1, wy0, wy1, wz0, wz1
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
  data_dir = './128_mem_FFT/'

  cmd = 'mkdir -p ' // trim(data_dir)
  call system(cmd)

  ! ------------------------------------------------------------------------
  ! calculate grid parameters
  ! ------------------------------------------------------------------------
  n = 2**(ngrids) + 1
  h = twopi/real(n - 1) !box length?


  nk = n/2 + 1
  nkb = n/2 + 1
  nkt = n/2 + 1

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

  time = 0.

  tot_time = omp_get_wtime()
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
  enddo

  do l = 1, ngrids

    allocate (mgrid(l)%drx(n,n,n))
    allocate (mgrid(l)%dry(n,n,n))
    allocate (mgrid(l)%drz(n,n,n)) !3d
    
    !only allocate once
    if (l==1) then

      !set global magnetic field with perturb - only filtered afterwards
      mgrid(1)%bx(:,:,:) = bx(:,:,:) !:,:,: 3d
      mgrid(1)%by(:,:,:) = by(:,:,:)
      mgrid(1)%bz(:,:,:) = bz(:,:,:) !3d


      !allocating second set
      allocate (mgrid(1)%dbx(n,n,n))
      allocate (mgrid(1)%dby(n,n,n))
      allocate (mgrid(1)%dbz(n,n,n)) !3d
      allocate (mgrid(1)%etz_x(n,n,n))
      allocate (mgrid(1)%etz_y(n,n,n))
      allocate (mgrid(1)%etz_z(n,n,n))
    
    endif


    ! ------------------------------------------------------------------------
    ! build bx and by in each grid level, using spectral filtering
    ! ------------------------------------------------------------------------
    if (l==1) then

      do q = 2, ngrids
      
        kmax = sqrt(2.d0)*(2.d0**(ngrids - q))

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
        f(:,:,:) = mgrid(q-1)%bx(1:m,1:m,1:m) !1:m for all dimensions?
        call fftw_execute_dft_r2c(plan1, f, bxk)
        !r2c = real to complex
        f(:,:,:) = mgrid(q-1)%by(1:m,1:m,1:m)
        call fftw_execute_dft_r2c(plan1, f, byk)

        f(:,:,:) = mgrid(q-1)%bz(1:m,1:m,1:m) !3d again do i need bzk?
        call fftw_execute_dft_r2c(plan1, f, bzk)

        call fftw_destroy_plan(plan1)
#else
        f(:,:,:) = mgrid(q-1)%bx(1:m,1:m,1:m)
        call fftwf_execute_dft_r2c(plan1, f, bxk)

        f(:,:,:) = mgrid(q-1)%by(1:m,1:m,1:m)
        call fftwf_execute_dft_r2c(plan1, f, byk)

        f(:,:,:) = mgrid(q-1)%bz(1:m,1:m,1:m)
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
        print*, "performing spectral filtering to find next B_i dist."
        do kk = min((-m/2 + 1), 0), m/2
          if (kk >= 0) then
            k = kk + 1
          else
            k = m + kk + 1
          endif  
          !print*, kk
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
        print*, ngrids, l, l+1
        mgrid(q)%bx(1:m,1:m,1:m) = f(:,:,:) !inner cube

        mgrid(q)%bx(m+1,1:m,1:m) = f(1,:,:) !outer faces
        mgrid(q)%bx(1:m,m+1,1:m) = f(:,1,:)
        mgrid(q)%bx(1:m,1:m,m+1) = f(:,:,1)
        
        mgrid(q)%bx(1:m,m+1,m+1) = f(:,1,1) !outer edges
        mgrid(q)%bx(m+1,m+1,1:m) = f(1,1,:)
        mgrid(q)%bx(m+1,1:m,m+1) = f(1,:,1)
        
        mgrid(q)%bx(m+1,m+1,m+1) = f(1,1,1) !last point

        call fftw_execute_dft_c2r(plan2, byk, f)
        mgrid(q)%by(1:m,1:m,1:m) = f(:,:,:) !inner cube

        mgrid(q)%by(m+1,1:m,1:m) = f(1,:,:) !outer faces
        mgrid(q)%by(1:m,m+1,1:m) = f(:,1,:)
        mgrid(q)%by(1:m,1:m,m+1) = f(:,:,1)
        
        mgrid(q)%by(1:m,m+1,m+1) = f(:,1,1) !outer edges
        mgrid(q)%by(m+1,m+1,1:m) = f(1,1,:)
        mgrid(q)%by(m+1,1:m,m+1) = f(1,:,1)
        
        mgrid(q)%by(m+1,m+1,m+1) = f(1,1,1) !last point

        call fftw_execute_dft_c2r(plan2, bzk, f) !3d
        mgrid(q)%bz(1:m,1:m,1:m) = f(:,:,:) !inner cube

        mgrid(q)%bz(m+1,1:m,1:m) = f(1,:,:) !outer faces
        mgrid(q)%bz(1:m,m+1,1:m) = f(:,1,:)
        mgrid(q)%bz(1:m,1:m,m+1) = f(:,:,1)
        
        mgrid(q)%bz(1:m,m+1,m+1) = f(:,1,1) !outer edges
        mgrid(q)%bz(m+1,m+1,1:m) = f(1,1,:)
        mgrid(q)%bz(m+1,1:m,m+1) = f(1,:,1)
        
        mgrid(q)%bz(m+1,m+1,m+1) = f(1,1,1) !last point

        call fftw_destroy_plan(plan2)
#else
        call fftw_execute_dft_c2r(plan2, bxk, f)
        !c2r = complex to real
        !need to add the slices for periodic boundaries
        mgrid(q)%bx(1:m,1:m,1:m) = f(:,:,:) !inner cube

        mgrid(q)%bx(m+1,1:m,1:m) = f(1,:,:) !outer faces
        mgrid(q)%bx(1:m,m+1,1:m) = f(:,1,:)
        mgrid(q)%bx(1:m,1:m,m+1) = f(:,:,1)
        
        mgrid(q)%bx(1:m,m+1,m+1) = f(:,1,1) !outer edges
        mgrid(q)%bx(m+1,m+1,1:m) = f(1,1,:)
        mgrid(q)%bx(m+1,1:m,m+1) = f(1,:,1)
        
        mgrid(q)%bx(m+1,m+1,m+1) = f(1,1,1) !last point

        call fftw_execute_dft_c2r(plan2, byk, f)
        mgrid(q)%by(1:m,1:m,1:m) = f(:,:,:) !inner cube

        mgrid(q)%by(m+1,1:m,1:m) = f(1,:,:) !outer faces
        mgrid(q)%by(1:m,m+1,1:m) = f(:,1,:)
        mgrid(q)%by(1:m,1:m,m+1) = f(:,:,1)
      
        mgrid(q)%by(1:m,m+1,m+1) = f(:,1,1) !outer edges
        mgrid(q)%by(m+1,m+1,1:m) = f(1,1,:)
        mgrid(q)%by(m+1,1:m,m+1) = f(1,:,1)
        
        mgrid(q)%by(m+1,m+1,m+1) = f(1,1,1) !last point

        call fftw_execute_dft_c2r(plan2, bzk, f) !3d
        mgrid(q)%bz(1:m,1:m,1:m) = f(:,:,:) !inner cube

        mgrid(q)%bz(m+1,1:m,1:m) = f(1,:,:) !outer faces
        mgrid(q)%bz(1:m,m+1,1:m) = f(:,1,:)
        mgrid(q)%bz(1:m,1:m,m+1) = f(:,:,1)
        
        mgrid(q)%bz(1:m,m+1,m+1) = f(:,1,1) !outer edges
        mgrid(q)%bz(m+1,m+1,1:m) = f(1,1,:)
        mgrid(q)%bz(m+1,1:m,m+1) = f(1,:,1)
        
        mgrid(q)%bz(m+1,m+1,m+1) = f(1,1,1) !last point
        
        call fftw_destroy_plan(plan2)

#endif

        ! deallocate auxiliary matrices
        deallocate (f)
        deallocate (bxk)
        deallocate (byk)
        deallocate (bzk)

      enddo

    endif

  
    ! ------------------------------------------------------------------------
    ! calculate vector field db in each grid level
    ! ------------------------------------------------------------------------
    print*, "calculating db"
    do k = 1, n !3d
      !print*, k
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

          mgrid(1)%dbx(i,j,k) = b(1) - mgrid(l)%bx(i,j,k)
          mgrid(1)%dby(i,j,k) = b(2) - mgrid(l)%by(i,j,k)
          mgrid(1)%dbz(i,j,k) = b(3) - mgrid(l)%bz(i,j,k) !3d


        enddo
      enddo
    enddo

  !no longer need l+1
  


! ------------------------------------------------------------------------
! calculate etz (only non-null component) from equation 
! curl(et) = db in each grid level
! ------------------------------------------------------------------------

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
    f(:,:,:) = mgrid(1)%dbx(1:m,1:m,1:m)
    call fftw_execute_dft_r2c(plan1, f, dbxk)

    f(:,:,:) = mgrid(1)%dby(1:m,1:m,1:m)
    call fftw_execute_dft_r2c(plan1, f, dbyk)

    f(:,:,:) = mgrid(1)%dbz(1:m,1:m,1:m) !3d
    call fftw_execute_dft_r2c(plan1, f, dbzk)

    call fftw_destroy_plan(plan1)
#else
    f(:,:,:) = mgrid(1)%dbx(1:m,1:m,1:m)
    call fftwf_execute_dft_r2c(plan1, f, dbxk)

    f(:,:,:) = mgrid(1)%dby(1:m,1:m,1:m)
    call fftwf_execute_dft_r2c(plan1, f, dbyk)

    f(:,:,:) = mgrid(1)%dbz(1:m,1:m,1:m) !3d
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
    print*, "calculating etz components"
    do kk = min((-m/2 + 1), 0), m/2
      if (kk >= 0) then
        k = kk + 1
      else
        k = m + kk + 1
      endif
      !print*, kk    
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

          !testing requirement dot product is zero
          test_var = abs(etzk_x(i,j,k))*ki + abs(etzk_y(i,j,k))*kj + abs(etzk_z(i,j,k))*kk
          if (abs(test_var)>0.00001) then
            print*, test_var
          endif
        enddo
      enddo
    enddo

    deallocate (dbxk)
    deallocate (dbyk)
    deallocate (dbzk)


    ! transform each component of etzk to real space, destroy plan2
#ifdef DP
    call fftw_execute_dft_c2r(plan2, etzk_x, f)
    mgrid(1)%etz_x(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(1)%etz_x(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(1)%etz_x(1:m,m+1,1:m) = f(:,1,:)
    mgrid(1)%etz_x(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(1)%etz_x(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(1)%etz_x(m+1,m+1,1:m) = f(1,1,:)
    mgrid(1)%etz_x(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(1)%etz_x(m+1,m+1,m+1) = f(1,1,1) !last point

    call fftw_execute_dft_c2r(plan2, etzk_y, f)
    mgrid(1)%etz_y(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(1)%etz_y(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(1)%etz_y(1:m,m+1,1:m) = f(:,1,:)
    mgrid(1)%etz_y(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(1)%etz_y(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(1)%etz_y(m+1,m+1,1:m) = f(1,1,:)
    mgrid(1)%etz_y(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(1)%etz_y(m+1,m+1,m+1) = f(1,1,1) !last point

    call fftw_execute_dft_c2r(plan2, etzk_z, f)
    mgrid(1)%etz_z(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(1)%etz_z(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(1)%etz_z(1:m,m+1,1:m) = f(:,1,:)
    mgrid(1)%etz_z(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(1)%etz_z(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(1)%etz_z(m+1,m+1,1:m) = f(1,1,:)
    mgrid(1)%etz_z(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(1)%etz_z(m+1,m+1,m+1) = f(1,1,1) !last point
    call fftw_destroy_plan(plan2)
#else
    call fftw_execute_dft_c2r(plan2, etzk_x, f)
    mgrid(1)%etz_x(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(1)%etz_x(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(1)%etz_x(1:m,m+1,1:m) = f(:,1,:)
    mgrid(1)%etz_x(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(1)%etz_x(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(1)%etz_x(m+1,m+1,1:m) = f(1,1,:)
    mgrid(1)%etz_x(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(1)%etz_x(m+1,m+1,m+1) = f(1,1,1) !last point
    
    call fftw_execute_dft_c2r(plan2, etzk_y, f)
    mgrid(1)%etz_y(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(1)%etz_y(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(1)%etz_y(1:m,m+1,1:m) = f(:,1,:)
    mgrid(1)%etz_y(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(1)%etz_y(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(1)%etz_y(m+1,m+1,1:m) = f(1,1,:)
    mgrid(1)%etz_y(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(1)%etz_y(m+1,m+1,m+1) = f(1,1,1) !last point

    call fftw_execute_dft_c2r(plan2, etzk_z, f)
    mgrid(1)%etz_z(1:m,1:m,1:m) = f(:,:,:) !inner cube

    mgrid(1)%etz_z(m+1,1:m,1:m) = f(1,:,:) !outer faces
    mgrid(1)%etz_z(1:m,m+1,1:m) = f(:,1,:)
    mgrid(1)%etz_z(1:m,1:m,m+1) = f(:,:,1)
    
    mgrid(1)%etz_z(1:m,m+1,m+1) = f(:,1,1) !outer edges
    mgrid(1)%etz_z(m+1,m+1,1:m) = f(1,1,:)
    mgrid(1)%etz_z(m+1,1:m,m+1) = f(1,:,1)
    
    mgrid(1)%etz_z(m+1,m+1,m+1) = f(1,1,1) !last point
    
    call fftwf_destroy_plan(plan2)
#endif

    ! deallocate auxiliary matrices
    deallocate (f)
    deallocate (etzk_x)
    deallocate (etzk_y)
    deallocate (etzk_z)




  ! ------------------------------------------------------------------------
  ! calculate displacement field dr from equation 
  ! et = dr x b in each grid level
  ! (need to assume dr . b = 0)
  ! 3d so cross product
  ! ------------------------------------------------------------------------
    print*, "calculating dr"

    do k = 1, n
      do j = 1, n
        do i = 1, n

          b(1) = mgrid(l)%bx(i,j,k) + mgrid(l)%dbx(i,j,k)
          b(2) = mgrid(l)%by(i,j,k) + mgrid(l)%dby(i,j,k)
          b(3) = mgrid(l)%bz(i,j,k) + mgrid(l)%dbz(i,j,k)

          b2 = b(1)**2 + b(2)**2+b(3)**2 !magnitude of 3d B field

          mgrid(l)%drx(i,j,k) = (b(2)*mgrid(1)%etz_z(i,j,k) - b(3)*mgrid(1)%etz_y(i,j,k))/b2
          mgrid(l)%dry(i,j,k) = (b(3)*mgrid(1)%etz_x(i,j,k) - b(1)*mgrid(1)%etz_z(i,j,k))/b2
          mgrid(l)%drz(i,j,k) = (b(1)*mgrid(1)%etz_y(i,j,k) - b(2)*mgrid(1)%etz_x(i,j,k))/b2

        enddo
      enddo
    enddo

    
    !only deallocate in last loop

    !need to keep drx for all grids to further calculations

    if (l==ngrids) then 
      
      deallocate (mgrid(1)%dbx)
      deallocate (mgrid(1)%dby)
      deallocate (mgrid(1)%dbz)
      deallocate (mgrid(1)%etz_x)
      deallocate (mgrid(1)%etz_y)
      deallocate (mgrid(1)%etz_z)

      do q = 1, ngrids
        deallocate (mgrid(q)%bx)
        deallocate (mgrid(q)%by)
        deallocate (mgrid(q)%bz)
      enddo

    endif
  
  enddo


! this is where the big loop stops - from here only need drx,dry,drz for each mgrid


! ------------------------------------------------------------------------
! calculate origin points rx0, ry0 and rz0
! ------------------------------------------------------------------------
  allocate (rx0(n,n,n))
  allocate (ry0(n,n,n))
  allocate (rz0(n,n,n))

  print*, "calculating origin (rx,ry,rz)"
  do k = 1, n
    !print*, k
    do j = 1, n
      !print*, j
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

          deallocate (mgrid(l)%drx)
          deallocate (mgrid(l)%dry)
          deallocate (mgrid(l)%drz)

        enddo

      enddo
    enddo
  enddo

  deallocate (mgrid)
  ! ------------------------------------------------------------------------
  ! calculate total displacement fields drx and dry and write to file - UNSURE FOR 3D
  ! ------------------------------------------------------------------------
  allocate (drx(n,n,n))
  allocate (dry(n,n,n))
  allocate (drz(n,n,n))
  print*, "calculating total displacement"
  do k = 1, n
    do j = 1, n
      do i = 1, n
        drx(i,j,k) = rx0(i,j,k) - real(i - 1)*h
        dry(i,j,k) = ry0(i,j,k) - real(j - 1)*h
        drz(i,j,k) = rz0(i,j,k) - real(k - 1)*h
      enddo
    enddo
  enddo
 

  ! ------------------------------------------------------------------------
  ! build scalar field phi0 and write to file
  ! ------------------------------------------------------------------------
  m = n - 1

  allocate (phi0(n,n,n)) 
  
  allocate (phi0k((m/2 + 1), m, m))
  allocate (fk((m/2 + 1), m, m))
  allocate (f(m, m, m))

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
  
  print*, omp_get_max_threads(), "calculating phi0k"
  wtime = omp_get_wtime()

  call dfftw_plan_dft_c2r_3d(dftplan, m,m,m, fk, f, FFTW_ESTIMATE)

  phi0k(:,:,:) = 0
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

        kmod = sqrt(real(ki)**2 + real(kj)**2 + real(kk)**2)

        k_para = abs(ki)
        k_perp = sqrt(max((kmod**2 - k_para**2), 0.))

        ! GS95
        if (k_perp > 0.) then
          E_coeff = k_perp**(-10./3.)*exp(-k_para/k_perp**(2./3.))  ! 3D
        else
          E_coeff = 0.
        endif

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
  !lsum_power = .false.

  !tmp3d(:,:,:) = phi0(:,:,:)
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
        rz0(i,j,k) = rzp

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
  !lsum_power = .false.

  !tmp3d(:,:,:) = phi(:,:,:)
  !call power_kk(tmp3d, ps_k, ps_kb, ps_kt, ps_kk, lsum_power, lx, ly, lz)
  !call write_power_spectra('PHI.DAT', 400) - no DAT files yet


  ! ------------------------------------------------------------------------
  ! data for Reinaldo's structure function code 
  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------
  ! ATTENTION: SAVE FIELD IN DN (TEST-PURPOSE ONLY!)
  ! ------------------------------------------------------------------------


  tot_time = omp_get_wtime() - tot_time
  print *, tot_time, "Total Time"
  ! ------------------------------------------------------------------------
  ! deallocate memory
  ! ------------------------------------------------------------------------
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

  deallocate (phi0)
  deallocate (phi)


  stop


end program main
