program main

  use spectrum, only: power_kk, dft_filter
  use, intrinsic :: iso_c_binding

  implicit none


  ! ------------------------------------------------------------------------
  ! include FFTW
  ! ------------------------------------------------------------------------
  include 'fftw3.f03'


  ! ------------------------------------------------------------------------
  ! define fixed parameters
  ! ------------------------------------------------------------------------
#ifdef DP
  integer, parameter :: sp = kind(1.d0)
#else
  integer, parameter :: sp = kind(1.0)
#endif
  integer, parameter :: dp = kind(1.d0)

  real(dp), parameter :: pi = acos(-1.d0)
  real(dp), parameter :: twopi = 2.d0*pi


  ! ------------------------------------------------------------------------
  ! define and initialize problem parameters
  ! ------------------------------------------------------------------------
integer :: ngrids = 7
integer :: nblocks_min = 1       
integer :: nblocks_max = 1!1024
integer :: nwaves_block = 4
real(sp) :: bx0 = 1.
real(sp) :: by0 = 0.
real(sp) :: bz0 = 0.
real(sp) :: rho0_avg = 1.
real(sp) :: kmod_start = 1.
real(sp) :: kmod_end = 256.
real(sp) :: anis = 1.
  real(sp) :: lx = twopi
  real(sp) :: ly = twopi
  real(sp) :: lz = twopi


  ! ------------------------------------------------------------------------
  ! define variables
  ! ------------------------------------------------------------------------
  integer :: n

  real(sp), dimension(:,:), allocatable :: rho_f
  real(sp), dimension(:,:), allocatable :: vx_f, vy_f, vz_f
  real(sp), dimension(:,:), allocatable :: bx_f, by_f, bz_f
  real(sp), dimension(:), allocatable :: x, y, z
  real(sp) :: time, dx, dy, dz

  real(sp), dimension(:,:), allocatable :: phi
  real(sp), dimension(:,:), allocatable :: field

  type(C_PTR) :: plan
  real(sp), dimension(:,:), allocatable :: f
  complex(sp), dimension(:,:), allocatable :: fk
  complex(sp), dimension(:,:), allocatable :: phik


  ! ------------------------------------------------------------------------
  ! define auxiliary variables
  ! ------------------------------------------------------------------------
  integer :: nblocks, m
  real(sp) :: kmod_start_block, kmod_end_block
  real(sp) :: kmod
  real(sp) :: b(2), b2
  real(sp) :: b1(2), ph1, ph2, aux1, aux2

  integer :: nk, nkb, nkt
  real(dp), dimension(:), allocatable :: ps_k
  real(dp), dimension(:), allocatable :: ps_kb
  real(dp), dimension(:), allocatable :: ps_kt
  real(dp), dimension(:,:), allocatable :: ps_kk
  real(dp), dimension(:,:,:), allocatable :: tmp3d

  logical :: lsum_power

  real(sp) :: k_para, k_perp, E_coeff, ph

  real(sp) :: fdummy
  integer, dimension(:), allocatable :: shift

  character(len=400) :: file_name
  character(len=400) :: data_dir


  ! ------------------------------------------------------------------------
  ! define dummy variables
  ! ------------------------------------------------------------------------
  integer :: lun
  integer :: i, j, k
  integer :: ki, kj
  integer :: block_i, block_j
  integer :: ib, jb
  integer :: i_waveset
  character(len=3) :: iblock_string
  character(len=400) :: cmd

 
  ! ------------------------------------------------------------------------
  ! specify folder for output data
  ! ------------------------------------------------------------------------
  data_dir = './128run'
  cmd = 'mkdir -p ' // trim(data_dir)
  call system(cmd)


  ! ------------------------------------------------------------------------
  ! calculate grid parameters
  ! ------------------------------------------------------------------------
  n = 2**(ngrids)

  nk = n/2 + 1
  nkb = n/2 + 1
  nkt = n/2 + 1


  ! ------------------------------------------------------------------------
  ! allocate arrays
  ! ------------------------------------------------------------------------
  allocate (rho_f(n,n))
  allocate (vx_f(n,n))
  allocate (vy_f(n,n))
  allocate (vz_f(n,n))
  allocate (bx_f(n,n))
  allocate (by_f(n,n))
  allocate (bz_f(n,n))

  allocate(x(n))
  allocate(y(n))
  allocate(z(1))

  allocate (shift(n))
  allocate (tmp3d(n,n,1))

  allocate (phi(n,n))
  allocate (field(n,n))

  allocate (ps_k(nk))
  allocate (ps_kb(nkb))
  allocate (ps_kt(nkt))
  allocate (ps_kk(nkb,nkt))


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
  ! set initial background density, velocity, and magnetic field
  ! ------------------------------------------------------------------------
  rho_f(:,:) = rho0_avg

  vx_f(:,:) = 0.
  vy_f(:,:) = 0.
  vz_f(:,:) = 0.

  bx_f(:,:) = bx0
  by_f(:,:) = by0
  bz_f(:,:) = bz0

  ! build an monochromatic field
  ! do j = 1, n
  !  do i = 1, n
  !    by_f(i,j) = by_f(i,j) + 0.5*sin(2.0*x(i))
  !  enddo
  ! enddo

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

  !         bx_f(i,j) = bx_f(i,j) + aux1
  !         by_f(i,j) = by_f(i,j) + aux2

  !       enddo
  !     enddo

  !   enddo ! ki
  ! enddo ! kj


  ! ------------------------------------------------------------------------
  ! set initial scalar field
  ! ------------------------------------------------------------------------
  field(:,:) = 0.


! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! START OF MAIN CALCULATION
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
!   - loop over wave sets (every set contains same number of waves)
!     - loop over blocks in x and y direction (number of blocks increases for larger wavenumbers)
!       - loop over grid points of an individual block 

  ! initial parameters
  i_waveset = 0
  nblocks = nblocks_min

  print*, 'Start wavenumber  --  ', 'End wavenumber  --  ', 'Num blocks  '


  ! ------------------------------------------------------------------------
  ! loop over sets of waves: k inside [kmod_start, kmod_end)
  ! ------------------------------------------------------------------------
  loop_iwaveset: do

    ! update i_waveset, kmod_start_block
    i_waveset = i_waveset + 1

    if (i_waveset == 1) then
      kmod_start_block = kmod_start
    else
      ! use kmod_end_block from previous i_waveset
      kmod_start_block = kmod_end_block
    endif


    ! check exit condition; skip nyquist frequency
    if ((kmod_start_block >= kmod_end) .or. (kmod_start_block >= real(nk - 1))) then
      exit loop_iwaveset
    endif


    ! number of waves insithe the set: logarithm spacing is a good option here
    kmod_end_block = min((kmod_start_block + kmod_start_block), real(nk - 1))


    ! number of blocks: use the maximum number of blocks 
    ! but trying to make sure that at least ~ NWAVES_BLOCK waves 
    ! (of minimum wavenumber) fit inside each block
    loop_nblocks: do
      if (nblocks >= nblocks_max) then
        nblocks = nblocks_max
        exit loop_nblocks
      endif
      if (nint(kmod_start_block/real(nblocks)) <= nwaves_block) then
        ! if less than NWAVES_BLOCK waves (of minimum wavenumber) fit inside the boxsize 
        ! from previous i_waveset, do not change nblocks
        exit loop_nblocks
      endif
      if (nint(kmod_start_block/real(2*nblocks)) >= nwaves_block) then
        ! double nblocks if at least NWAVES_BLOCK waves of each wavenumber will fit inside the new boxsize 
        ! otherwise, do not change nblocks
        nblocks = 2*nblocks
      else
        exit loop_nblocks
      endif
    enddo loop_nblocks


    ! perform circular shift on the elements of matrices 
    ! rho_f, vx_f, vy_f, vz_f, bx_f, by_f, bz_f, field 
    ! in order to reduce cummulative features at the blocks borders
     
    ! i-direction
    call random_number(fdummy)
    fdummy = fdummy*real(n/nblocks)
    shift(:) = nint(fdummy)

    rho_f = cshift(rho_f, SHIFT=shift, dim=1)
    vx_f = cshift(vx_f, SHIFT=shift, dim=1)
    vy_f = cshift(vy_f, SHIFT=shift, dim=1)
    vz_f = cshift(vz_f, SHIFT=shift, dim=1)
    bx_f = cshift(bx_f, SHIFT=shift, dim=1)
    by_f = cshift(by_f, SHIFT=shift, dim=1)
    bz_f = cshift(bz_f, SHIFT=shift, dim=1)
    field = cshift(field, SHIFT=shift, dim=1)

    ! j-direction
    call random_number(fdummy)
    fdummy = fdummy*real(n/nblocks)
    shift(:) = nint(fdummy)

    rho_f = cshift(rho_f, SHIFT=shift, dim=2)
    vx_f = cshift(vx_f, SHIFT=shift, dim=2)
    vy_f = cshift(vy_f, SHIFT=shift, dim=2)
    vz_f = cshift(vz_f, SHIFT=shift, dim=2)
    bx_f = cshift(bx_f, SHIFT=shift, dim=2)
    by_f = cshift(by_f, SHIFT=shift, dim=2)
    bz_f = cshift(bz_f, SHIFT=shift, dim=2)
    field = cshift(field, SHIFT=shift, dim=2)


    if (i_waveset < 10) write (iblock_string, "(A1, I1)") '_', i_waveset
    if (i_waveset >= 10) write (iblock_string, "(A1, I2)") '_', i_waveset 


    ! write current local background magnetic field to file - local background includes waves added during previous interation
    file_name = trim(data_dir) // '/' // 'MAGNETIC0' // trim(iblock_string) // '.DAT'
    open(unit=400, file=trim(file_name), form='formatted', status='replace', action='write')
    do i = 1, n
      do j = 1, n
        write(400,'(5(es24.16, 1x))') x(i), y(j), bx_f(i,j), by_f(i,j), bz_f(i,j)
      enddo
      write(400,*)
    enddo
    close(400)


    print*, kmod_start_block, kmod_end_block, nblocks, i_waveset


    ! ------------------------------------------------------------------------
    ! loop over blocks and build scalar field phi as a mosaic
    ! ------------------------------------------------------------------------
    m = n/nblocks

    allocate (phik((m/2 + 1), m))

    allocate (fk((m/2 + 1), m))
    allocate (f(m, m))

    ! create plan
#ifdef DP
    plan = fftw_plan_dft_c2r_2d (m, m, fk, f, FFTW_ESTIMATE)
#else
    plan = fftwf_plan_dft_c2r_2d (m, m, fk, f, FFTW_ESTIMATE)
#endif

    do block_j = 1, nblocks
      do block_i = 1, nblocks

        ! find the mean magnetic field inside the current block
        b(:) = 0.
        do j = 1, m
          do i = 1, m
            ib = (block_i - 1)*m + i
            jb = (block_j - 1)*m + j
            b(1) = b(1) + bx_f(ib, jb)
            b(2) = b(2) + by_f(ib, jb)
          enddo
        enddo 
        b(:) = b(:)/real(m*m)
        b2 = b(1)**2 + b(2)**2

        ! ! find the mean magnetic field inside the current block
        ! b(:) = 0.
        ! do j = 1, m
        !   ib = (block_i - 1)*m + 1
        !   jb = (block_j - 1)*m + j
        !   b(1) = b(1) + (bx_f(ib + 1, jb) + bx_f(ib + m-2, jb))/2.
        !   b(2) = b(2) + (by_f(ib + 1, jb) + by_f(ib + m-2, jb))/2.
        ! enddo 

        ! do i = 1, m
        !   ib = (block_i - 1)*m + i
        !   jb = (block_j - 1)*m + 1
        !   b(1) = b(1) + (bx_f(ib, jb + 1) + bx_f(ib, jb + m-2))/2.
        !   b(2) = b(2) + (by_f(ib, jb + 1) + by_f(ib, jb + m-2))/2.
        ! enddo 

        ! b(:) = b(:)/real(2.*m)
        ! b2 = b(1)**2 + b(2)**2

        ! ------------------------------------------------------------------------
        ! loop over every wavenumber in the set and build scalar field phi
        ! as a mosaic
        ! ------------------------------------------------------------------------
        phik(:,:) = 0.

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

            if ((ki == 0) .and. (kj == 0)) cycle

            kmod = sqrt(real(ki*nblocks)**2 + real(kj*nblocks)**2)

            ! only continue if kmod is inside interval 
            ! [kmod_start_block - kmod_start_block/2, kmod_end_block + kmod_end_block]
            ! if ( (kmod < (kmod_start_block - kmod_start_block/2)) &
            !      .or. (kmod > (kmod_end_block + kmod_end_block)) ) cycle

            if ( (kmod < (kmod_start_block)) &
                 .or. (kmod > (kmod_end_block)) ) cycle

            if (b2 > 0.) then
              k_para = abs(real(ki*nblocks)*b(1) + real(kj*nblocks)*b(2))/sqrt(b2)
            else
              k_para = 0.
            endif
            k_perp = sqrt(max((kmod**2 - k_para**2), 0.))


            ! print*, ki, kj, nblocks, kmod, k_para, k_perp

            ! GS95
            if (k_perp > 0.) then
             E_coeff = k_perp**(-7./3.)*exp(-k_para/k_perp**(2./3.))  ! 2D
             !E_coeff = k_perp**(-10./3.)*exp(-k_para/k_perp**(2./3.))  ! 3D
            else
             E_coeff = 0.
            endif

            ! Kolmogorov
            !E_coeff = kmod**(-8./3.) ! 2D
            !E_coeff = kmod**(-11./3.) ! 3D

            ! Kolmogorov + anisotropy
            ! E_coeff = (k_para**2 + (k_perp/(1.d0 + anis))**2)**(-8./6.) ! 2D
            !E_coeff = (k_para**2 + (k_perp/(1.d0 + anis))**2)**(-11./6.) ! 3D

            ! print*, kmod, k_para, k_perp

            ! sort random phase
            call random_number(ph)
            ph = ph*twopi

! if (k_para .eq. 0.) then
!   E_coeff = (k_perp**2)**(-8./6.)
! else
!   E_coeff = 0.
! endif

            phik(i,j) = sqrt(E_coeff)*(cos(ph) + (0., 1.)*sin(ph))

          enddo ! ki
        enddo ! kj

        ! execute inverse DFT
        ib = (block_i - 1)*m + 1
        jb = (block_j - 1)*m + 1
        ! attention with the normalization of the DFT
        fk(:,:) = phik(:,:)*real(nblocks)
#ifdef DP
        call fftw_execute_dft_c2r(plan, fk, f)
#else
        call fftwf_execute_dft_c2r(plan, fk, f)
#endif
        phi(ib:ib+m-1, jb:jb+m-1) = f(:,:)

        ! print*, 'End block_i'
        ! print*, ''
      enddo  ! block_i
      ! print*, 'End block_j'
      ! print*, ''
    enddo  ! block_j

    ! destroy plan
#ifdef DP
    call fftw_destroy_plan(plan)
#else
    call fftwf_destroy_plan(plan)
#endif

    deallocate (phik)

    deallocate (fk)
    deallocate (f)

    !----------------------------------------------------------------
    ! calculate power spectrum of phi and write to file
    !----------------------------------------------------------------
    lsum_power = .false.

    tmp3d(:,:,1) = phi(:,:)
    call power_kk(tmp3d, ps_k, ps_kb, ps_kt, ps_kk, lsum_power, lx, ly, lz)
    call write_power_spectra('PHI0' // trim(iblock_string) // '.DAT', 400)

    lun = 400
    file_name = trim(data_dir) // '/' // 'PHI0' // trim(iblock_string) // '.DAT'
    open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
    do i = 1, n
      do j = 1, n
        write(lun,'(3(es24.16, 1x))') x(i), y(j), phi(i,j)
      enddo
      write(lun,*)
    enddo
    close(lun)


    !----------------------------------------------------------------
    ! apply spectral filtering on phi
    !----------------------------------------------------------------
    if (nblocks > 1) then

      tmp3d(:,:,1) = phi(:,:)
      !call dft_filter(tmp3d, kmin=(kmod_start_block - 1.), kmax=(kmod_end_block + 1.), sigmak=1.d0, lx=lx, ly=ly, lz=lz)

      ! call dft_filter(tmp3d, kmin=kmod_start_block, kmax=kmod_end_block, sigmak=1.d0, lx=lx, ly=ly, lz=lz)
      call dft_filter(tmp3d, kmin=kmod_start_block, kmax=kmod_end_block, sigmak=-1.d0, lx=lx, ly=ly, lz=lz)
      phi(:,:) = tmp3d(:,:,1)

    endif


    !----------------------------------------------------------------
    ! calculate power spectrum of phi and write to file
    !----------------------------------------------------------------
    lsum_power = .false.

    tmp3d(:,:,1) = phi(:,:)
    call power_kk(tmp3d, ps_k, ps_kb, ps_kt, ps_kk, lsum_power, lx, ly, lz)
    call write_power_spectra('PHI' // trim(iblock_string) // '.DAT', 400)

    lun = 404
    file_name = trim(data_dir) // '/' // 'PHI' // trim(iblock_string) // '.DAT'
    open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
    do i = 1, n
      do j = 1, n
        write(lun,'(3(es24.16, 1x))') x(i), y(j), phi(i,j)
      enddo
      write(lun,*)
    enddo
    close(lun)


    !----------------------------------------------------------------
    ! add phi to field
    !----------------------------------------------------------------
    field(:,:) = field(:,:) + phi(:,:)


    !----------------------------------------------------------------
    ! calculate power spectrum of field and write to file
    !----------------------------------------------------------------
    lsum_power = .false.

    tmp3d(:,:,1) = field(:,:)
    call power_kk(tmp3d, ps_k, ps_kb, ps_kt, ps_kk, lsum_power, lx, ly, lz)
    call write_power_spectra('FIELD' // trim(iblock_string) // '.DAT', 400)

    lun = 405
    file_name = trim(data_dir) // '/' // 'FIELD' // trim(iblock_string) // '.DAT'
    open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
    do i = 1, n
      do j = 1, n
        write(lun,'(3(es24.16, 1x))') x(i), y(j), field(i,j)
      enddo
      write(lun,*)
    enddo
    close(lun)


  enddo loop_iwaveset  ! set of waves


  ! ------------------------------------------------------------------------
  ! ATTENTION: SAVE FIELD IN DN (TEST-PURPOSE ONLY!)
  ! ------------------------------------------------------------------------
  rho_f(:,:) = field(:,:)


  ! ------------------------------------------------------------------------
  ! data for Reinaldo's structure function code
  ! ------------------------------------------------------------------------
  lun = 700
  file_name = trim(data_dir) // '/' // 'DN0.BIN'
  open(unit=lun, file=trim(file_name), form='unformatted', status='replace', action='write')
  write(lun) n, n, 1

  tmp3d(:,:,1) = rho_f(:,:)
  write(lun) tmp3d(:,:,:)

  write(lun) time, x, y, z, dx, dy, dz
  close(lun)


  lun = 701
  file_name = trim(data_dir) // '/' // 'BB0.BIN'
  open(unit=lun, file=trim(file_name), form='unformatted', status='replace', action='write')
  write(lun) n, n, 1

  tmp3d(:,:,1) = bx_f(:,:)
  write(lun) tmp3d(:,:,:)

  tmp3d(:,:,1) = by_f(:,:)
  write(lun) tmp3d(:,:,:)

  tmp3d(:,:,1) = bz_f(:,:)
  write(lun) tmp3d(:,:,:)

  write(lun) time, x, y, z, dx, dy, dz
  close(lun)


  lun = 702
  file_name = trim(data_dir) // '/' // 'VV0.BIN'
  open(unit=lun, file=trim(file_name), form='unformatted', status='replace', action='write')
  write(lun) n, n, 1

  tmp3d(:,:,1) = vx_f(:,:)
  write(lun) tmp3d(:,:,:)

  tmp3d(:,:,1) = vy_f(:,:)
  write(lun) tmp3d(:,:,:)

  tmp3d(:,:,1) = vz_f(:,:)
  write(lun) tmp3d(:,:,:)

  write(lun) time, x, y, z, dx, dy, dz
  close(lun)


  ! ------------------------------------------------------------------------
  ! deallocate memory
  ! ------------------------------------------------------------------------
  deallocate (rho_f)
  deallocate (vx_f)
  deallocate (vy_f)
  deallocate (vz_f)
  deallocate (bx_f)
  deallocate (by_f)
  deallocate (bz_f)

  deallocate (x)
  deallocate (y)
  deallocate (z)

  deallocate (shift)
  deallocate (tmp3d)

  deallocate (phi)
  deallocate (field)

  deallocate (ps_k)
  deallocate (ps_kb)
  deallocate (ps_kt)
  deallocate (ps_kk)

  stop


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


      file_name = trim(data_dir) // '/' // 'PS_K_' // trim(partial_file_name)
      open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
      write(lun, "(a1, a4, 1x, 1(a24, 1x))") '#', 'k', &
        '|phi(k)|^2'
      do k = 1, nk
        write(lun, '(i5, 1x, 1(es24.16, 1x))') k-1, &
          ps_k(k)
      enddo
      close(lun)

      file_name = trim(data_dir) // '/' // 'PS_KB_' // trim(partial_file_name)
      open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
      write(lun, "(a1, a4, 1x, 1(a24, 1x))") '#', 'kx', &
        '|phi(kx)|^2'
      do k = 1, nkb
        write(lun, '(i5, 1x, 1(es24.16, 1x))') k-1, &
          ps_kb(k)
      enddo
      close(lun)

      file_name = trim(data_dir) // '/' // 'PS_KT_' // trim(partial_file_name)
      open(unit=lun, file=trim(file_name), form='formatted', status='replace', action='write')
      write(lun, "(a1, a4, 1x, 1(a24, 1x))") '#', 'kyz', &
        '|phi(kyz)|^2'
      do k = 1, nkt
        write(lun, '(i5, 1x, 1(es24.16, 1x))') k-1, &
          ps_kt(k)
      enddo
      close(lun)

      file_name = trim(data_dir) // '/' // 'PS_KK_' // trim(partial_file_name)
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
