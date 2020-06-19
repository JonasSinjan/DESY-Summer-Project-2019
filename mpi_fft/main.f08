program main
    
    use, intrinsic :: iso_c_binding
    ! include 'fftw3.f03'
    include 'fftw3-mpi.f03'

    integer, parameter :: sp = kind(1.0)
    real(dp), parameter :: pi = acos(-1.d0)
    real(dp), parameter :: twopi = 2.d0*pi
    

    real(sp), dimension(:,:,:), allocatable :: phi0
    complex(sp), dimension(:,:,:), allocatable :: phi0k, fk
    type(C_PTR) :: cdata
    integer*8 :: dftplan
    real(sp), dimension(:,:,:), allocatable :: f
    integer(C_INTPTR_T) :: i, j, alloc_local, local_M, local_j_offset

    integer :: m, num_seed
    integer, dimension(1) :: rand_seed
    integer :: ngrids = 2
    integer(C_INTPTR_T), parameter :: n = 2**(ngrids) + 1
    real(sp) :: kmax, kmod, mag
    real(sp) :: k_para, k_perp, E_coeff, ph
    integer :: i, j, k
    integer :: ii, jj, k_k
    integer :: ki, kj, kk 



    ! n = 2**(ngrids) + 1
    m = n - 1

    call random_seed(size = num_seed)
    allocate(seed(num_seed))

    do i = 1, num_seed
        seed(i) = i*4251
    enddo
    
    call random_seed(put=seed)

    allocate (phi0(n,n,n))
    allocate (phi0k((m/2 + 1), m, m))
    allocate (fk((m/2 + 1), m, m))
    allocate (f(m, m, m))

    call MPI_init()
    call fftw_mpi_init()    

    phi0(:,:,:) = 0 ! initialise all entries to zero

    alloc_local = fftw_mpi_local_size_2d(M, L, MPI_COMM_WORLD, &
                                       local_M, local_j_offset)
    cdata = fftw_alloc_complex(alloc_local)
    call c_f_pointer(cdata, phi0, [L,local_M])

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
                ! if (kmod > kmax) then from michaels method
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
    call dfftw_free(cdata)

    deallocate (phi0k)
    deallocate (fk)
    deallocate (f)
    deallocate (phi0)

    print*, 'The loop has successfully completed'


  stop

end program main