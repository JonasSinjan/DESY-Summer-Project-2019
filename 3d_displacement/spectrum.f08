module spectrum

use, intrinsic :: iso_c_binding

implicit none

include 'fftw3.f03'

#ifdef DP
integer, parameter :: sp = kind(1.d0)
#else
integer, parameter :: sp = kind(1.0)
#endif
integer, parameter :: dp = kind(1.d0)

private
public :: power_kk
public :: dft_filter

contains

subroutine power_kk(field, ps_k, ps_kx, ps_kyz, ps_kxkyz, lsum_power, lx, ly, lz)

  real(sp), dimension(:,:,:), intent(inout) :: field
  real(dp), dimension(:), intent(inout) :: ps_k, ps_kx, ps_kyz
  real(dp), dimension(:,:), intent(inout) :: ps_kxkyz
  logical, intent(in) :: lsum_power
  real(dp), intent(in) :: lx, ly, lz

  type(C_PTR) :: plan
  real(sp), dimension(:,:,:), allocatable :: f
  complex(sp), dimension(:,:,:), allocatable :: fk

  integer :: nx, ny, nz
  integer :: ki, kj, kk
  integer :: nk, nkx, nkyz, nkxkyz1, nkxkyz2
  integer :: k_int, kx_int, kyz_int
  real(dp) :: kx, ky, kz, kmod
  real(dp) :: tmp


  nx = size(field, dim=1)
  ny = size(field, dim=2)
  nz = size(field, dim=3)

  nk = size(ps_k)
  nkx = size(ps_kx)
  nkyz = size(ps_kyz)
  nkxkyz1 = size(ps_kxkyz, dim=1)
  nkxkyz2 = size(ps_kxkyz, dim=2)

  allocate(f(nx, ny, nz))
  allocate(fk((nx/2 + 1), ny, nz))

#ifdef DP
  plan = fftw_plan_dft_r2c_3d(nz, ny, nx, f, fk, FFTW_ESTIMATE)
#else
  plan = fftwf_plan_dft_r2c_3d(nz, ny, nx, f, fk, FFTW_ESTIMATE)
#endif

  ! make copy of matrix (will be destroyed by the fft)
  f(:,:,:) = field(:,:,:)

#ifdef DP
  call fftw_execute_dft_r2c(plan, f, fk)
  call fftw_destroy_plan(plan)
#else
  call fftwf_execute_dft_r2c(plan, f, fk)
  call fftwf_destroy_plan(plan)
#endif

  ! normalize fk
  fk(:,:,:) = fk(:,:,:) / real(nx*ny*nz)

  if (lsum_power .eqv. .false.) then
    ! initialize ps
    ps_k(:) = 0.
    ps_kx(:) = 0.
    ps_kyz(:) = 0.
    ps_kxkyz(:,:) = 0.
  endif

  ! accumulate power
  do kk = 0, nz/2
    kz = real(kk)*(ly/lz)

    do kj = 0, ny/2
      ky = real(kj)

      do ki = 0, nx/2
        kx = real(ki)*(ly/lx)

        kmod = sqrt(kx**2 + ky**2 + kz**2)
        !k_int = nint(kmod)
        k_int = floor(kmod)

        !kx_int = floor(kx)
        kx_int = ki

        kmod = sqrt(ky**2 + kz**2)
        !kyz_int = nint(kmod)
        kyz_int = floor(kmod)

        ! component (ki,kj,kk)
        tmp =   realpart(fk(ki+1, kj+1, kk+1))**2 &
              + imagpart(fk(ki+1, kj+1, kk+1))**2

        if (kj > 0) then
          ! component (ki,-kj,kk)
          tmp = tmp + realpart(fk(ki+1, ny-kj+1, kk+1))**2 &
                    + imagpart(fk(ki+1, ny-kj+1, kk+1))**2
        endif

        if (kk > 0) then
          ! component (ki,kj,-kk)
          tmp = tmp + realpart(fk(ki+1, kj+1, nz-kk+1))**2 &
                    + imagpart(fk(ki+1, kj+1, nz-kk+1))**2
        endif

        if ((kj > 0) .and. (kk > 0)) then
          ! component (ki,-kj,-kk)
          tmp = tmp + realpart(fk(ki+1, ny-kj+1, nz-kk+1))**2 &
                    + imagpart(fk(ki+1, ny-kj+1, nz-kk+1))**2
        endif

        if (ki > 0) then
          ! component (-ki,-kj,-kk)
          tmp = 2.*tmp
        endif

        if (k_int < nk) then
          ps_k(k_int+1) = ps_k(k_int+1) + tmp
        endif

        if (kx_int < nkx) then
          ps_kx(kx_int+1) = ps_kx(kx_int+1) + tmp
        endif

        if (kyz_int < nkyz) then
          ps_kyz(kyz_int+1) = ps_kyz(kyz_int+1) + tmp
        endif

        if ((kx_int < nkxkyz1) .and. (kyz_int < nkxkyz2)) then
          ps_kxkyz(kx_int+1, kyz_int+1) = ps_kxkyz(kx_int+1, kyz_int+1) + tmp
        endif

      enddo
    enddo
  enddo

  deallocate(f)
  deallocate(fk)

  return
end subroutine power_kk


subroutine dft_filter(a, kmin, kmax, sigmak, lx, ly, lz)

  real(sp), dimension(:,:,:), intent(inout) :: a
  real(dp), intent(in) :: kmin, kmax, sigmak
  real(dp), intent(in) :: lx, ly, lz

  ! size of the array a
  integer :: nx, ny, nz

  ! dft variables
  type(C_PTR) :: plan1, plan2
  real(sp), dimension(:,:,:), allocatable :: f
  complex(sp), dimension(:,:,:), allocatable :: fk

  ! auxiliary variables
  real(dp) :: kx, ky, kz, kmod
  real(dp) :: aux, damp

  ! dummy variables
  integer :: ki, kj, kk


  nx = size(a, dim=1)
  ny = size(a, dim=2)
  nz = size(a, dim=3)

  allocate(f(nx, ny, nz))
  allocate(fk((nx/2 + 1), ny, nz))


#ifdef DP
  plan1 = fftw_plan_dft_r2c_3d(nz, ny, nx, f, fk, FFTW_ESTIMATE)
  plan2 = fftw_plan_dft_c2r_3d(nz, ny, nx, fk, f, FFTW_ESTIMATE)
#else
  plan1 = fftwf_plan_dft_r2c_3d(nz, ny, nx, f, fk, FFTW_ESTIMATE)
  plan2 = fftwf_plan_dft_c2r_3d(nz, ny, nx, fk, f, FFTW_ESTIMATE)
#endif

  f(:,:,:) = a(:,:,:)

#ifdef DP
  call fftw_execute_dft_r2c(plan1, f, fk)
#else
  call fftwf_execute_dft_r2c(plan1, f, fk)
#endif

  aux = 2.d0*sigmak**2

  do kk = 0, nz/2
    kz = real(kk)*(ly/lz)

    do kj = 0, ny/2
      ky = real(kj)

      do ki = 0, nx/2
        kx = real(ki)*(ly/lx)

        ! do not remove background component
        if ((ki==0) .and. (kj==0) .and. (kk==0)) cycle

        kmod = sqrt(kx**2 + ky**2 + kz**2)

        if ((kmod < kmin) .or. (kmod > kmax)) then

          if (sigmak > 0.) then
            if (kmod < kmin) then
              damp = exp(-(kmin - kmod)**2/aux)
            else
              damp = exp(-(kmod - kmax)**2/aux)
            endif
          else
            damp = 0.
          endif

          ! component (ki,kj,kk)
          fk(ki+1, kj+1, kk+1) = damp*fk(ki+1, kj+1, kk+1)

          if (kj > 0) then
            ! component (ki,-kj,kk)
            fk(ki+1, ny-kj+1, kk+1) = damp*fk(ki+1, ny-kj+1, kk+1)
          endif

          if (kk > 0) then
            ! component (ki,kj,-kk)
            fk(ki+1, kj+1, nz-kk+1) = damp*fk(ki+1, kj+1, nz-kk+1)
          endif

          if ((kj > 0) .and. (kk > 0)) then
            ! component (ki,-kj,-kk)
            fk(ki+1, ny-kj+1, nz-kk+1) = damp*fk(ki+1, ny-kj+1, nz-kk+1)
          endif

        endif

      enddo
    enddo
  enddo


  ! perform inverse discrete fourier transform with pack fftw3
  ! and destroy plans
#ifdef DP
  call fftw_execute_dft_c2r(plan2, fk, f)

  call fftw_destroy_plan(plan1)
  call fftw_destroy_plan(plan2)
#else
  call fftwf_execute_dft_c2r(plan2, fk, f)

  call fftwf_destroy_plan(plan1)
  call fftwf_destroy_plan(plan2)
#endif


  ! normalize dft
  a(:,:,:) = f(:,:,:) / real(nx*ny*nz)


  deallocate(f)
  deallocate(fk)

  return
end subroutine dft_filter

end module spectrum
