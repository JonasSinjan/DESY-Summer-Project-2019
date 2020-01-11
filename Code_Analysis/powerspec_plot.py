"""
TODO:
    1.  check which dimension is kperp and kpara
    2.  integrate up kperp and kpara - get spectrum in each (just by adding up each for constant other value)
    3.  plot them against each other
    4.  do this for multiple dimensions/methods
    5.  save the plots in new folder
FIXME:
    1. for some reason 129 points in phi0 file
    2. have to transpose - now get same numbers
    3. but plots look quite odd right now
"""

import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib import ticker
from matplotlib import gridspec
from k_perp_3d_spectrum import k_perp_calculator
import scipy.interpolate as spint

def plot_power2d(dir1, dir2, n):
    # 2D
    nx = n + 1 #displacement
    ny = n + 1

    # Reading in PHI0 binary file
    filename = dir1 + 'PHI0' + '.BIN'  # PHI for displacement, RHO for squares
    # print(filename)
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    tmp = np.reshape(abx, (nx, ny))
    phi0 = tmp.transpose()

    # Reading in PHI0 binary file
    filename = dir1 + 'PHI' + '.BIN'
    # print(filename)
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

    tmp2 = np.reshape(abx, (nx, ny))
    phi = tmp2.transpose()

    # print(phi0[0, :], len(phi0[0, :]))

    # FFT
    phi0k_noshift = fft.fft2(phi0)
    phik_noshift = fft.fft2(phi)

    # Perform Shift
    phi0k = fft.fftshift(phi0k_noshift)
    phik = fft.fftshift(phik_noshift)

    # plotting phik at particular k_perp to see exponential drop off at centre
    # var = abs(phi0k) ** 2
    # plt.figure()
    # plt.plot(np.log(var[:, 2]), label='2')
    # plt.plot(np.log(var[:, 80]), label='80')
    # plt.legend()
    # plt.title('Phi0k at particular k_perp to see exp. drop from centre')
    #plt.show()

    # test array for index
    # tmp_arr = np.zeros((20, 20))
    # for i in range(20):
    #     tmp_arr[i, :] = i

    # print(tmp_arr)

    # fig = plt.figure()
    # fig = plt.figure(figsize=(5.0, 5.0))
    # gs = gridspec.GridSpec(2, 1, hspace=0.5, wspace=0.2)

    # ax0 = plt.subplot(gs[0], aspect='equal')

    # plt.imshow(np.log(abs(phi0k) ** 2), cmap='seismic',
    #            interpolation='nearest', origin='lower')
    # # plt.imshow(tmp_arr)
    # fig = plt.gcf()  # get the current figure
    # plt.clim()  # clamp the color limits
    # plt.colorbar()
    # plt.ylabel('K_parallel')
    # plt.xlabel('K_perp')
    # plt.title('FFT of Phi0 2D')

    # ax1 = plt.subplot(gs[1], aspect='equal')

    # plt.imshow(np.log(abs(phik) ** 2), cmap='seismic',
    #            interpolation='nearest', origin='lower')
    # fig = plt.gcf()  # get current figure
    # plt.clim()  # clamp the color limits
    # plt.colorbar()
    # plt.ylabel('K_parallel')
    # plt.xlabel('K_perp')
    # plt.title('FFT of Phi 2D')
    #plt.show()

    perp_spectrum = np.zeros(nx)
    para_spectrum = np.zeros(nx)
    for i in range(nx):
        perp_spectrum[i] = np.sum(abs(phik[:, i]) ** 2)
        para_spectrum[i] = np.sum(abs(phik[i, :]) ** 2)

    para_first = para_spectrum[:int(n / 2)]
    para_second = para_spectrum[int(n / 2 + 1):]
    para_flipped = np.flip(para_first)

    perp_first = perp_spectrum[:int(n / 2)]  # skip nyquist
    perp_second = perp_spectrum[int(n / 2 + 1):]
    perp_flipped = np.flip(perp_first)

    perp_total = (perp_second + perp_flipped) / 2.0
    para_total = (para_second + para_flipped) / 2.0

    start = 10
    end = 50

    logk = np.log(range(int((n - 1) / 2)))

    slope, intercept, rval, p, err = linregress(logk[start:end], np.log(perp_total[start:end]))
    slope_par, intercept_par, rval_pa, p_para, err_para = linregress(logk[start:end],
                                                                     np.log(para_total[start:end]))

    lin_perb = [slope * i + intercept for i in np.log(range(start, end))]
    lin_par = [slope_par * i + intercept_par for i in np.log(range(start, end))]

    plt.figure(figsize=(6.5,6.0), dpi=200)
    plt.scatter(logk[start:end], np.log(perp_total[start:end]), label='$log(E(K_{\perp}))$', color='orange')
    plt.scatter(logk[start:end], np.log(para_total[start:end]), label='$log(E(K_{\parallel}))$', color='green')
    plt.plot(logk[start:end], lin_perb, color='orange', label='Slope $K_{\perp}$: %s err: %s' % (round(slope, 3), round(err,3)))
    plt.plot(logk[start:end], lin_par, color='green', label='Slope $K_{\parallel}$: %s err: %s' % (round(slope_par, 3), round(err_para,3)))
    plt.xlabel('Log K')
    plt.ylabel('Log E(k)')
    plt.legend(loc='lower left', fontsize=14)
    plt.title('Phi Power Spectrum 2D %s' %n)
    #plt.show()

    print(slope, slope_par, "np.flip method")

    # Plotting PS_KT files from power_kk method in spectrum.f08

    # perpendicular
    # filename = dir2 + 'PS_KT_PHI0.DAT'
    # ps_kperp = np.loadtxt(filename)
    # log_ps_kperp = [np.log(i) for i in ps_kperp[:, 1] if i.any() != 0]
    #
    # # parallel
    # filename = dir2 + 'PS_KB_PHI0.DAT'
    # ps_kpar = np.loadtxt(filename)
    # log_ps_kpar = [np.log(i) for i in ps_kpar[:, 1] if i.any() != 0]
    #
    # # print(log_ps_kpar[1:])
    # # print(log_ps_kperp[1:])
    #
    # logk = np.log(range(int((n - 1) / 2)))
    #
    # # print(logk)
    #
    # slope, intercept, rval, p, err = linregress(logk[start:end], log_ps_kperp[start:end])
    # slope_par, intercept_par, rval_pa, p_para, err_para = linregress(logk[start:end], log_ps_kpar[start:end])
    #
    # lin_perb = [slope * i + intercept for i in logk[start:end]]
    # lin_par = [slope_par * i + intercept_par for i in logk[start:end]]
    #
    # print(slope, slope_par, "power spectra dir")
    #
    # plt.figure()
    # plt.scatter(logk[start:end], log_ps_kperp[start:end], label='Kperp', color='blue')
    # plt.scatter(logk[start:end], log_ps_kpar[start:end], label='Kpara', color='red')
    # plt.plot(logk[start:end], lin_perb, label='Slope K_perp: %f' % round(slope, 3))
    # plt.plot(logk[start:end], lin_par, label='Slope K_para: %f' % round(slope_par, 3))
    #
    # plt.xlabel('Log K')
    # plt.ylabel('Log E(k)')
    # plt.title('Phi0 Spectrum 2D from spectrum.f08')
    # plt.legend()
    # #plt.show()

def plot_power3d(dirphi0, dir1, n):
    # 3D
    nx = n + 1 #displacement only
    ny = n + 1
    nz = n + 1

    # Reading in PHI0 binary file
    filename = dirphi0 + 'PHI0' + '.BIN'  # PHI for displacement, RHO for squares
    # print(filename)
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    tmp = np.reshape(abx, (nx, ny, nz))
    phi0 = tmp.transpose()

    # Reading in PHI0 binary file
    filename = dir1 + 'PHI' + '.BIN'
    # print(filename)
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    tmp2 = np.reshape(abx, (nx, ny, nz))
    phi = tmp2.transpose()

    # print(phi0[0, :], len(phi0[0, :]))

    # FFT
    phi0k_noshift = fft.fftn(phi0)
    phik_noshift = fft.fftn(phi)

    # Perform Shift
    phi0k = fft.fftshift(phi0k_noshift)
    phik = fft.fftshift(phik_noshift)

    # plotting phik at particular k_perp to see exponential drop off at centre
    # var = abs(phi0k) ** 2
    # plt.figure()
    # plt.plot(np.log(var[:, 2, 2]), label=':2,2')
    # plt.plot(np.log(var[:, 40, 40]), label=':40,40')
    # plt.legend()
    # plt.title('Phi0k at one k_perp to see exp. drop from centre')
    #plt.show()

    # fig = plt.figure()
    # fig = plt.figure(figsize=(5.0, 5.0))
    # gs = gridspec.GridSpec(2, 1, hspace=0.5, wspace=0.2)

    # ax0 = plt.subplot(gs[0], aspect='equal')

    # plt.imshow(np.log(abs(phi0k) ** 2)[:,2,:], cmap='seismic',
    #            interpolation='nearest', origin='lower')

    # fig = plt.gcf()  # get the current figure
    # plt.clim()  # clamp the color limits
    # plt.colorbar()
    # plt.ylabel('K_y')
    # plt.xlabel('K_x')
    # plt.title('FFT of Phi0 3D')

    # ax1 = plt.subplot(gs[1], aspect='equal')

    # plt.imshow(np.log(abs(phik) ** 2)[:,2,:], cmap='seismic',
    #            interpolation='nearest', origin='lower')

    # fig = plt.gcf()  # get current figure
    # plt.clim()  # clamp the color limits
    # plt.colorbar()
    # plt.ylabel('K_z')
    # plt.xlabel('K_x')
    # plt.title('FFT of Phi 3D')
    #plt.show()

    para_spectrum = np.zeros(nx)
    para_total = np.zeros(int(n/2))
    for i in range(nx):
        para_spectrum[i] = np.sum(abs(phik[i, :, :]) ** 2)
    # for w in range(1,n/2): 
    #     para_total[w] = 0.5*(para_spectrum[nx/2-1+w]+para_spectrum[nx/2-1-w])

    para_first = para_spectrum[:int(n / 2)]
    para_second = para_spectrum[int(n / 2 + 1):]
    para_flipped = np.flip(para_first)

    para_total = (para_second + para_flipped) / 2.0
   
    start = 10
    end = 40

    logk = np.log(range(int((n - 1) / 2)))

    perp_total = k_perp_calculator(n,phi,phi0,dir1)

    slope, intercept, rval, p, err = linregress(logk[start:end], np.log(perp_total[start:end]))
    slope_par, intercept_par, rval_pa, p_para, err_para = linregress(logk[start:end],
                                                                     np.log(para_total[start:end]))
    print(slope, slope_par,"3d kperp calc")

    lin_perb = [slope * i + intercept for i in logk[start:end]]
    lin_par = [slope_par * i + intercept_par for i in logk[start:end]]
    #for i in np.log(range(1, int(n / 2)))]

    plt.figure(figsize=(6.5,6.0), dpi=200)
    plt.scatter(logk[start:end], 1.4 * np.log(perp_total[start:end]), label='$log(E(K_{\perp}))$')
    plt.scatter(logk[start:end], np.log(para_total[start:end]), label='$log(E(K_{\parallel}))$')
    plt.plot(logk[start:end], 1.4*np.array(lin_perb), label='Slope $K_{\perp}$: %s err: %s' % (round(slope, 3), round(err,3)))
    plt.plot(logk[start:end], lin_par, label='Slope $K_{\parallel}$: %s err: %s' % (round(slope_par, 3), round(err_para,3)))
    plt.legend(loc='lower left',fontsize=14)
    plt.xlabel('Log K')
    plt.ylabel('Log E(k)')
    plt.title('Phi Power Spectrum 3D B M_A = 1.97 %s' % n)
    plt.show()


if __name__ == "__main__":

    # n = 128
    # dir_data1 = "c:/Users/jonas/DESY/2d_displacement/Runs/128run2D_73_test/"  # this one allows kx == 0
    # dir_data2 = "c:/Users/jonas/DESY/2d_displacement/Runs/128run2D_73_test/power_spectra/"
    # plot_power(dir_data1, dir_data2, n)

    #dir_data3 = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/Runs/128run2D_73_frac/
    #dir_data4 = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/Runs/128run2D_73_frac/power_spectra/"
    #plot_power2d(dir_data3, dir_data4, n)

    #n = 512
    #dir_sf = "/lustre/fs23/group/that/jonas/Github_repo/DESY/final_data/2d/512run2D_disp_real/"
    # dir_data2 = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_vs_3d_data/256run2D_FFT/power_spectra/"
    #dir_sf = '/home/jonas/Documents/VSCode/DESY/final_data/2d/512run2D_disp_real/'
    #dir_data2 = 2
    #plot_power2d(dir_sf, dir_data2, n)

    n = 512
    #dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/3d_disp_mem/Runs/512_amp02"
    dir_phi0 =  "/home/jonas/Documents/VSCode/DESY/phi0init/Runs/512_test/"
    dir_data = '/home/jonas/Documents/VSCode/DESY/3d_disp_mem/Runs/512_amp04/'
    plot_power3d(dir_phi0, dir_data, n)

# #dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/Runs/128run2D_73/power_spectra/"
# dir__data = "c:/Users/jonas/DESY/2d_displacement/Runs/128run2D_73_test/power_spectra/"
