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

# data files
# dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/128run2D_73_test/"
dir_data = "c:/Users/jonas/DESY/2d_displacement/128run2D_73_test/"

# 2D
n = 128
nx = n + 1
ny = n + 1

# Reading in PHI0 binary file
filename = dir_data + 'PHI0' + '.BIN'  # PHI for displacement, RHO for squares
print(filename)
fd = open(filename, 'rb')

abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

tmp = np.reshape(abx, (nx, ny))
phi0 = tmp.transpose()

# Reading in PHI0 binary file
filename = dir_data + 'PHI' + '.BIN'
print(filename)
fd = open(filename, 'rb')

abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny)

tmp2 = np.reshape(abx, (nx, ny))
phi = tmp2.transpose()

print(phi0[0, :], len(phi0[0, :]))

# FFT
phi0k_noshift = fft.fft2(phi0)
phik_noshift = fft.fft2(phi)

# Perform Shift
phi0k = fft.fftshift(phi0k_noshift)
phik = fft.fftshift(phik_noshift)

# print(phi0k)
# print('~')
# print(phik)
tmp_arr = np.zeros((20, 20))
for i in range(20):
    tmp_arr[i, :] = i

# print(tmp_arr)


fig = plt.figure(1)
fig = plt.figure(figsize=(5.0, 5.0))
gs = gridspec.GridSpec(2, 1, hspace=0.5, wspace=0.2)

ax0 = plt.subplot(gs[0], aspect='equal')

plt.imshow(np.log(abs(phi0k) ** 2), cmap='seismic',
           interpolation='nearest', origin='lower')
# plt.imshow(tmp_arr)
fig = plt.gcf()  # get the current figure
plt.clim()  # clamp the color limits
plt.colorbar()
plt.ylabel('K_parallel')
plt.xlabel('K_perp')
plt.title('FFT of Phi0')

ax1 = plt.subplot(gs[1], aspect='equal')

plt.imshow(np.log(abs(phik) ** 2), cmap='seismic',
           interpolation='nearest', origin='lower')
fig = plt.gcf()  # get current figure
plt.clim()  # clamp the color limits
plt.colorbar()
plt.ylabel('K_parallel')
plt.xlabel('K_perp')
plt.title('FFT of Phi')
plt.savefig('test')
plt.show()

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

start = 5
end = 40

slope, intercept, rval, p, err = linregress(np.log(range(start, end)), np.log(perp_total[start:end]))
slope_par, intercept_par, rval_pa, p_para, err_para = linregress(np.log(range(start, end)),
                                                                 np.log(para_total[start:end]))

lin_perb = [slope * i + intercept for i in np.log(range(1, int(n / 2)))]
lin_par = [slope_par * i + intercept_par for i in np.log(range(1, int(n / 2)))]

plt.figure(2)
plt.plot(np.log(range(start, end)), np.log(perp_total[start:end]), label='Perp')
plt.plot(np.log(range(start, end)), np.log(para_total[start:end]), label='Para')
plt.legend()
# plt.show()

print(slope, slope_par)

# Plotting PS_KT files from power_kk method in spectrum.f08

n = 128
#
# #dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/Runs/128run2D_73/power_spectra/"  # data files
dir_data = "c:/Users/jonas/DESY/2d_displacement/128run2D_73_test/power_spectra/"

# perpendicular
filename = dir_data + 'PS_KT_PHI.DAT'
ps_kperp = np.loadtxt(filename)
log_ps_kperp = [np.log(i) for i in ps_kperp[:, 1] if i.any() != 0]

# parallel
filename = dir_data + 'PS_KB_PHI.DAT'
ps_kpar = np.loadtxt(filename)
log_ps_kpar = [np.log(i) for i in ps_kpar[:, 1] if i.any() != 0]

# print(log_ps_kpar[1:])
# print(log_ps_kperp[1:])

logk = np.log(range(int((n - 1) / 2)))

# print(logk)

slope, intercept, rval, p, err = linregress(logk[start:end], log_ps_kperp[start:end])
slope_par, intercept_par, rval_pa, p_para, err_para = linregress(logk[start:end], log_ps_kpar[start:end])

lin_perb = [slope * i + intercept for i in logk[start:end]]
lin_par = [slope_par * i + intercept_par for i in logk[start:end]]

print(slope, slope_par)

plt.figure(2)
plt.scatter(logk[start:end], log_ps_kperp[start:end], label='Kperp', color='blue')
plt.scatter(logk[start:end], log_ps_kpar[start:end], label='Kpara', color='red')
plt.plot(logk[start:end], lin_perb, label='Slope K_perp: %f' % round(slope, 3))
plt.plot(logk[start:end], lin_par, label='Slope K_para: %f' % round(slope_par, 3))

plt.xlabel('Log K')
plt.ylabel('Log E(k)')

plt.legend()
plt.show()
