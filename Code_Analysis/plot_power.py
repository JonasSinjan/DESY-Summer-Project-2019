import matplotlib.pyplot as plt
import numpy as np

n = 128

dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/128run2D_73_paratest/power_spectra/"  # data files
dir_output = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/128run2D_73_paratest/power_spectra/"

# perpendicular
filename = dir_data + 'PS_KT_PHI.DAT'
ps_kperp = np.loadtxt(filename)
log_ps_kperp = np.log(ps_kperp)

# parallel
filename = dir_data + 'PS_KB_PHI.DAT'
ps_kpar = np.loadtxt(filename)
log_ps_kpar = np.log(ps_kpar)

start = -((n-1)/2 + 1)
logk = [0]*len(range(n-2))

for i in range(n-2):
    tmp = -start + i
    if tmp == 0:
        continue
    logk[i] = np.log(tmp)

plt.figure()
plt.plot(logk, log_ps_kperp)
plt.plot(logk, log_ps_kpar)
plt.show()