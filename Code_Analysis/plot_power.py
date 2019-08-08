import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

n = 128

dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/128run2D_73_paratest/power_spectra/"  # data files
dir_output = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/128run2D_73_paratest/power_spectra/"

# perpendicular
filename = dir_data + 'PS_KT_PHI.DAT'
ps_kperp = np.loadtxt(filename)
log_ps_kperp = [np.log(i) for i in ps_kperp[:,1] if i.any() != 0]

# parallel
filename = dir_data + 'PS_KB_PHI.DAT'
ps_kpar = np.loadtxt(filename)
log_ps_kpar = [np.log(i) for i in ps_kpar[:,1] if i.any() != 0]

#print(len(log_ps_kperp))S
#print(len(log_ps_kpar))

#print(log_ps_kpar)
#print(log_ps_kperp)

logk = np.log(range(n/2))

#print(len(logk))

slope, intercept, rval, p, err = linregress(logk, log_ps_kperp[1:])
slope_par, intercept_par, rval_pa, p_para, err_para = linregress(logk, log_ps_kpar[1:])

print(slope, slope_par)

plt.figure()
plt.plot(logk, log_ps_kperp[1:], label='Kperp')
plt.plot(logk, log_ps_kpar[1:], label ='Kpara')
plt.legend()
plt.show()