"""
TODO:
    1.  check which dimension is kperp and kpara
    2.  integrate up kperp and kpara - get spectrum in each (just by adding up each for constant other value)
    3.  plot them against each other
    4.  do this for multiple dimensions/methods
    5.  save the plots in new folder
FIXME:
    
"""

import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import gridspec

# data files
dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/128run2D_73_test/"  

#2D
n = 128
nx = n
ny = n

#Reading in PHI0 binary file
filename=dir_data+'PHI0'+'.BIN' # PHI for displacement, RHO for squares
print(filename)  
fd = open(filename, 'rb')

abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny)
  
tmp = np.reshape(abx,(nx,ny))
phi0 = tmp.transpose()

#Reading in PHI0 binary file
filename=dir_data+'PHI'+'.BIN'
print(filename)  
fd = open(filename, 'rb')

abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny)
  
tmp2 = np.reshape(abx,(nx,ny))
phi = tmp2.transpose()

print(phi0[0,:], len(phi0[0,:]))

#FFT
phi0k_noshift = fft.fft2(phi0)
phik_noshift = fft.fft2(phi)

#Perform Shift
phi0k = fft.fftshift(phi0k_noshift)
phik = fft.fftshift(phik_noshift)

#print(phi0k)
#print('~')
#print(phik)

fig=plt.figure()
fig = plt.figure(figsize=(5.0, 5.0))
gs = gridspec.GridSpec(2, 1, hspace=0.2, wspace=0.2)

ax0 = plt.subplot(gs[0],aspect='equal')

plt.imshow(np.log(abs(phi0k)**2), cmap='seismic',
           interpolation='nearest', origin='lower')
fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()

ax1 = plt.subplot(gs[1],aspect='equal')

plt.imshow(np.log(abs(phik)**2), cmap='seismic',
           interpolation='nearest', origin='lower')
fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()

plt.show()
#fig.savefig('test.eps',dpi=300,bbox_inches='tight')


