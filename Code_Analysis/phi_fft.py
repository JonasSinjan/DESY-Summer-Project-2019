import numpy as np
import numpy.fft as fft

# data files
dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/Runs/128run2D_73/"  
dir_output = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/Runs/128run2D_73/"

#2D
n = 128
nx = n
ny = n

#Reading in PHI0 binary file
filename=dir_data+'PHI0'+'.BIN' # PHI for displacement, RHO for squares
print(filename)  
fd = open(filename, 'rb')

abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny)
  
phi0 = np.reshape(abx,(nx,ny))

#Reading in PHI0 binary file
filename=dir_data+'PHI'+'.BIN'
print(filename)  
fd = open(filename, 'rb')

abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny)
  
phi = np.reshape(abx,(nx,ny))

#FFT
phi0k_noshift = fft.fft2(phi0)
phik_noshift = fft.fft2(phi)

#Perform Shift
phi0k = fft.fftshift(phi0k_noshift)
phik = fft.fftshift(phik_noshift)

print(phi0k)
print('~')
print(phik)



