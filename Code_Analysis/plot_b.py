import numpy as np
import matplotlib.pyplot as plt

dir_data = r"lustre/fs23/group/that/jonas/Github_repo/DESY/localB/Runs/512_B_amp05/"

n = 513
nx, ny, nz = n, n, n

filename=dir_data+'BX'+'.BIN'
print(filename)
fd = open(filename, 'rb')

abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)

bx = np.reshape(abx,(nx,ny,nz))

filename=dir_data+'BY'+'.BIN'
print(filename)
fd = open(filename, 'rb')

abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)

by = np.reshape(abx,(nx,ny,nz))

filename=dir_data+'BZ'+'.BIN'
print(filename)
fd = open(filename, 'rb')

abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)

bz = np.reshape(abx,(nx,ny,nz))

b = bx + by + bz

plt.figure()
plt.imshow(b[:,:,23])
plt.show()