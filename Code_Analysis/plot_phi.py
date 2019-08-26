import struct
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import math
import scipy.interpolate as spint
from numpy.random import seed
from numpy.random import randint
from numpy.random import rand
from matplotlib import ticker
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D

res = 512
nx=res
ny=res
nz=res



def read_phi_2d(dir_data, n):
        res=n
        nx=res
        ny=res
        nz=res
        filename=dir_data+'PHI'+'.BIN'
        print(filename)  
        fd = open(filename, 'rb')

        abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny)
        
        temp1 = np.reshape(abx,(nx,ny))

        filename=dir_data+'PHI0'+'.BIN'
        print(filename)
        fd = open(filename, 'rb')

        abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny)

        temp2 = np.reshape(abx,(nx,ny))

        return temp1, temp2

def read_phi_3d(dir_data, n):
        res=n
        nx=res
        ny=res
        nz=res
        filename=dir_data+'PHI'+'.BIN'
        print(filename)  
        fd = open(filename, 'rb')

        abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)
        
        temp1 = np.reshape(abx,(nx,ny,nz))

        filename=dir_data+'PHI0'+'.BIN'
        print(filename)
        fd = open(filename, 'rb')

        abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)

        temp2 = np.reshape(abx,(nx,ny,nz))

        filename=dir_data+'BY'+'.BIN'
        print(filename)
        fd = open(filename, 'rb')

        abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)

        temp3 = np.reshape(abx,(nx,ny,nz))

        return temp1, temp2, temp3

dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_vs_3d_data/256run2D_FFT/"
n=256
phi_2d,phi0_2d = read_phi_2d(dir_data,n)

dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/3d_displacement/128_test/"
n=128
phi_3d,phi0_3d, by = read_phi_3d(dir_data,n)

# #2d
# fig=plt.figure(1)
# fig = plt.figure(figsize=(5.0, 5.0))
# gs = gridspec.GridSpec(2, 1, hspace=0.2, wspace=0.2)
# ax0 = plt.subplot(gs[0],aspect='equal')

# z = phi0_2d


# plt.imshow(z, cmap='seismic', extent=[0, 1, 0, 1],
#         interpolation='nearest', origin='lower')
# fig = plt.gcf()
# plt.clim()   # clamp the color limits
# plt.colorbar()

# ax1 = plt.subplot(gs[1],aspect='equal')

# z=phi_2d

# plt.imshow(z, cmap='seismic', extent=[0, 1, 0, 1],
#         interpolation='nearest', origin='lower')

# fig = plt.gcf()
# plt.clim()   # clamp the color limits
# plt.colorbar()


#3d
#slice_index = randint(0,res)

slice_index = 23

fig=plt.figure(2)
fig = plt.figure(figsize=(5.0, 5.0))
gs = gridspec.GridSpec(2, 1, hspace=0.2, wspace=0.2)
ax0 = plt.subplot(gs[0],aspect='equal')

z=phi0_3d[:,:,slice_index] #one slice


plt.imshow(z, cmap='seismic', extent=[0, 1, 0, 1],
        interpolation='nearest', origin='lower')
fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()

ax1 = plt.subplot(gs[1],aspect='equal')

z=phi_3d[:,:,slice_index] #one slice

plt.imshow(z, cmap='seismic', extent=[0, 1, 0, 1],
        interpolation='nearest', origin='lower')

fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()

slice_index = 64

fig=plt.figure(3)
fig = plt.figure(figsize=(5.0, 5.0))
gs = gridspec.GridSpec(2, 1, hspace=0.2, wspace=0.2)
ax0 = plt.subplot(gs[0],aspect='equal')

z=by[:,slice_index,:] #one slice

plt.imshow(z, cmap='seismic', extent=[0, 1, 0, 1],
        interpolation='nearest', origin='lower')
fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()

ax1 = plt.subplot(gs[1],aspect='equal')

z=phi_3d[:,slice_index,:] #one slice

plt.imshow(z, cmap='seismic', extent=[0, 1, 0, 1],
        interpolation='nearest', origin='lower')

fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()

plt.show()
