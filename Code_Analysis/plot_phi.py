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

res = 64
nx=res
ny=res
nz=res

dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/3d_displacement/64run3D/"

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

#print(temp2[:,22])

for i in range(5):
    slice_index = rand(res)

    fig=plt.figure(i)
    fig = plt.figure(figsize=(5.0, 5.0))
    gs = gridspec.GridSpec(2, 1, hspace=0.2, wspace=0.2)
    ax0 = plt.subplot(gs[0],aspect='equal',projection='3d')

    z_arr=temp1[slice_index,:,:] #one slice

    plt.imshow(z, cmap='seismic', extent=[0, 1, 0, 1],
                interpolation='nearest', origin='lower')
    fig = plt.gcf()
    plt.clim()   # clamp the color limits
    plt.colorbar()

    ax1 = plt.subplot(gs[1],aspect='equal')

    z=temp2[slice_index,:,:] #one slice

    plt.imshow(z, cmap='seismic', extent=[0, 1, 0, 1],
            interpolation='nearest', origin='lower')

    fig = plt.gcf()
    plt.clim()   # clamp the color limits
    plt.colorbar()
    plt.title('%f'% slice_index)

plt.show()
