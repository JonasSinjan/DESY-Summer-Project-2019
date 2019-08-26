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
	temp1 = temp1.transpose()

        filename=dir_data+'PHI0'+'.BIN'
        print(filename)
        fd = open(filename, 'rb')

        abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)

        temp2 = np.reshape(abx,(nx,ny,nz))
	temp2 = temp2.transpose()

        filename=dir_data+'BY'+'.BIN'
        print(filename)
        fd = open(filename, 'rb')

        abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)

        temp3 = np.reshape(abx,(nx,ny,nz))

        return temp1, temp2, temp3

dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_vs_3d_data/128run2D_FFT/"
n=129
phi_2d,phi0_2d = read_phi_2d(dir_data,n)

dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_vs_3d_data/128run3D_FFT/"
n=129
phi_3d,phi0_3d, by = read_phi_3d(dir_data,n)


slice_index = 23
plt.figure(constrained_layout=True)
fig, axs = plt.subplots(nrows=2, ncols=2)

z=[phi0_2d, phi_2d, phi0_3d[:,:,slice_index], phi_3d[:,:,slice_index]]#one slice

def phi_plot(ax,z, title):
    ax.imshow(z, cmap='seismic', extent=[0, 1, 0, 1],
        interpolation='nearest', origin='lower')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title('%s' % (title), fontsize=14)
    fig = plt.gcf()
    plt.clim()   # clamp the color limits
    plt.colorbar()

title=['2D Phi0', '2D Phi', '3D Phi0', '3D Phi']

for ax, count in enumerate(axs.flat):
        phi_plot(ax,z[count],title[count])


plt.show()

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
