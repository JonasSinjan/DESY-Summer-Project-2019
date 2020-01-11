import struct
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.interpolate as spint
from numpy.random import seed
from numpy.random import randint
from numpy.random import rand
from matplotlib import ticker
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D

def read_phi_2d(dir_data, n):
        res=n
        nx=res
        ny=res
        filename=dir_data+'PHI'+'.BIN'
        print(filename)  
        fd = open(filename, 'rb')

        abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny)
        
        temp1 = np.reshape(abx,(nx,ny))
        #temp1 = temp1.transpose()

        filename=dir_data+'PHI0'+'.BIN'
        print(filename)
        fd = open(filename, 'rb')

        abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny)

        temp2 = np.reshape(abx,(nx,ny))
        #temp2 = temp2.transpose()

        return temp1, temp2

def read_phi_3d(dir_data, n):
        nx=n
        ny=n
        nz=n
        filename=dir_data+'PHI'+'.BIN'
        print(filename)  
        fd = open(filename, 'rb')

        abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)
        
        temp1 = np.reshape(abx,(nx,ny,nz))
        #temp1 = temp1.transpose()

        # filename=dir_data+'PHI0'+'.BIN'
        # print(filename)
        # fd = open(filename, 'rb')

        # abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)

        # temp2 = np.reshape(abx,(nx,ny,nz))
        #temp2 = temp2.transpose()

        # filename=dir_data+'BY'+'.BIN'
        # print(filename)
        # fd = open(filename, 'rb')

        # abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)

        # temp3 = np.reshape(abx,(nx,ny,nz))

        return temp1#, temp2#, temp3

working_dir = r"/home/jonas/Documents/VSCode/DESY/"#r"/lustre/fs23/group/that/jonas/Github_repo/DESY/"

n = 513
nx=n
ny=n
nz=n
filename= working_dir + r"phi0init/Runs/512_test/PHI0.BIN"
print(filename)  
fd = open(filename, 'rb')

abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)
phi0_512 = np.reshape(abx,(nx,ny,nz))


#for 512 amp 05 phi.bin file
<<<<<<< HEAD
dir_data = r"/lustre/fs23/group/that/jonas/Github_repo/DESY/3d_disp_mem/Runs/512_fix_r/"
=======
dir_data = working_dir + r"3d_disp_mem/Runs/512_amp05/"
>>>>>>> 575ba075eb0f8c3083c8297109c5d1a13c449f47
#dir_phi0 = "/home/jonas/Documents/VSCode/DESY/3d_displacement/Runs/128_FFT_mem_test_compar/"
n=513 #(n+1)_
phi_amp05 = read_phi_3d(dir_data,n)

dir_data = working_dir + r"3d_disp_mem/Runs/512_amp1/"
#dir_data = "/home/jonas/Documents/VSCode/DESY/3d_disp_mem/Runs/128_mem_method_test_FFT/"
#(n+1)_
phi_amp1 = read_phi_3d(dir_data,n)

dir_data = working_dir + r"final_data/2d/512run2D_disp_real/"

phi_2d, phi0_2d = read_phi_2d(dir_data, 513)
slice_index = 23

fig = plt.figure(figsize=(8.0, 8.0))

gs0 = gridspec.GridSpec(1,2,hspace=0.2, wspace=0.2)

gs1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0])


#z=[phi0_2d, phi_2d, phi0_3d[:,:,slice_index], phi_3d[:,:,slice_index]]#one slice
title=['2D Phi0', '2D Phi', '3D Phi0', '3D Phi']
#title=['3D MEM Phi0', '3D MEM Phi', '3D Phi0', '3D Phi']
#title = ['3D Phi0', '3D PHI', 'PHI-PHI0', ' 3D PHI']
z=[phi0_2d, phi_2d, phi0_512[:,:,slice_index], phi_amp05[:,:,slice_index]]#one slice

ax0 = fig.add_subplot(gs1[0])

plt.imshow(z[0], cmap='seismic', extent=[0, 1, 0, 1], interpolation='nearest', origin='lower')
ax0.set_xlabel('x', fontsize=12)
ax0.set_ylabel('y', fontsize=12)
ax0.set_title(title[0], fontsize=14)
fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()

ax1 = fig.add_subplot(gs1[1])

plt.imshow(z[1], cmap='seismic', extent=[0, 1, 0, 1], interpolation='nearest', origin='lower')
ax1.set_xlabel('x', fontsize=12)
ax1.set_ylabel('y', fontsize=12)
ax1.set_title(title[1], fontsize=14)
fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()

gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1])

ax2 = fig.add_subplot(gs2[0])

plt.imshow(z[2], cmap='seismic', extent=[0, 1, 0, 1], interpolation='nearest', origin='lower')
ax2.set_xlabel('x', fontsize=12)
ax2.set_ylabel('y', fontsize=12)
ax2.set_title(title[2], fontsize=14)
fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()

ax3 = fig.add_subplot(gs2[1])

plt.imshow(z[3], cmap='seismic', extent=[0, 1, 0, 1], interpolation='nearest', origin='lower')
ax3.set_xlabel('x', fontsize=12)
ax3.set_ylabel('y', fontsize=12)
ax3.set_title(title[3], fontsize=14)
fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()

plt.show()

# fig = plt.figure(figsize=(8.0, 8.0))

# gs0 = gridspec.GridSpec(2,1,hspace=0.2, wspace=0.2)

# ax0 = fig.add_subplot(gs0[0])

# img = z[1]-z[0]

# plt.imshow(img, cmap='seismic', extent=[0, 1, 0, 1], interpolation='nearest', origin='lower')
# ax0.set_xlabel('x', fontsize=12)
# ax0.set_ylabel('y', fontsize=12)
# ax0.set_title('Amp05 Phi - Phi0', fontsize=14)
# plt.gcf()
# plt.clim()   # clamp the color limits
# plt.colorbar()

# ax1 = fig.add_subplot(gs0[1])

# img = z[2]-z[0]

# plt.imshow(img, cmap='seismic', extent=[0, 1, 0, 1], interpolation='nearest', origin='lower')
# ax1.set_xlabel('x', fontsize=12)
# ax1.set_ylabel('y', fontsize=12)
# ax1.set_title('Amp1 Phi - Phi0', fontsize=14)
# plt.gcf()
# plt.clim()   # clamp the color limits
# plt.colorbar()
# plt.show()

# for i in range(3):
#         slice_index = randint(0,60)
#         #z=[phi0_2d, phi_2d, phi0_3d[:,:,slice_index], phi_3d[:,:,slice_index]]#one slice
        

#         fig = plt.figure(figsize=(12.0, 10.0))

#         gs0 = gridspec.GridSpec(2,1,hspace=0.2, wspace=0.2)

#         ax0 = fig.add_subplot(gs0[0])

#         img = z[1]-z[0]

#         plt.imshow(img, cmap='seismic', extent=[0, 1, 0, 1], interpolation='nearest', origin='lower')
#         ax0.set_xlabel('x', fontsize=12)
#         ax0.set_ylabel('y', fontsize=12)
#         ax0.set_title('2D Phi Change', fontsize=14)
#         plt.gcf()
#         plt.clim()   # clamp the color limits
#         plt.colorbar()

#         ax1 = fig.add_subplot(gs0[1])

#         img = z[3]-z[2]

#         plt.imshow(img, cmap='seismic', extent=[0, 1, 0, 1], interpolation='nearest', origin='lower')
#         ax1.set_xlabel('x', fontsize=12)
#         ax1.set_ylabel('y', fontsize=12)
#         ax1.set_title('3D Phi Change , slice: %s' % (slice_index), fontsize=14)
#         plt.gcf()
#         plt.clim()   # clamp the color limits
#         plt.colorbar()

# plt.show()

