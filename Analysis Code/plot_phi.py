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


nx=129 #129 as 128 boxes but data comes with points for each edge
ny=129
dir_data = "../2d_displacement/data/" # into parent directory

filename=dir_data+'PHI'+'.BIN'
print(filename)  
fd = open(filename, 'rb')

abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny)
  
temp = np.reshape(abx,(nx,ny))

fig=plt.figure()
fig = plt.figure(figsize=(5.0, 5.0))
gs = gridspec.GridSpec(1, 1, hspace=0.2, wspace=0.2)
ax0 = plt.subplot(gs[0],aspect='equal')
z=temp
plt.imshow(z, cmap='seismic', extent=[0, 1, 0, 1],
           interpolation='nearest', origin='lower')

fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()
plt.show()
fig.savefig('test.eps',dpi=300,bbox_inches='tight')
  
  



