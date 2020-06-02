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


dir_data = "/lustre/fs23/group/that/jonas/Github_repo/DESY/3d_displacement/128_test/"
n=129
nx = n
ny = n
nz = n
filename=dir_data+'BY'+'.BIN'
print(filename)
fd = open(filename, 'rb')
abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)
#temp1 = abx.astype('float64')
temp3 = np.reshape(abx,(nx,ny,nz))

print(temp3[12,:,34])
