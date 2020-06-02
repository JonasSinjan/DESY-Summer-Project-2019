import struct
import numpy as np
#import matplotlib.pyplot as plt
from multiprocessing import Pool
import math
import scipy.interpolate as spint
from numpy.random import seed
from numpy.random import randint
from numpy.random import rand
#from pylab import *



def read_files3D(dir_data, global_var):

  filename = dir_data + 'PHI0' + '.BIN'
  fd = open(filename, 'rb')

  abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

  temp = np.reshape(abx, (nx, ny, nz))
  phi0 = temp.transpose()

  if global_var == False:
    """
    filename = dir_data + 'PHI' + '.BIN'
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(abx, (nx, ny, nz))
    phi = temp.transpose()
    """
    filename = dir_data + 'BX' + '.BIN' 
    fd = open(filename, 'rb')

    abx = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(abx, (nx, ny, nz))
    bx = temp.transpose()

    filename = dir_data + 'BY' + '.BIN'
    fd = open(filename, 'rb')

    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(aby, (nx, ny, nz))
    by = temp.transpose()

    filename = dir_data + 'BZ' + '.BIN'
    fd = open(filename, 'rb')

    aby = np.fromfile(file=fd, dtype=np.float64, count=nx * ny * nz)

    temp = np.reshape(aby, (nx, ny, nz))
    bz = temp.transpose()
  
  elif global_var:
    bx = np.ones((nx, ny, nz))
    by = np.zeros((nx, ny, nz))
    bz = np.zeros((nx, ny, nz))

  return phi0, bx, by, bz

def struct2D_funk(ipar, phi, bx, by, bz):

  lpar = ipar*1.0
  print('ipar = ', ipar)

  sf_lperp = np.zeros(lent/4)

  #looping over lperp for each lpar
  for iperp in range (0,lent/4) :

    lperp = iperp*1.0 

    ssff = 0.0
    numpt = 0 
    for kup in range (0,nrandpts) : 

      #choose a random point
      ri = randint(0,lent,size=3)
 
      #radius of sphere over which average b field should be taken
      sp_rad = np.sqrt(lpar*lpar+lperp*lperp)/2.0

      #calculate the average b field direction in a sphere of radius sp_rad around random point (xi,yi,zi)
      #do this by looping over npts_avg_field
      #print(rand(n_avg_bfield_pts))
      lr = rand(n_avg_bfield_pts)*sp_rad
      theta = rand(n_avg_bfield_pts)*np.pi
      phi_ang   = rand(n_avg_bfield_pts)*2.0*np.pi
      lx = lr*np.sin(theta)*np.cos(phi_ang)
      ly = lr*np.sin(theta)*np.sin(phi_ang)
      lz = lr*np.cos(theta)
      xis = np.int_(np.floor(ri[0]+lx))
      yis = np.int_(np.floor(ri[1]+ly))
      zis = np.int_(np.floor(ri[2]+lz))

      #renormalize indexes if they are going out of box
      xis = xis%lent
      yis = yis%lent
      zis = zis%lent

      bhat = np.array([ np.mean(bx[xis,yis,zis]) , np.mean(by[xis,yis,zis])  , np.mean(bz[xis,yis,zis]) ])
      bhat = bhat/(np.sqrt(np.sum(bhat*bhat)))

      #make a unit vector perpendicular to bhat
      #generate another random vector
      vrand = rand(3)
      #taking cross product with bhat to generate perpendicular vector
      perpb_x = bhat[1]*vrand[2]-bhat[2]*vrand[1]
      perpb_y = bhat[2]*vrand[0]-bhat[0]*vrand[2]
      perpb_z = bhat[0]*vrand[1]-bhat[1]*vrand[0]
      perpb = np.array([ perpb_x, perpb_y, perpb_z ])
      perpb = perpb/(np.sqrt(np.sum(perpb*perpb)))
      
      #now take 2 points separated by distance ll along the perpendicular direction with center at ri
      r1 = np.int_(ri+(lperp/2.0)*perpb+(lpar/2.0)*bhat)
      r2 = np.int_(ri-(lperp/2.0)*perpb-(lpar/2.0)*bhat)

      #renormalize indices
      r1 = r1%lent
      r2 = r2%lent

      #phi here?
      # b1 = np.array([  bx[r1[0],r1[1],r1[2]]  , by[r1[0],r1[1],r1[2]] , bz[r1[0],r1[1],r1[2]]   ])
      # b2 = np.array([  bx[r2[0],r2[1],r2[2]]  , by[r2[0],r2[1],r2[2]] , bz[r2[0],r2[1],r2[2]]   ])
      b1 = np.array([phi[r1[0], r1[1], r1[2]]])
      b2 = np.array([phi[r2[0], r2[1], r2[2]]])

      diff = b1-b2
      summb = np.sum(diff*diff)
      if (math.isinf(summb) or math.isnan(summb)):
        print('isinf isnan')
        continue
  
      ssff = ssff + summb
      numpt = numpt + 1
   
    sf_lperp[iperp] = ssff/numpt

  return sf_lperp

  #print(np.amax(bx),np.amax(by),np.amax(bz))

if __name__ == '__main__': 

  lent=512
  nx=lent
  ny=lent
  nz=lent
  Lx=1.0
  Ly=1.0
  Lz=1.0

  seed(1)
  n_avg_bfield_pts = 5
  nrandpts = 1000

  mode = 'F'

  ntstp = 0
  sf2D_array=np.zeros((lent/4,lent/4))

  #working_dir_path = '/home/jonas/Documents/VSCode/DESY/'
  working_dir_path = '/lustre/fs23/group/that/jonas/Github_repo/DESY/'
  
  dir_data = working_dir_path + 'phi0init/Runs/512_test/'#final_data/3d/256run3D_FFT/'#'3d_disp_mem/Runs/256_2nd_B/'
  dir_output = working_dir_path + 'phi0init/Runs/512_test/'#'final_data/3d/256run3D_FFT/'#'3d_disp_mem/Runs/256_2nd_B/'

  sf2D_list=[0]*int(lent / 4)
  phi0,bx,by,bz = read_files3D(dir_data, True) 

  input_var = phi0
  for i in range(int(lent / 4)): #varying lpar
    sf2D_list[i] = struct2D_funk(i, input_var, bx, by, bz) #output sf_lperp
  
  sf2D_arr = np.asarray(sf2D_list)
  # print(np.shape(sf2D_arr))
  # print(sf2D_arr[40:1])


sf2D_array = sf2D_array + sf2D_arr
ntstp = ntstp + 1

#sf2D_array = sf2D_array/ntstp

np.save(dir_output+'sf2D_' + 'phi0_test_' + str(nrandpts) + '.npy',sf2D_array)

"""
sf_par = sf_par/npts
sf_perp=sf_perp/npts

#writing the spectra to a file
f=open(dir_output+'sf_par_perp_F_phi.txt','w')
for i in range(0,lent/2):
  value=str(i*1.0)+" "+str(sf_par[i])+" "+str(sf_perp[i])
  f.write(value+"\n")
f.close()
"""










