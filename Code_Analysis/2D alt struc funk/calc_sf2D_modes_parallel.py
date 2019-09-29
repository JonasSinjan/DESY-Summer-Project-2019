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

lent=1024
xpt=lent
ypt=lent
zpt=lent
Lx=1.0
Ly=1.0
Lz=1.0
t_start = 2
t_stop = 6
step = 1
dir_data = "../1024_runs/C4/data/decomped_modes/"
dir_output = "../1024_runs/C4/data/decomped_modes/"
seed(1)
n_avg_bfield_pts = 5
nrandpts = 2000
nprocs = 24
mode= 'F'

ntstp = 0
sf2D_array=np.zeros((lent/4,lent/4))

def struct2D_funk(ipar):

  lpar = ipar*1.0
  print("ipar = ", ipar)

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
      lr = rand(n_avg_bfield_pts)*sp_rad
      theta = rand(n_avg_bfield_pts)*np.pi
      phi   = rand(n_avg_bfield_pts)*2.0*np.pi
      lx = lr*np.sin(theta)*np.cos(phi)
      ly = lr*np.sin(theta)*np.sin(phi)
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

      b1 = np.array([  bx[r1[0],r1[1],r1[2]]  , by[r1[0],r1[1],r1[2]] , bz[r1[0],r1[1],r1[2]]   ])
      b2 = np.array([  bx[r2[0],r2[1],r2[2]]  , by[r2[0],r2[1],r2[2]] , bz[r2[0],r2[1],r2[2]]   ])
      
      diff = b1-b2
      summb = np.sum(diff*diff)
      if (math.isinf(summb) or math.isnan(summb)):
        print("isinf isnan")
        continue
  
      ssff = ssff + summb
      numpt = numpt + 1
   
    sf_lperp[iperp] = ssff/numpt

  return sf_lperp

 
  """

    #bhat contains the unit vector along the local magnetic field
    #now take 2 points separated by distance ll along this direction with center at ri
    r1 = np.int_(ri+(ll/2.0)*bhat)
    r2 = np.int_(ri-(ll/2.0)*bhat)

    #renormalize indices
    r1 = r1%lent
    r2 = r2%lent

    b1 = np.array([  bx[r1[0],r1[1],r1[2]]  , by[r1[0],r1[1],r1[2]] , bz[r1[0],r1[1],r1[2]]   ])
    b2 = np.array([  bx[r2[0],r2[1],r2[2]]  , by[r2[0],r2[1],r2[2]] , bz[r2[0],r2[1],r2[2]]   ])

    sf_pare = sf_pare + np.sum((b1-b2)*(b1-b2))

    #next calculating the perpendicular structure function
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
    r1 = np.int_(ri+(ll/2.0)*perpb)
    r2 = np.int_(ri-(ll/2.0)*perpb)

    #renormalize indices
    r1 = r1%lent
    r2 = r2%lent

    b1 = np.array([  bx[r1[0],r1[1],r1[2]]  , by[r1[0],r1[1],r1[2]] , bz[r1[0],r1[1],r1[2]]   ])
    b2 = np.array([  bx[r2[0],r2[1],r2[2]]  , by[r2[0],r2[1],r2[2]] , bz[r2[0],r2[1],r2[2]]   ])

    sf_perpe = sf_perpe + np.sum((b1-b2)*(b1-b2))

    numpt = numpt + 1.0

  #print(ll,numpt,sf_pare,sf_perpe)
  return [numpt,sf_pare,sf_perpe]

  """

for t in range (t_start,t_stop+1,step) : #the time loop 
  
  filename=dir_data+'B'+mode+str(t)+'.BIN'
  print(filename)
  fd = open(filename, 'rb')
  fd.read(4)
  nx = np.fromfile(file=fd,dtype=np.int32,count=1)[0]
  ny = np.fromfile(file=fd,dtype=np.int32,count=1)[0]
  nz = np.fromfile(file=fd,dtype=np.int32,count=1)[0]
  print(nx,ny,nz)
  fd.read(4)
  fd.read(4)
  abx = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)
  fd.read(4)
  fd.read(4)
  aby = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)
  fd.read(4)
  fd.read(4)
  abz = np.fromfile(file=fd,dtype=np.float64,count=nx*ny*nz)
  fd.read(4)

  temp = np.reshape(abx,(lent,lent,lent))
  bx = temp.transpose()
  temp = np.reshape(aby,(lent,lent,lent))
  by = temp.transpose()
  temp = np.reshape(abz,(lent,lent,lent))
  bz = temp.transpose()

  #print(np.amax(bx),np.amax(by),np.amax(bz))

  if __name__ == '__main__':  
    pool = Pool(processes=nprocs) 
    sf2D_list = pool.map(struct2D_funk, range(lent/4))
    sf2D_arr = np.asarray(sf2D_list)
    pool.terminate() 

  sf2D_array = sf2D_array + sf2D_arr
  ntstp = ntstp + 1

sf2D_array = sf2D_array/ntstp

np.save(dir_output+'sf2D_'+mode+'.npy',sf2D_array)

"""
sf_par = sf_par/npts
sf_perp=sf_perp/npts

#writing the spectra to a file
f=open(dir_output+'sf_par_perp_F.txt','w')
for i in range(0,lent/2):
  value=str(i*1.0)+" "+str(sf_par[i])+" "+str(sf_perp[i])
  f.write(value+"\n")
f.close()
"""










