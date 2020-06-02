import struct
import numpy as np
#import matplotlib.pyplot as plt
import math
import scipy.interpolate as spint
from multiprocessing import Pool
#from pylab import *

RGI = spint.RegularGridInterpolator

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
n_theta_pts = 50
nzslices = 20
dir_data = "../1024_runs/C4/data/decomped_modes/"
mode='F'
nprocs = 24
nprocsfft = 24
tinynum = 1.0E-12

magkperp = np.zeros(lent/2)
mag_power_spec = np.zeros(lent/2)
kin_power_spec = np.zeros(lent/2)
mag_pspec_tavg = np.zeros(lent/2)
kin_pspec_tavg = np.zeros(lent/2)
f_power_spec = np.zeros(lent/2)
pfolded = np.zeros((xpt/2,ypt/2))
logpfolded = np.zeros((xpt/2,ypt/2))
fieldx = np.zeros((xpt,ypt,zpt))
fieldy = np.zeros((xpt,ypt,zpt))
fieldz = np.zeros((xpt,ypt,zpt))


def kspec_funk(ff):

  print(ff)

  #looping over theta and phi angles in k space
  theta_array = np.linspace(0.0,np.pi/2.0,n_theta_pts)
  count = 0
  pwr_f = 0.0
 
  for theta in theta_array :
    kxx = magkperp[ff]*np.cos(theta)
    kyy = magkperp[ff]*np.sin(theta) 
    pont=(kxx,kyy)
    pwr_f = pwr_f + np.exp(logp_interp(pont))
    count=count+1
  
  #taking shell average
  pwr_f = pwr_f/count
 
  #multiplying by shell area to get 1D power 
  f_pow = pwr_f*2.0*np.pi*magkperp[ff]

  return [f_pow]


def slice_fft(sn):
  fx2d = fieldx[:,:,sn]
  fy2d = fieldy[:,:,sn]
  fz2d = fieldz[:,:,sn]

  #calculating the FFT
  fxk = np.fft.rfftn(fx2d)
  fyk = np.fft.rfftn(fy2d)
  fzk = np.fft.rfftn(fz2d)

  #shifting the zero frequency to center
  fxk_shifted = np.fft.fftshift(fxk,axes=0) 
  fyk_shifted = np.fft.fftshift(fyk,axes=0) 
  fzk_shifted = np.fft.fftshift(fzk,axes=0) 

  #computing the power
  pfxk = np.real(fxk_shifted*np.conjugate(fxk_shifted))
  pfyk = np.real(fyk_shifted*np.conjugate(fyk_shifted))
  pfzk = np.real(fzk_shifted*np.conjugate(fzk_shifted))

  for ii in range (0,np.size(pfxk,0)) :
    for jj in range (0,np.size(pfxk,1)) :
      if (math.isinf(pfxk[ii,jj]) or math.isnan(pfxk[ii,jj]) or pfxk[ii,jj]<tinynum):
        pfxk[ii,jj] = tinynum
      if (math.isinf(pfyk[ii,jj]) or math.isnan(pfyk[ii,jj]) or pfyk[ii,jj]<tinynum):
        pfyk[ii,jj] = tinynum
      if (math.isinf(pfzk[ii,jj]) or math.isnan(pfzk[ii,jj]) or pfzk[ii,jj]<tinynum):
        pfzk[ii,jj] = tinynum

  pfxyzk = pfxk+pfyk+pfzk        

  #removing the Nyquist component
  pfxyzk_wn = pfxyzk[1:xpt,0:ypt/2] 

  return(pfxyzk_wn)


ntstp = 0
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
  fieldx = temp.transpose()
  temp = np.reshape(aby,(lent,lent,lent))
  fieldy = temp.transpose()
  temp = np.reshape(abz,(lent,lent,lent))
  fieldz = temp.transpose()
 
  #fieldx = bx-np.mean(bx)
  #bx = None
  #fieldy = by-np.mean(by)
  #by = None
  #fieldz = bz-np.mean(bz)
  #bz = None

  #removing the mean fields
  #fieldx = bx - np.mean(bx)
  #bx = None
  #fieldy = by - np.mean(fieldy)
  #fieldz = fieldz - np.mean(fieldz)

  #fieldx = bx
  #fieldy = by
  #fieldz = bz
  #temp = None
  #abx = None
  #aby = None
  #abz = None
 
  #taking slices over the parallel direction
  pk_slice = np.zeros((xpt-1,ypt/2)) 
  """ 
  count = 0
  for sn in range (0,lent,lent/nzslices) :
    
    fx2d = fieldx[:,:,sn]
    fy2d = fieldy[:,:,sn]
    fz2d = fieldz[:,:,sn]

    #calculating the FFT
    fxk = np.fft.rfftn(fx2d)
    fyk = np.fft.rfftn(fy2d)
    fzk = np.fft.rfftn(fz2d)

    #shifting the zero frequency to center
    fxk_shifted = np.fft.fftshift(fxk,axes=0) 
    fyk_shifted = np.fft.fftshift(fyk,axes=0) 
    fzk_shifted = np.fft.fftshift(fzk,axes=0) 
 
    #computing the power
    pfxk = np.real(fxk_shifted*np.conjugate(fxk_shifted))
    pfyk = np.real(fyk_shifted*np.conjugate(fyk_shifted))
    pfzk = np.real(fzk_shifted*np.conjugate(fzk_shifted))
    pfxyzk = pfxk+pfyk+pfzk        

    #removing the Nyquist component
    pfxyzk_wn = pfxyzk[1:xpt,0:ypt/2] 

    pk_slice = pk_slice+pfxyzk_wn
    count = count+1

  pk_slice = pk_slice/count
  """
  
  if __name__ == '__main__':  
    pool = Pool(processes=nprocsfft) 
    oter = pool.map(slice_fft, np.arange(0,lent,lent/nzslices))
    fftslices = np.asarray(oter)
    pool.terminate()

  pk_slice = np.mean(fftslices,axis=0)
  

  #folding along the y axis
  pfolded[0,:] = pk_slice[xpt/2-1,:]
  for w in range (1,xpt/2):
    pfolded[w,:] = 0.5*(pk_slice[xpt/2-1+w,:]+pk_slice[xpt/2-1-w,:]) 
    
  #taking the log for easier interpolation
  logpfolded = np.log(pfolded)

  #defining the axes of the 3d spec array
  kxarr = np.linspace(0.0,(xpt/2-1)*2.0*np.pi,num=xpt/2)
  kyarr = np.linspace(0.0,(ypt/2-1)*2.0*np.pi,num=ypt/2)
  
  #making the interpolation function
  logp_interp = RGI(points=[kxarr,kyarr], values=logpfolded)

  #making the spectrum now
  #zero wavenumber is special case
  magkperp[0] = 0.0
  pont = (0.0,0.0)
  f_power_spec[0] = np.exp(logp_interp(pont))
  
  #looping over the higher wavenumbers
  for i in range (1,lent/2) :
    magkperp[i] = i*2.0*np.pi
  
  if __name__ == '__main__':  
    pool = Pool(processes=nprocs) 
    kspecs = pool.map(kspec_funk, np.arange(1,lent/2))
    sff = np.asarray(kspecs)
    pool.terminate()
 
  f_power_spec[1:] = sff[:,0]
  
  #adding up the spectra at different times
  mag_pspec_tavg = mag_pspec_tavg+f_power_spec  

  #bx = None
  #by = None
  #bz = None 

  #now processing the velocity field
  filename=dir_data+'V'+mode+str(t)+'.BIN'
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
  fieldx = temp.transpose()
  temp = np.reshape(aby,(lent,lent,lent))
  fieldy = temp.transpose()
  temp = np.reshape(abz,(lent,lent,lent))
  fieldz = temp.transpose()

  """
  temp = np.reshape(abx,(lent,lent,lent))
  abx = None
  bx = temp.transpose() 
  temp = None
  fieldx = bx-np.mean(bx)
  bx = None

  temp = np.reshape(aby,(lent,lent,lent))
  aby = None
  by = temp.transpose() 
  temp = None
  fieldy = by-np.mean(by)
  by = None

  temp = np.reshape(abz,(lent,lent,lent))
  abz = None
  bz = temp.transpose() 
  temp = None
  fieldz = bz-np.mean(bz)
  bz = None
  """

  #print('sample velocity fiel=',fieldx[514,83,145],fieldy[6,415,415],fieldz[76,546,936])

  """  
  temp = np.reshape(avx,(lent,lent,lent))
  fieldx = temp.transpose()
  temp = np.reshape(avy,(lent,lent,lent))
  fieldy = temp.transpose()
  temp = np.reshape(avz,(lent,lent,lent))
  fieldz = temp.transpose()

  #removing the mean fields
  fieldx = fieldx - np.mean(fieldx)
  fieldy = fieldy - np.mean(fieldy)
  fieldz = fieldz - np.mean(fieldz)
  
  #fieldx = vx
  #fieldy = vy
  #fieldz = vz
  temp = None
  avx = None
  avy = None
  avz = None 
  """

  #taking slices over the parallel direction
  pk_slice = np.zeros((xpt-1,ypt/2))
  """
  count = 0
  for sn in range (0,lent,lent/nzslices) :
    
    fx2d = fieldx[:,:,sn]
    fy2d = fieldy[:,:,sn]
    fz2d = fieldz[:,:,sn]

    #calculating the FFT
    fxk = np.fft.rfftn(fx2d)
    fyk = np.fft.rfftn(fy2d)
    fzk = np.fft.rfftn(fz2d)

    #shifting the zero frequency to center
    fxk_shifted = np.fft.fftshift(fxk,axes=0) 
    fyk_shifted = np.fft.fftshift(fyk,axes=0) 
    fzk_shifted = np.fft.fftshift(fzk,axes=0) 
 
    #computing the power
    pfxk = np.real(fxk_shifted*np.conjugate(fxk_shifted))
    pfyk = np.real(fyk_shifted*np.conjugate(fyk_shifted))
    pfzk = np.real(fzk_shifted*np.conjugate(fzk_shifted))
    pfxyzk = pfxk+pfyk+pfzk        

    #removing the Nyquist component
    pfxyzk_wn = pfxyzk[1:xpt,0:ypt/2] 

    pk_slice = pk_slice+pfxyzk_wn
    count = count+1

  pk_slice = pk_slice/count
  """
  
  if __name__ == '__main__':  
    pool = Pool(processes=nprocsfft) 
    oter = pool.map(slice_fft, np.arange(0,lent,lent/nzslices))
    fftslices = np.asarray(oter)
    pool.terminate()

  pk_slice = np.mean(fftslices,axis=0)
  

  #folding along the y axis
  pfolded[0,:] = pk_slice[xpt/2-1,:]
  for w in range (1,xpt/2):
    pfolded[w,:] = 0.5*(pk_slice[xpt/2-1+w,:]+pk_slice[xpt/2-1-w,:]) 
    
  #taking the log for easier interpolation
  logpfolded = np.log(pfolded)
 
  #making the interpolation function
  logp_interp = RGI(points=[kxarr,kyarr], values=logpfolded)

  #making the spectrum now
  #zero wavenumber is special case
  pont = (0.0,0.0)
  f_power_spec[0] = np.exp(logp_interp(pont))
  
  if __name__ == '__main__':  
    pool = Pool(processes=nprocs) 
    kspecs = pool.map(kspec_funk, np.arange(1,lent/2))
    sff = np.asarray(kspecs)
    pool.terminate()
 
  f_power_spec[1:] = sff[:,0]

  #adding up the spectra at different times
  kin_pspec_tavg = kin_pspec_tavg+f_power_spec  

  #vx = None
  #vy = None
  #vz = None 

  ntstp = ntstp+1


#taking the time average
mag_pspec_tavg = mag_pspec_tavg/ntstp
kin_pspec_tavg = kin_pspec_tavg/ntstp

#writing the spectra to a file
f=open(dir_data+mode+'_kperp_spec.txt','w')
for i in range(0,lent/2):
  value=str(magkperp[i])+" "+str(mag_pspec_tavg[i])+" "+str(kin_pspec_tavg[i])
  f.write(value+"\n")
f.close()




