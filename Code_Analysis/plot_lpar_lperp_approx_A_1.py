import matplotlib.pyplot as plt
import numpy as np
#import read,amrplot
from matplotlib import ticker
from matplotlib import gridspec
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.stats import linregress

def smoothing(xarr) :
  lls = xarr.size
  smoothx = np.zeros(lls)
  smoothx[0] = xarr[0]
  smoothx[lls-1] = xarr[lls-1]
  for i in range (1,lls-1) :
    smoothx[i] = 0.25*xarr[i-1]+0.5*xarr[i]+0.25*xarr[i+1]
  return smoothx

def lppcorr(llv,sfpar,sfperp) :
  lls = sfpar.size
  lperp_arr = np.zeros(lls)
  lpar_arr = np.zeros(lls)
  count = 0
  for i in range (1,lls) :
    lpar = i*1.0
    stfc = sfpar[i]
    if (stfc > np.amax(sfperp) or stfc < np.amin(sfperp)) : 
      continue
    #trying to find the ll for which sf_par matches the value of stfc
    for j in range (0,lls-1) :
      if ((sfperp[j] < stfc) and (sfperp[j+1] >= stfc)) : # this is where the 2d_squares struc funk fails because the sfperp[j+1] >= sfpar[i] for all (all decrease)
        xll = llv[j]+((llv[j+1]-llv[j])/(sfperp[j+1]-sfperp[j]))*(stfc-sfperp[j])
        break
    lperp_arr[count] = xll
    lpar_arr[count]  = lpar
    count = count+1  
  return [lperp_arr,lpar_arr]

####################################################################################
#max_size = 2048.0
####################################################################################

def read_sf(dir_data, n):
  filename = dir_data
  lentf= n
  data = np.loadtxt(filename,skiprows=1)
  ll = data[:,0]
  sf_par = data[:,1]
  sf_perp= data[:,2]
  valid = ~np.isnan(sf_perp)
  sf_perp = sf_perp[valid]
  ll = ll[valid]
  sf_par = sf_par[valid]
  lent = np.size(ll)
  sf_par_smoothed = smoothing(sf_par)
  sf_perp_smoothed= smoothing(sf_perp)
  [lperpe,lpyare] = lppcorr(ll,sf_par_smoothed,sf_perp_smoothed)
  lpar = lpyare/lentf
  lperp = lperpe/lentf

  return lpar, lperp

#2d displacement sf phi
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/128run2D_FFT/sf_par_perp_v_phiF.txt'
lpar1, lperp1 = read_sf(dir_sf, 128.0)
#
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/256run2D_FFT/sf_par_perp_v_phiF.txt'
lpar2, lperp2 = read_sf(dir_sf, 256.0)
#
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/512run2D_FFT/sf_par_perp_v_phiF.txt'
lpar3, lperp3 = read_sf(dir_sf, 512.0)
#
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/1024run2D_FFT/sf_par_perp_v_phiF.txt'
lpar8, lperp8 = read_sf(dir_sf, 1024.0)
#
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/2048run2D_FFT/sf_par_perp_v_phiF.txt'
lpar9, lperp9 = read_sf(dir_sf, 2048.0)
#
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/4096run2D_FFT/sf_par_perp_v_phiF.txt'
lpar16, lperp16 = read_sf(dir_sf, 4096.0)

#2d_squares sf rho(=phi) 
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/square_128run2D/sf_par_perp_v_phiF.txt'
lpar10, lperp10 = read_sf(dir_sf, 128.0)
#
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/square_256run2D/sf_par_perp_v_phiF.txt'
lpar11, lperp11 = read_sf(dir_sf, 256.0)
#
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/square_512run2D/sf_par_perp_v_phiF.txt'
lpar12, lperp12 = read_sf(dir_sf, 512.0)
#
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/square_1024run2D/sf_par_perp_v_phiF.txt'
lpar13, lperp13 = read_sf(dir_sf, 1024.0)
#
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/square_2048run2D/sf_par_perp_v_phiF.txt'
lpar14, lperp14 = read_sf(dir_sf, 2048.0)
#
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/square_4096run2D/sf_par_perp_v_phiF.txt'
lpar15, lperp15 = read_sf(dir_sf, 4096.0)

#3d displacement phi0
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_vs_3d_data/128run2D_FFT/sf_par_perp_v_phi0F.txt'
lpar17, lperp17 = read_sf(dir_sf, 128.0)
#
dir_sf = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_vs_3d_data/256run2D_FFT/sf_par_perp_v_phi0F.txt'
lpar18, lperp18 = read_sf(dir_sf, 256.0)
#
filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_vs_3d_data/512run2D_FFT/sf_par_perp_v_phi0F.txt'
lpar19, lperp19 = read_sf(dir_sf, 512.0)

#perp and para - data becomes 0 at some point, run into errors, so want to find point that they become 0 and stop at that point

#2d displacement phi
for count_128disp, i in enumerate(lperp1):
   if  i <= 0.0001:
     print(count_128disp) #returns the index where 0 starts
     break

for count_256disp, i in enumerate(lperp2):
   if  i <= 0.0001:
     print(count_256disp)
     break

for count_512disp, i in enumerate(lperp3):
   if  i <= 0.0001:
     print(count_512disp)
     break

for count_1024disp, i in enumerate(lperp8):
   if  i <= 0.0001:
     print(count_1024disp)
     break

for count_2048disp, i in enumerate(lperp9):
   if  i <= 0.0001:
     print(count_2048disp)
     break

for count_4096disp, i in enumerate(lperp16):
   if  i <= 0.0001:
     print(count_4096disp)
     break

#2d squares rho
for count_128sq, i in enumerate(lperp10):
   if  i <= 0.0001:
     print(count_128sq)
     break

for count_256sq, i in enumerate(lperp11):
   if  i <= 0.0001:
     print(count_256sq)
     break

for count_512sq, i in enumerate(lperp12):
  if  i <= 0.0001:
    print(count_512sq)
    break

for count_1024sq, i in enumerate(lperp13):
  if  i <= 0.0001:
    print(count_1024sq)
    break  

for count_2048sq, i in enumerate(lperp14):
  if  i <= 0.0001:
    print(count_2048sq)
    break      

for count_4096sq, i in enumerate(lperp15):
  if  i <= 0.0001:
    print(count_4096sq)
    break      

#2d displacement phi0
for count_128disp_phi0, i in enumerate(lperp17):
   if  i <= 0.0001:
     print(count_128disp_phi0) #returns the index where 0 starts
     break

for count_256disp_phi0, i in enumerate(lperp18):
   if  i <= 0.0001:
     print(count_256disp_phi0)
     break

for count_512disp_phi0, i in enumerate(lperp19):
   if  i <= 0.0001:
     print(count_512disp_phi0)
     break

#3d displacement phi0
# for count_128disp_3dphi0, i in enumerate(lperp20):
#    if  i <= 0.0001:
#      print(count_128disp_3dphi0) #returns the index where 0 starts
#      break

# for count_256disp_3dphi0, i in enumerate(lperp21):
#    if  i <= 0.0001:
#      print(count_256disp_3dphi0)
#      break

#3d displacement phi
# for count_128disp_3dphi, i in enumerate(lperp20):
#    if  i <= 0.0001:
#      print(count_128disp_3dphi) #returns the index where 0 starts
#      break

# for count_256disp_3dphi, i in enumerate(lperp21):
#    if  i <= 0.0001:
#      print(count_256disp_3dphi0)
#      break

perp_arr,para_arr = [], []

def linfit(perp_arr, para_arr, count):
  slope, intercept, rval, p, err = linregress(np.log(perp_arr[:count]), np.log(para_arr[:count]))
  tmp_slop = round(slope,3)
  tmp_r = round(rval,3)
  tmp_err = round(err,3)
  return tmp_slop, tmp_r, tmp_err

#slopes linefitting for 2d displacement phi
slope_128_disp, rval_128_disp, err_128_disp = linfit(lperp1,lpar1, count_128disp)
slope_256_disp, rval_256_disp, err_256_disp = linfit(lperp2,lpar2, count_256disp)
slope_512_disp, rval_512_disp, err_512_disp = linfit(lperp3,lpar3, count_512disp)
slope_1024_disp, rval_1024_disp, err_1024_disp = linfit(lperp8,lpar8, count_1024disp)
slope_2048_disp, rval_2048_disp, err_2048_disp = linfit(lperp9,lpar9, count_2048disp)
slope_4096_disp, rval_4096_disp, err_4096_disp = linfit(lperp16,lpar16, count_4096disp)

#slope linefitting for 2d squares rho (=phi)
slope_512_sq, rval_512_sq, err_512_sq = linfit(lperp12,lpar12, count_512sq)
slope_1024_sq, rval_1024_sq, err_1024_sq = linfit(lperp13,lpar13, count_1024sq)
slope_2048_sq, rval_2048_sq, err_2048_sq = linfit(lperp14,lpar14, count_2048sq)
slope_4096_sq, rval_4096_sq, err_4096_sq = linfit(lperp15,lpar15, count_4096sq)


#Reference slopes
ref_slope_2_3 = lpar16[100]*(np.power(lperp16[:count_4096disp],(2.0/3.0))/np.power(lperp16[100],(2.0/3.0)))
#ref_slope_1 = lpar4[3]*(np.power(lperp4,(3.0/3.0))/np.power(lpar4[3],(3.0/3.0)))
slope_ref, rval_ref, err_ref = linfit(lperp16, ref_slope_2_3, count_4096disp)


#plot for 2d squares vs displacement phi 
fig=plt.figure(1)
fig = plt.figure(figsize=(16.0, 10.0))
gs = gridspec.GridSpec(1, 1, hspace=0.0, wspace=0.0)

ax0 = plt.subplot(gs[0])

#2D displacement phi
ax0.plot(lperp3[:count_512disp], lpar3[:count_512disp], lw=3, ls = "-", label="512_2D_disp grad: %s R^2: %s  Err: %s" % (slope_512_disp, rval_512_disp, err_512_disp))
ax0.plot(lperp8[:count_1024disp], lpar8[:count_1024disp], lw=3, ls = "-", label="1024_2D_disp grad: %s R^2: %s  Err: %s" % (slope_1024_disp, rval_1024_disp, err_1024_disp))
ax0.plot(lperp9[:count_2048disp], lpar9[:count_2048disp], lw=3, ls = "-", label="2048_2D_disp grad: %s R^2: %s  Err: %s" % (slope_2048_disp, rval_2048_disp, err_2048_disp))
ax0.plot(lperp16[:count_4096disp], lpar16[:count_4096disp], lw=3, ls = "-", label="4096_2D_disp grad: %s R^2: %s  Err: %s" % (slope_4096_disp, rval_4096_disp, err_4096_disp))

#2D squares rho(=phi)
ax0.plot(lperp12[:count_512sq], lpar12[:count_512sq], lw=5, ls = ":", label="512_2D_sq grad: %s R^2: %s  Err: %s" % (slope_512_sq, rval_512_sq, err_512_sq))
ax0.plot(lperp13[:count_1024sq], lpar13[:count_1024sq], lw=5, ls = ":", label="1024_2D_sq grad: %s R^2: %s  Err: %s" % (slope_1024_sq, rval_1024_sq, err_1024_sq))
ax0.plot(lperp14[:count_2048sq], lpar14[:count_2048sq], lw=5, ls = ":", label="2048_2D_sq grad: %s R^2: %s  Err: %s" % (slope_2048_sq, rval_2048_sq, err_2048_sq))
ax0.plot(lperp14[:count_4096sq], lpar14[:count_4096sq], lw=5, ls = ":", label="4096_2D_sq grad: %s R^2: %s  Err: %s" % (slope_4096_sq, rval_4096_sq, err_4096_sq))

ax0.plot(lperp16[:count_4096disp], ref_slope_2_3, lw=6, color = "black", ls = "-", label="GS95 grad: %s R^2: %s  Err: %s" % (slope_ref, rval_ref, err_ref))

ax0.set_xscale('log')
ax0.set_yscale('log')
ax0.set_xlabel(r'$l_{\perp}/ L $ perpendicular',fontsize=18)
ax0.set_ylabel('$l_{\parallel}/L $ parallel',fontsize=18)
ax0.set_title('Structure Function 2D Squares vs Displacement')
ax0.legend(loc='lower right',ncol=2,fontsize=14)

plt.show()

#plot for 2d vs 3d displacement method both phi and phi0
fig=plt.figure(2)
fig = plt.figure(figsize=(16.0, 10.0))
gs = gridspec.GridSpec(2, 1, hspace=0.0, wspace=0.0)

ax0 = plt.subplot(gs[0],aspect='equal')

#2D displacement PHI0
ax0.plot(lperp17[:count_128disp_phi0], lpar17[:count_128disp_phi0], lw=3, ls = "-", label="128_2D_disp_PHI0 grad: %s R^2: %s  Err: %s" % (slope_128_disp_phi0, rval_128_disp_phi0, err_128_disp_phi0))
ax0.plot(lperp18[:count_256disp_phi0], lpar18[:count_256disp_phi0], lw=3, ls = "-", label="256_2D_disp_PHI0 grad: %s R^2: %s  Err: %s" % (slope_256_disp_phi0, rval_256_disp_phi0, err_256_disp_phi0))
ax0.plot(lperp19[:count_512disp_phi0], lpar19[:count_512disp_phi0], lw=3, ls = "-", label="512_2D_disp_PHI0 grad: %s R^2: %s  Err: %s" % (slope_512_disp_phi0, rval_512_disp_phi0, err_512_disp_phi0))

ax0.set_xlabel(r'$l_{\perp}/ L $ perpendicular',fontsize=18)
ax0.set_ylabel('$l_{\parallel}/L $ parallel',fontsize=18)
ax0.set_title('Struc Funk 2D vs 3D Displacement PHI0')
ax0.legend(loc='lower right',ncol=2,fontsize=14)




ax1 = plt.subplot(gs[1],aspect='equal')

#2D displacement PHI
ax1.plot(lperp3[:count_128disp], lpar3[:count_128disp], lw=3, ls = "-", label="128_2D_disp_PHI grad: %s R^2: %s  Err: %s" % (slope_128_disp, rval_128_disp, err_128_disp))
ax1.plot(lperp8[:count_256disp], lpar8[:count_256disp], lw=3, ls = "-", label="256_2D_disp_PHI grad: %s R^2: %s  Err: %s" % (slope_256_disp, rval_256_disp, err_256_disp))
ax1.plot(lperp3[:count_512disp], lpar3[:count_512disp], lw=3, ls = "-", label="512_2D_disp_PHI grad: %s R^2: %s  Err: %s" % (slope_512_disp, rval_512_disp, err_512_disp))

ax1.set_xlabel(r'$l_{\perp}/ L $ perpendicular',fontsize=18)
ax1.set_ylabel('$l_{\parallel}/L $ parallel',fontsize=18)
ax1.set_title('Struc Funk 2D vs 3D Displacement PHI')
ax1.legend(loc='lower right',ncol=2,fontsize=14)




#plt.show()

"""
ax1 = plt.subplot(gs[1])
ax1.plot(lperp_v_1, lpar_v_1, lw=3,label="S1")
ax1.plot(lperp_v_2, lpar_v_2, lw=3,label="S2")
ax1.plot(lperp_v_3, lpar_v_3, lw=3,label="S3")
ax1.plot(lperp_v_4, lpar_v_4, lw=3,label="S4")
ax1.plot(lperp_v_5, lpar_v_5, lw=3,label="C1")
ax1.plot(lperp_v_6, lpar_v_6, lw=3,label="C2")
ax1.plot(lperp_v_7, lpar_v_7, lw=3,label="C4")
ax1.plot(lperp_B_4, ref_slope_1, lw=4,ls=":",label="Isotropic",color="blue")
ax1.plot(lperp_B_4, ref_slope_2_3, lw=4,ls=":",label="GS",color="red")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(xmin=2.0)
ax1.set_xlim(xmax=120)
ax1.set_ylim(ymin=2.0)
ax1.set_ylim(ymax=130.0)
ax1.set_xlabel(r'$l_{\perp}$',fontsize=18)
ax1.set_ylabel('$l_{\parallel}$',fontsize=18)
ax1.text(3,100,'(b)',fontsize=18)
ax1.legend(loc='lower right',ncol=2,fontsize=14)
"""




