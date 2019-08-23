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
max_size = 2048.0
####################################################################################

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/128run2D_FFT/sf_par_perp_v_phiF.txt'
lentf= 128.0
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
lpar1 = lpyare/lentf
lperp1 = lperpe/lentf
#
filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/256run2D_FFT/sf_par_perp_v_phiF.txt'
lentf= 256.0
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
lpar2 = lpyare/lentf
lperp2 = lperpe/lentf
#
filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/512run2D_FFT/sf_par_perp_v_phiF.txt'
lentf= 512.0
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
lpar3 = lpyare/lentf
lperp3 = lperpe/lentf

#filename = 'c:/Users/jonas/DESY/2d_displacement/128run2D_73_frac/sf_par_perp_v_phi0F_re.txt'
filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/1024run2D_FFT/sf_par_perp_v_phiF.txt'
lentf=1024.0
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
lpar8 = lpyare/lentf
lperp8 = lperpe/lentf

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/2048run2D_FFT/sf_par_perp_v_phiF.txt'
lentf=2048.0
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
lpar9 = lpyare/lentf
lperp9 = lperpe/lentf

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/square_128run2D/sf_par_perp_v_phiF.txt'
lentf=128.0
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
lpar10 = lpyare/lentf
lperp10 = lperpe/lentf

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/square_256run2D/sf_par_perp_v_phiF.txt'
lentf=256.0
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
lpar11 = lpyare/lentf
lperp11 = lperpe/lentf

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/square_512run2D/sf_par_perp_v_phiF.txt'
lentf=512.0
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
lpar12 = lpyare/lentf
lperp12 = lperpe/lentf

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/square_1024run2D/sf_par_perp_v_phiF.txt'
lentf=1024.0
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
lpar13 = lpyare/lentf
lperp13 = lperpe/lentf

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_sq_vs_disp_data/square_2048run2D/sf_par_perp_v_phiF.txt'
lentf=2048.0
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
lpar14 = lpyare/lentf
lperp14 = lperpe/lentf

for count_128disp, i in enumerate(lperp1):
   if  i <= 0.0001:
     print(count_128disp)
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

perp_arr,para_arr = [], []

def linfit(perp_arr, para_arr, count):
  slope, intercept, rval, p, err = linregress(np.log(perp_arr[:count]), np.log(para_arr[:count]))
  tmp_slop = round(slope,3)
  tmp_r = round(rval,3)
  tmp_err = round(err,3)
  return tmp_slop, tmp_r, tmp_err

slope_512_disp, rval_512_disp, err_512_disp = linfit(lperp3,lpar3, count_512disp)
slope_1024_disp, rval_1024_disp, err_1024_disp = linfit(lperp8,lpar8, count_1024disp)
slope_2048_disp, rval_2048_disp, err_2048_disp = linfit(lperp9,lpar9, count_2048disp)

slope_512_sq, rval_512_sq, err_512_sq = linfit(lperp12,lpar12, count_512sq)
slope_1024_sq, rval_1024_sq, err_1024_sq = linfit(lperp13,lpar13, count_1024sq)
slope_2048_sq, rval_2048_sq, err_2048_sq = linfit(lperp14,lpar14, count_2048sq)


#reference slopes

ref_slope_2_3 = lpar9[100]*(np.power(lperp9[:count_2048disp],(2.0/3.0))/np.power(lperp9[100],(2.0/3.0)))
# ref_slope_1 = lpar4[3]*(np.power(lperp4,(3.0/3.0))/np.power(lpar4[3],(3.0/3.0)))

slope_ref, rval_ref, err_ref = linfit(lperp9, ref_slope_2_3, count_2048disp)
 
fig=plt.figure()
fig = plt.figure(figsize=(16.0, 10.0))
gs = gridspec.GridSpec(1, 1, hspace=0.0, wspace=0.0)

ax0 = plt.subplot(gs[0])

#ax0.plot(lperp1[:count_128disp], lpar1[:count_128disp], lw=3, ls = "-", label="128_2D_disp")
#ax0.plot(lperp2[:count_256disp], lpar2[:count_256disp], lw=3, ls = "-", label="256_2D_disp")
ax0.plot(lperp3[:count_512disp], lpar3[:count_512disp], lw=3, ls = "-", label="512_2D_disp grad: %s R^2: %s  Err: %s" % (slope_512_disp, rval_512_disp, err_512_disp))
ax0.plot(lperp8[:count_1024disp], lpar8[:count_1024disp], lw=3, ls = "-", label="1024_2D_disp grad: %s R^2: %s  Err: %s" % (slope_1024_disp, rval_1024_disp, err_1024_disp))
ax0.plot(lperp9[:count_2048disp], lpar9[:count_2048disp], lw=3, ls = "-", label="2048_2D_disp grad: %s R^2: %s  Err: %s" % (slope_2048_disp, rval_2048_disp, err_2048_disp))

#ax0.plot(lperp10[:count_128sq], lpar10[:count_128sq], lw=3, ls = ":", label="128_2D_sq")
#ax0.plot(lperp11[:count_256sq], lpar11[:count_256sq], lw=3, ls = ":", label="256_2D_sq")
ax0.plot(lperp12[:count_512sq], lpar12[:count_512sq], lw=5, ls = ":", label="512_2D_sq grad: %s R^2: %s  Err: %s" % (slope_512_sq, rval_512_sq, err_512_sq))
ax0.plot(lperp13[:count_1024sq], lpar13[:count_1024sq], lw=5, ls = ":", label="1024_2D_sq grad: %s R^2: %s  Err: %s" % (slope_1024_sq, rval_1024_sq, err_1024_sq))
ax0.plot(lperp14[:count_2048sq], lpar14[:count_2048sq], lw=5, ls = ":", label="2048_2D_sq grad: %s R^2: %s  Err: %s" % (slope_2048_sq, rval_2048_sq, err_2048_sq))

ax0.plot(lperp9[:count_2048disp], ref_slope_2_3, lw=6, color = "black", ls = "-", label="GS95 Slopegrad: %s R^2: %s  Err: %s" % (slope_ref, rval_ref, err_ref))

ax0.set_xscale('log')
ax0.set_yscale('log')
#ax0.set_xlim(xmin=0.002)
#ax0.set_xlim(xmax=0.4)
#ax0.set_ylim(ymax=0.4)
#ax0.set_ylim(ymin=0.002)
ax0.set_xlabel(r'$l_{\perp}/ L $ perpendicular',fontsize=18)
ax0.set_ylabel('$l_{\parallel}/L $ parallel',fontsize=18)
ax0.set_title('Structure Function 2D Squares vs Displacement')
#ax0.text(3,100,'(a)',fontsize=18)
ax0.legend(loc='lower right',ncol=2,fontsize=14)

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

#plt.savefig('test.eps', bbox_inches='tight',transparent=False)
plt.show()

