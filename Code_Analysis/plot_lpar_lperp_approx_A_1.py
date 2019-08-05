import matplotlib.pyplot as plt
import numpy as np
#import read,amrplot
from matplotlib import ticker
from matplotlib import gridspec
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline

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
max_size = 512.0
####################################################################################

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/128run2D/sf_par_perp_v_F.txt'
lentf = 128.0
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

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/256run2D/sf_par_perp_v_F.txt'
lentf = 256.0
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

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/512run2D/sf_par_perp_v_F.txt'
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
lpar3 = lpyare/lentf
lperp3 = lperpe/lentf

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_squares/128run_sq/sf_par_perp_v_F.txt'
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
lpar4 = lpyare/lentf
lperp4 = lperpe/lentf

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_squares/256run_sq/sf_par_perp_v_F.txt'
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
lpar8 = lpyare/lentf
lperp8 = lperpe/lentf

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_squares/512run_sq/sf_par_perp_v_F.txt'
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
lpar9 = lpyare/lentf
lperp9 = lperpe/lentf

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/128run2D_73/sf_par_perp_v_F.txt'
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
lpar5 = lpyare/lentf
lperp5 = lperpe/lentf

filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/256run2D_73/sf_par_perp_v_F.txt'
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
lpar6 = lpyare/lentf
lperp6 = lperpe/lentf


filename = '/lustre/fs23/group/that/jonas/Github_repo/DESY/2d_displacement/512run2D_73/sf_par_perp_v_F.txt'
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
lpar7 = lpyare/lentf
lperp7 = lperpe/lentf

# in lper and lpar arrays they stop and become 0 - unsure why - this is a filter to slice the array for plotting due to log scale errors with zero otherwise
for count, i in enumerate(lperp5):
  if  i <= 0.0001:
    print(count)
    break

for count_256, i in enumerate(lperp6):
  if  i <= 0.0001:
    print(count_256)
    break

for count_512, i in enumerate(lperp7):
  if  i <= 0.0001:
    print(count_512)
    break

for count_128sq, i in enumerate(lperp4):
   if  i <= 0.0001:
     print(count_128sq)
     break

for count_256sq, i in enumerate(lperp8):
   if  i <= 0.0001:
     print(count_256sq)
     break

for count_512sq, i in enumerate(lperp9):
   if  i <= 0.0001:
     print(count_512sq)
     break

#reference slopes

ref_slope_2_3 = lpar7[10]*(np.power(lperp7[:count_512],(2.0/3.0))/np.power(lperp7[12],(2.0/3.0)))
# ref_slope_1 = lpar4[3]*(np.power(lperp4,(3.0/3.0))/np.power(lpar4[3],(3.0/3.0)))
 
fig=plt.figure()
fig = plt.figure(figsize=(10.0, 6.0))
gs = gridspec.GridSpec(1, 1, hspace=0.0, wspace=0.0)

ax0 = plt.subplot(gs[0])
#ax0.plot(lperp1[:23], lpar1[:23], lw=3,label="128")
#ax0.plot(lperp2[:42], lpar2[:42], lw=3,label="256")
#ax0.plot(lperp3[:108], lpar3[:108], lw=3,label="512")
ax0.plot(lperp4[:count_128sq], lpar4[:count_128sq], lw=3, ls = "-", label="128_SQ")
ax0.plot(lperp8[:count_256sq], lpar8[:count_256sq], lw=3, ls = "-", label="256_SQ")
ax0.plot(lperp9[:count_512sq], lpar9[:count_512sq], lw=3, ls = "-", label="512_SQ")
ax0.plot(lperp5[:count], lpar5[:count], ls = "--", lw=4,label="128_73")
ax0.plot(lperp6[:count_256], lpar6[:count_256], ls = "--", lw=4, label="256_73")
ax0.plot(lperp7[:count_512], lpar7[:count_512], lw=4, ls = "--", label="512_73")
# ax0.plot(lperp4, ref_slope_1, lw=5,ls=":",label="Isotropic",color="blue")
ax0.plot(lperp7[:count_512], ref_slope_2_3, lw=4, label="GS95", ls = "-", color="black")
ax0.set_xscale('log')
ax0.set_yscale('log')
ax0.set_xlim(xmin=0.002)
ax0.set_xlim(xmax=0.3)
ax0.set_ylim(ymax=0.3)
ax0.set_ylim(ymin=0.002)
ax0.set_xlabel(r'$l_{\perp}/ L $ perpendicular',fontsize=18)
ax0.set_ylabel('$l_{\parallel}/L $ parallel',fontsize=18)
ax0.set_title('Structure Function')
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

