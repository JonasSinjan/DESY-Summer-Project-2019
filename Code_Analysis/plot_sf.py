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


def read_sf(dir_data, n):
  filename = dir_data
  lentf = n
  data = np.loadtxt(filename,skiprows=1)
  ll = data[:,0]
  sf_par = data[:,1]
  sf_perp = data[:,2]
  valid = ~np.isnan(sf_perp)
  sf_perp = sf_perp[valid]
  ll = ll[valid]
  sf_par = sf_par[valid]
  lent = np.size(ll)
  sf_par_smoothed = smoothing(sf_par)
  sf_perp_smoothed = smoothing(sf_perp)
  [lperpe,lpyare] = lppcorr(ll,sf_par_smoothed,sf_perp_smoothed)
  lpar = lpyare/lentf
  lperp = lperpe/lentf

  return lpar, lperp

def find_index(arr): #normally pass through the perp array: lperp
  for count, i in enumerate(arr):
    if  i <= 0.0001:
      return count #returns the index where 0 starts
      break

#perp_arr,para_arr = [], []

def linfit(perp_arr, para_arr, start, end):
  slope, intercept, rval, p, err = linregress(np.log(perp_arr[start:end]), np.log(para_arr[start:end]))
  tmp_slop = round(slope,3)
  tmp_r = round(rval,3)
  tmp_err = round(err,3)
  return tmp_slop, tmp_r, tmp_err 

def process(dir_sf, resolution, start):
  lpar_tmp, lperp_tmp = read_sf(dir_sf, resolution)
  count_tmp = find_index(lperp_tmp)
  slope_tmp, rval_tmp, err_tmp = linfit(lperp_tmp, lpar_tmp, start, count_tmp)
  return [lpar_tmp, lperp_tmp, slope_tmp, err_tmp, rval_tmp, count_tmp]

def plot(process_return, name):
  count = process_return[-1]
  lpar = process_return[0]
  lperp = process_return[1]
  slope = process_return[2]
  err = process_return[3]
  rval = process_return[4]
  plt.plot(lperp[0:count], lpar[0:count], label = "%s grad: %s +/-  %s" % (name, slope, err))

#reading in and processing the structure function files
working_dir_path = '/lustre/fs23/group/that/jonas/Github_repo/DESY/'#r'/home/jonas/Documents/VSCode/DESY/' #

#phi0 wrt global
tmp = working_dir_path + r'phi0init/Runs/512_test/sf_par_perp_v_phi0_wrt_global_10_kpara_2F.txt'
phi0_wrt_global_10kpara2 = process(tmp, 512.0, 4)
#phi0 wrt local amp05
# tmp = working_dir_path + r"localB/Runs/sf_par_perp_v_phi0_wrt_local_amp05F.txt"
# phi0_wrt_local_amp05 = process(tmp, 512.0, 0)
# #phi0 wrt local amp1
# dir_sf = working_dir_path + r'localB/Runs/sf_par_perp_v_phi0_wrt_local_amp1F.txt'
# phi0_wrt_local_amp1 = process(dir_sf, 512.0, 0)
# #phi wrt local amp05
# dir_sf = working_dir_path + r'localB/Runs/512_B_amp05/sf_par_perp_v_phi_wrt_local_fixF.txt'
# phi_wrt_local_amp05 = process(dir_sf, 512.0, 0)
# #phi wrt local amp1
# dir_sf = working_dir_path + r'localB/Runs/512_B_amp1/sf_par_perp_v_phi_wrt_localF.txt'
# phi_wrt_local_amp1 = process(dir_sf, 512.0, 0)

# #phi wrt local sf amp02,3,4
# dir_sf = working_dir_path + r'localB/Runs/sf_par_perp_v_phi_wrt_local_amp02F.txt'
# phi_wrt_local_amp02 = process(dir_sf, 512.0, 0)

# dir_sf = working_dir_path + r'localB/Runs/sf_par_perp_v_phi_wrt_local_amp03F.txt'
# phi_wrt_local_amp03 = process(dir_sf, 512.0, 0)

# dir_sf = working_dir_path + r'localB/Runs/sf_par_perp_v_phi_wrt_local_amp04F.txt'
# phi_wrt_local_amp04 = process(dir_sf, 512.0, 0)

# #phi0 wrt local sf amp02,3,4
# dir_sf = working_dir_path + r'localB/Runs/sf_par_perp_v_phi0_wrt_local_amp02F.txt'
# phi0_wrt_local_amp02 = process(dir_sf, 512.0, 0)

# dir_sf = working_dir_path + r'localB/Runs/sf_par_perp_v_phi0_wrt_local_amp03F.txt'
# phi0_wrt_local_amp03 = process(dir_sf, 512.0, 0)

# dir_sf = working_dir_path + r'localB/Runs/sf_par_perp_v_phi0_wrt_local_amp04F.txt'
# phi0_wrt_local_amp04 = process(dir_sf, 512.0, 0)

# #2D 512 original sf

# #rea;
# dir_sf = working_dir_path + r'final_data/2d/512run2D_disp_real/sf_par_perp_v_phi0_wrt_globalF.txt' #same as old sf file
# phi0_wrt_globalalt = process(dir_sf, 512.0, 0)

# dir_sf = working_dir_path + r'final_data/2d/512run2D_disp_real/sf_par_perp_v_phi0_wrt_localF.txt'
# phi0_wrt_localalt = process(dir_sf, 512.0, 0)

# dir_sf = working_dir_path + r'final_data/2d/512run2D_disp_real/sf_par_perp_v_phiF.txt' #same as old sf file
# phi_wrt_localalt = process(dir_sf, 512.0, 0)

# #FFT
# dir_sf = working_dir_path + r'final_data/2d/512run2D_disp_FFT/sf_par_perp_v_phi0_wrt_globalF.txt' #same as old sf file
# phi0_wrt_global2D = process(dir_sf, 512.0, 0)

# dir_sf = working_dir_path + r'final_data/2d/512run2D_disp_FFT/sf_par_perp_v_phi0_wrt_localF.txt'
# phi0_wrt_local2D = process(dir_sf, 512.0, 0)

# dir_sf = working_dir_path + r'final_data/2d/512run2D_disp_FFT/sf_par_perp_v_phi_wrt_localF.txt' #same as old sf file
# phi_wrt_local2D = process(dir_sf, 512.0, 0)

#PHI0 IMPROVEMENT
dir_sf = working_dir_path + r'phi0init/Runs/512_no_kpara/sf_par_perp_v_phi0_wrt_global_10_kparaF.txt'
phi0_nokpara = process(dir_sf, 512.0, 0)

dir_sf = working_dir_path + r'phi0init/Runs/512_3_2/sf_par_perp_v_phi0_wrt_global_10_kparaF.txt'
phi0_3_2 = process(dir_sf, 512.0, 0)

dir_sf = working_dir_path + r'phi0init/Runs/512_5_3/sf_par_perp_v_phi0_wrt_global_10_kparaF.txt'
phi0_5_3 = process(dir_sf, 512.0, 0)

dir_sf = working_dir_path + r'phi0init/Runs/1024_test/sf_par_perp_v_phi0_wrt_global_10_kparaF.txt'
phi0_1024 = process(dir_sf, 512.0, 0)

#reference straight line for GS95
lpar_temp =phi0_1024[0] #phi0_wrt_global_10kpara2[0]
lperp_temp =phi0_1024[1] #phi0_wrt_global_10kpara2[1]
count_temp =phi0_1024[-1] #phi0_wrt_global_10kpara2[-1]
ref_slope_3d_512_f = lpar_temp[0]*(np.power(lperp_temp[:count_temp],(2.0/3.0)))/(np.power(lperp_temp[0],(2.0/3.0)))

# lpar_temp = phi0_wrt_global2D[0]
# lperp_temp = phi0_wrt_global2D[1]
# count_temp = phi0_wrt_global2D[-1]
# ref_slope_2d = lpar_temp[5]*(np.power(lperp_temp[:count_temp],(2.0/3.0)))/(np.power(lperp_temp[6],(2.0/3.0)))

#plotting the structure functions
plt.figure(figsize=(5.0, 2.0), dpi=100)

plot(phi0_wrt_global_10kpara2, '512 Phi0 10kpara^-2')

plot(phi0_nokpara, '512 Phi0 initial')

plot(phi0_3_2,' 512 Phi0 10kpara^-(3/2)')

plot(phi0_5_3,' 512 Phi0 10kpara^-(5/3)')

plot(phi0_1024, '1024 Phi0 10kpara^-2')

#plot(phi0_wrt_global2D, 'PHI0 WRT GLOBAL FFT')

#plot(phi0_wrt_local2D, 'PHI0 WRT LOCAL FFT')

#plot(phi_wrt_local2D, 'PHI WRT LOCAL FFT')

#plot(phi_wrt_localalt, 'PHI WRT LOCAL REAL')

#plot(phi0_wrt_localalt, 'PHI0 WRT LOCAL REAL')

#plot(phi0_wrt_globalalt, 'PHI0 WRT GLOBAL REAL')


#plot(phi0_wrt_local_amp05, 'B M_A = 2.46')

#plot(phi0_wrt_local_amp1, 'B M_A = 4.93')

#plot(phi_wrt_local_amp05, 'PHI M_A = 2.46')

#plot(phi_wrt_local_amp1, 'PHI M_A = 4.93')

#plot(phi0_wrt_local_amp02, 'PHI0 B M_A = 0.985')
#plot(phi0_wrt_local_amp03, 'PHI0 B M_A = 1.477')
#plot(phi0_wrt_local_amp04, 'PHI0 B M_A = 1.970')

#plot(phi_wrt_local_amp02, 'PHI M_A = 0.985')
#plot(phi_wrt_local_amp03, 'PHI M_A = 1.477')
#plot(phi_wrt_local_amp04, 'PHI M_A = 1.970')

plt.plot(lperp_temp[:count_temp], 1.5*ref_slope_3d_512_f, lw=2.5, color = "black", ls = "-", label="GS95 2/3")
#plt.plot(lperp_temp[:count_temp], 1.5*ref_slope_2d, lw=2.5, color = "black", ls = "-", label="GS95 2/3")
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.002, 0.8)
plt.ylim(0.0015, 0.3)
plt.xlabel(r'$l_{\perp}/ L $ perpendicular',fontsize=11)
plt.ylabel(r'$l_{\parallel}/L $ parallel',fontsize=11)
plt.title('SF 3D PHI0 wrt global development')
plt.legend(loc='best',ncol=1,fontsize=11)

#plt.tight_layout()

plt.show()
