import struct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import gridspec


niter = 5

def smoothing(xarr) :
  lext = np.size(xarr,0) 
  sf2D_temp = np.zeros((lext,lext))
  for iter in range (0,niter) :
    for i in range (0,lext-1) :
      for j in range (0,lext-1) :
        if (i==0 and j==0) :
          sf2D_temp[i,j] = xarr[i,j]
          continue
        if (i==0) :
          sf2D_temp[i,j] = 0.5*xarr[i,j]+0.25*xarr[i,j-1]+0.25*xarr[i,j+1]
        elif (j==0) :
          sf2D_temp[i,j] = 0.5*xarr[i,j]+0.25*xarr[i-1,j]+0.25*xarr[i+1,j]
        else :
          sf2D_temp[i,j] = 0.5*xarr[i,j]+0.125*xarr[i-1,j]+0.125*xarr[i+1,j]+0.125*xarr[i,j-1]+0.125*xarr[i,j+1]
    xarr = sf2D_temp

  return xarr

working_dir_path = '/home/jonas/Documents/VSCode/DESY/'


sf2D_array = np.load(working_dir_path + 'final_data/3d/256run3D_FFT/sf2D_phi0_2000.npy')
#sf2D_array = np.load(working_dir_path + '3d_disp_mem/Runs/256_2nd_B/sf2D_phi0_2000.npy')
lent=256
lpar1 = 1.0*np.arange(lent/4)/lent
lperp1 = 1.0*np.arange(lent/4)/lent
sf2D1 = smoothing(sf2D_array)


sf2D_array = np.load(working_dir_path + '3d_disp_mem/Runs/256_2nd_B/sf2D_phi_2000.npy')
#sf2D_array = np.load(working_dir_path + '3d_disp_mem/Runs/256_2nd_B/sf2D_phi_2000.npy')
lent=256
lpar2 = 1.0*np.arange(lent/4)/lent
lperp2 = 1.0*np.arange(lent/4)/lent
sf2D2 = smoothing(sf2D_array)

sf2D_array = np.load(working_dir_path + 'final_data/3d/128run3D_FFT/sf2D_phi0_1000.npy')
#sf2D_array = np.load(working_dir_path + '3d_disp_mem/Runs/256_2nd_B/sf2D_phi0.npy')
#sf2D_array = np.load(working_dir_path + '3d_disp_mem/Runs/128_2nd_B/sf2D_phi0_1000.npy')
lent=128
lpar3 = 1.0*np.arange(lent/4)/lent
lperp3 = 1.0*np.arange(lent/4)/lent
sf2D3 = smoothing(sf2D_array)

#sf2D_array = np.load(working_dir_path + 'final_data/3d/256run3D_FFT/sf2D_phi.npy')
#sf2D_array = np.load(working_dir_path + '3d_disp_mem/Runs/256_2nd_B/sf2D_phi.npy')
sf2D_array = np.load(working_dir_path + 'final_data/3d/128run3D_FFT/sf2D_phi_1000.npy')
lent=128
lpar4 = 1.0*np.arange(lent/4)/lent
lperp4 = 1.0*np.arange(lent/4)/lent
sf2D4 = smoothing(sf2D_array)

plt.figure()
plt.imshow(sf2D2)
plt.title('256 2nd B PHI 2000')
fig = plt.gcf()
plt.clim()   # clamp the color limits
plt.colorbar()
plt.show()

# sf2D_array = np.load('../1024_runs/S1/data/decomped_modes/sf2D_v_A.npy')
# lent=1024
# lpar5 = 1.0*np.arange(lent/4)/lent
# lperp5 = 1.0*np.arange(lent/4)/lent
# sf2D5 = smoothing(sf2D_array)

# sf2D_array = np.load('../1024_runs/S4/data/decomped_modes/sf2D_v_A.npy')
# lent=1024
# lpar6 = 1.0*np.arange(lent/4)/lent
# lperp6 = 1.0*np.arange(lent/4)/lent
# sf2D6 = smoothing(sf2D_array)

# sf2D_array = np.load('../512_runs_4/C2/data/decomped_modes/sf2D_v_A.npy')
# lent=512
# lpar7 = 1.0*np.arange(lent/4)/lent
# lperp7 = 1.0*np.arange(lent/4)/lent
# sf2D7 = smoothing(sf2D_array)

# sf2D_array = np.load('../512_runs_4/CB1/data/decomped_modes/sf2D_v_A.npy')
# lent=512
# lpar8 = 1.0*np.arange(lent/4)/lent
# lperp8 = 1.0*np.arange(lent/4)/lent
# sf2D8 = smoothing(sf2D_array)


#fig=plt.figure()
fig = plt.figure(figsize=(20.0, 15.0))
# gs = gridspec.GridSpec(2, 4, hspace=0.0, wspace=0.1)
gs = gridspec.GridSpec(2, 2, hspace=0.4, wspace=0.0)

fig.suptitle('128 v 256 min 1000 1st B', size = 20)

ax0 = plt.subplot(gs[0],aspect='equal')
ax0.set_xlim(xmax=0.2)
#ax0.set_ylim(ymin=10.0)
ax0.set_ylim(ymax=0.2)
sff2d = sf2D1
lpar=lpar1
lperp=lperp1
#print(sff2d)
#cs=ax0.contour(lpar,lperp,np.transpose(sff2d),levels=[sff2d[3,3], sff2d[3,6], sff2d[3,9], sff2d[3,12],sff2d[3,15],sff2d[3,20], sff2d[3,25]],linewidths=2) #[5,5] just picks the sf value at that location - finds isocontours
cs=ax0.contour(lpar,lperp,np.transpose(sff2d),levels=[sff2d[3,1], sff2d[5,2],sff2d[5,4], sff2d[5,7],sff2d[5,10],sff2d[5,12], sff2d[5,15], sff2d[5,20], sff2d[5,25],sff2d[5,30], sff2d[5,35], sff2d[5,40]],linewidths=2)
plt.clabel(cs,inline=1,fontsize=10)
#ax0.axes.xaxis.set_ticklabels([])
ax0.set_title('256 Phi0 2000',fontsize=18)
ax0.set_ylabel(r'$l_{\perp}/L$',fontsize=18)
ax0.set_xlabel('$l_{\parallel}$',fontsize=16)
#ax0.legend(loc='lower left')


ax1 = plt.subplot(gs[1],aspect='equal')
ax1.set_xlim(xmax=0.2)
ax1.set_ylim(ymax=0.2)
sff2d = sf2D2
lpar=lpar2
lperp=lperp2
cs=ax1.contour(lpar,lperp,np.transpose(sff2d),levels=[sff2d[3,1], sff2d[5,2],sff2d[5,4], sff2d[5,7],sff2d[5,10],sff2d[5,12], sff2d[5,15], sff2d[5,20], sff2d[5,25],sff2d[5,30], sff2d[5,35], sff2d[5,40]], linewidths=2)
plt.clabel(cs,inline=1,fontsize=10)
#ax1.axes.xaxis.set_ticklabels([])
#ax1.axes.yaxis.set_ticklabels([])
ax1.set_title('256 Phi 2000',fontsize=18)
ax1.set_ylabel(r'$l_{\perp}$',fontsize=16)
ax1.set_xlabel('$l_{\parallel}$',fontsize=16)
# #ax0.legend(loc='bottom left')


ax2 = plt.subplot(gs[2],aspect='equal')
ax2.set_xlim(xmax=0.2)
#ax0.set_ylim(ymin=10.0)
ax2.set_ylim(ymax=0.2)
sff2d = sf2D3
lpar=lpar3
lperp=lperp3
cs=ax2.contour(lpar,lperp,np.transpose(sff2d),levels=[sff2d[3,1], sff2d[5,2],sff2d[5,4], sff2d[5,7],sff2d[5,10],sff2d[5,12], sff2d[5,15], sff2d[5,20]], linewidths=2)
plt.clabel(cs,inline=1,fontsize=10)
ax2.set_title('128 Phi0 1000',fontsize=18)
ax2.set_ylabel(r'$l_{\perp}$',fontsize=16)
ax2.set_xlabel('$l_{\parallel}$',fontsize=16)
#ax0.legend(loc='bottom left')
#plt.axvline(x=25.2)

ax3 = plt.subplot(gs[3],aspect='equal')
ax3.set_xlim(xmax=0.2)
#ax0.set_ylim(ymin=10.0)
ax3.set_ylim(ymax=0.2)
sff2d = sf2D4
lpar=lpar4
lperp=lperp4
cs=ax3.contour(lpar,lperp,np.transpose(sff2d),levels=[sff2d[3,1], sff2d[5,2],sff2d[5,4], sff2d[5,7],sff2d[5,10],sff2d[5,12], sff2d[5,15], sff2d[5,20]], linewidths=2)
plt.clabel(cs,inline=1,fontsize=10)
ax3.set_title('128 Phi 1000',fontsize=18)
ax3.set_ylabel(r'$l_{\perp}$',fontsize=16)
ax3.set_xlabel('$l_{\parallel}$',fontsize=16)
#ax0.legend(loc='bottom left')
#plt.axvline(x=25.2)



# ax4 = plt.subplot(gs[4],aspect='equal')
# ax4.set_xlim(xmax=0.2)
# #ax0.set_ylim(ymin=10.0)
# ax4.set_ylim(ymax=0.2)
# sff2d = sf2D5
# lpar=lpar5
# lperp=lperp5
# cs=ax4.contour(lpar,lperp,np.transpose(sff2d),levels=[sff2d[5,5],sff2d[5,10],sff2d[5,15],sff2d[5,30],sff2d[5,60],sff2d[5,100],sff2d[5,150]],linewidths=2)
# #plt.clabel(cs,inline=1,fontsize=10)
# ax4.set_title('(e) S1b v field',fontsize=18)
# ax4.set_ylabel(r'$l_{\perp}/L$',fontsize=18)
# ax4.set_xlabel('$l_{\parallel}/L$',fontsize=18)
# #ax0.legend(loc='bottom left')
# #plt.axvline(x=25.2)

# ax5 = plt.subplot(gs[5],aspect='equal')
# ax5.set_xlim(xmax=0.2)
# #ax0.set_ylim(ymin=10.0)
# ax5.set_ylim(ymax=0.2)
# sff2d = sf2D6
# lpar=lpar6
# lperp=lperp6
# cs=ax5.contour(lpar,lperp,np.transpose(sff2d),levels=[sff2d[5,5],sff2d[5,10],sff2d[5,15],sff2d[5,30],sff2d[5,60],sff2d[5,100],sff2d[5,150]],linewidths=2)
# #plt.clabel(cs,inline=1,fontsize=10)
# ax5.set_title('(f) S4b v field',fontsize=18)
# ax5.axes.yaxis.set_ticklabels([])
# #ax5.set_ylabel(r'$l_{\perp}$',fontsize=16)
# ax5.set_xlabel('$l_{\parallel}/L$',fontsize=18)
# #ax0.legend(loc='bottom left')
# #plt.axvline(x=25.2)

# ax6 = plt.subplot(gs[6],aspect='equal')
# ax6.set_xlim(xmax=0.2)
# #ax0.set_ylim(ymin=10.0)
# ax6.set_ylim(ymax=0.2)
# sff2d = sf2D7
# lpar=lpar7
# lperp=lperp7
# cs=ax6.contour(lpar,lperp,np.transpose(sff2d),levels=[sff2d[5,5],sff2d[5,10],sff2d[5,20],sff2d[5,40],sff2d[5,60],sff2d[5,80],sff2d[5,100]],linewidths=2)
# #plt.clabel(cs,inline=1,fontsize=10)
# ax6.set_title('(g) C2a v field',fontsize=18)
# #ax6.set_ylabel(r'$l_{\perp}$',fontsize=16)
# ax6.axes.yaxis.set_ticklabels([])
# ax6.set_xlabel('$l_{\parallel}/L$',fontsize=18)
# #ax0.legend(loc='bottom left')
# #plt.axvline(x=25.2)

# ax7 = plt.subplot(gs[7],aspect='equal')
# ax7.set_xlim(xmax=0.2)
# #ax0.set_ylim(ymin=10.0)
# ax7.set_ylim(ymax=0.2)
# sff2d = sf2D8
# lpar=lpar8
# lperp=lperp8
# cs=ax7.contour(lpar,lperp,np.transpose(sff2d),levels=[sff2d[5,5],sff2d[5,10],sff2d[5,20],sff2d[5,40],sff2d[5,60],sff2d[5,80],sff2d[5,100]],linewidths=2)
# #plt.clabel(cs,inline=1,fontsize=10)
# ax7.set_title('(h) CB1a v field',fontsize=18)
# ax7.axes.yaxis.set_ticklabels([])
# #ax7.set_ylabel(r'$l_{\perp}$',fontsize=16)
# ax7.set_xlabel('$l_{\parallel}/L$',fontsize=18)
# #ax0.legend(loc='bottom left')
# #plt.axvline(x=25.2)




#plt.savefig('test.eps', bbox_inches='tight',transparent=False)
plt.show()

