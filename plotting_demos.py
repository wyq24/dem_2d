import pickle
import matplotlib.pyplot as plt
import numpy as np

#cx, cy = [82,102]
cx, cy = [94,109]
cfile1 = '/Users/walterwei/Downloads/work/20221029/save/DEM_2022-10-29T19:09:09.352_fov_-250_320_-150_420.p'
cfile2 = '/Users/walterwei/Downloads/work/20221029/save/DEM_2022-10-29T18:54:57.630_fov_-250_320_-150_420.p'


with open(cfile1, 'rb') as csf1:
#with open('/Users/walterwei/Downloads/work/20221029/save/DEM_2022-10-29T18:54:57.630_fov_-250_320_-150_420.p','rb') as csf:
        dem_res1 = pickle.load(csf1)
csf1.close()
with open(cfile2,'rb') as csf2:
        dem_res2 = pickle.load(csf2)
csf2.close()
#plot======================================================
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(20, 6))
axs[1].imshow(dem_res2[4][:,:,3], origin='lower')
axs[2].imshow(dem_res1[4][:,:,3], origin='lower')
tlist = ['18:54:57(pre-jet)', '19:09:09(jet)']
clist = ['orange','r']
for axi in [1,2]:
        axs[axi].plot(cx,cy, marker='x', color=clist[axi-1])
        axs[axi].set_title('AIA193'+tlist[axi-1])
        axs[axi].set_ylabel('Y_pixel')
        axs[axi].set_xlabel('X_pixel')

axs[0].errorbar(np.log10(dem_res1[5][1:]), dem_res1[0][cy, cx,:], xerr=dem_res1[2][cy, cx,:], yerr = dem_res1[1][cy, cx,:], color = clist[1], label='Jet')
axs[0].errorbar(np.log10(dem_res2[5][1:]), dem_res2[0][cy, cx,:], xerr=dem_res2[2][cy, cx,:], yerr = dem_res2[1][cy, cx,:], color = clist[0], label = 'Pre- Jet')
axs[0].set_yscale("log")
axs[0].set_xlabel('T[k] log_10')
axs[0].set_ylabel('DEM')
axs[0].set_title('DEM result comparison')
axs[0].legend()
axs[0].set_ylim([1.e17, 5.e21])
plt.show()