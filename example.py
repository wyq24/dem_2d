from dem_2d import cal_dem_map
import pickle
import matplotlib.pyplot as plt
import numpy as np

# aia_file = '/Users/walterwei/Downloads/work/20221029/New/aia.lev1_euv_12s.2022-10-29T190910Z.171.image_lev1.fits'  # file list of all 6 aia files or just the path of one of them.
# aia_dir = '/Users/walterwei/Downloads/work/20221029/New/'
# save_path = '/Users/walterwei/Downloads/work/20221029/save/'
# #fov = [[905, -295], [965, -235]]
# fov = [[-250, 320], [-150, 420]]
# cres = cal_dem_map(save_path=save_path, target_file=aia_file, tar_dir=aia_dir,
#             fov=fov, plot_res=False)

with open('/Users/walterwei/Downloads/work/20221029/save/DEM_2022-10-29T19:09:09.352_fov_-250_320_-150_420.p', 'rb') as csf:
#with open('/Users/walterwei/Downloads/work/20221029/save/DEM_2022-10-29T18:54:57.630_fov_-250_320_-150_420.p','rb') as csf:
        dem_res = pickle.load(csf)
csf.close()
with open('/Users/walterwei/Downloads/work/20221029/save/DEM_2022-10-29T18:54:57.630_fov_-250_320_-150_420.p','rb') as csf2:
        dem_res2 = pickle.load(csf2)
csf2.close()

print('')
from utils import temperature_em_map
#tmap, em_map= temperature_em_map(dem_res, [8.e6, 2.5e7])
tmap, em_map= temperature_em_map(dem_res)
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 6),sharex=True, sharey=True)

axs[2].imshow(dem_res[3][:,:], origin='lower')
axs[0].imshow(tmap, origin='lower')
axs[1].imshow(em_map, origin='lower')

# kw_list = ['dem', 'edem', 'elogt', 'chisq', 'dn_reg', 'temps']
# for kwi, ckw in enumerate(kw_list):
#         axs.flat[kwi+2].imshow()

#==plot response=======================
# cx, cy = [82,102]
# cx, cy = [94,109]
# fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))
# axs[1].imshow(dem_res[4][:,:,3], origin='lower')
# axs[1].plot(cx,cy, marker='x', color='r')
# axs[1].set_title('AIA193')
# axs[1].set_ylabel('Y_pixel')
# axs[1].set_xlabel('X_pixel')
#
# axs[0].errorbar(np.log10(dem_res[5][1:]), dem_res[0][cy, cx,:], xerr=dem_res[2][cy, cx,:], yerr = dem_res[1][cy, cx,:], label='Jet')
# axs[0].errorbar(np.log10(dem_res2[5][1:]), dem_res2[0][cy, cx,:], xerr=dem_res2[2][cy, cx,:], yerr = dem_res2[1][cy, cx,:], label = 'Pre- Jet')
# axs[0].set_yscale("log")
# axs[0].set_xlabel('T[k] log_10')
# axs[0].set_ylabel('DEM')
# axs[0].set_title('DEM result comparison')
# axs[0].legend()
# axs[0].set_ylim([1.e17, 5.e21])
#plt.xscale("log")
#=============================
# fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(10, 6))
# wlist = [94, 131, 171, 193, 211, 335]
# for wi in range(6):
#         axs.flat[wi].imshow(dem_res[4][:,:,wi], origin='lower')
#         axs.flat[wi].set_title(str(wlist[wi]))

plt.show()