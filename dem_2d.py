import os.path
import numpy as np
import scipy.io as io
from astropy import units as u
import sunpy.map
from aiapy.calibrate import degradation
from aiapy.calibrate import register, update_pointing
import pickle
import time
import matplotlib.pyplot as plt
from demregpy import dn2dem
import utils as ut
import warnings
warnings.simplefilter('ignore')

def time_it(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        time_spent = end_time - start_time
        print(f"Time taken by '{func.__name__}': {time_spent:.6f} seconds")
        return result
    return wrapper

@time_it
def cal_dem_map(save_path, fov=None, target_file=None, tar_dir=None, savefile_name='DEM', scalefactor = 1.0, plot_res=False):
    """

    :param save_path:
    :param fov: fov = [[920, -280], [950, -250]] [[x,y],[x,y]]
    :param target_file: could be a file list that contain 6 aia files OR just one of them, the rest of the files will be
    found automaticlly
    :param tar_dir:
    :param savefile_name:
    :param scalefactor = 1.0(original resolution), #if scalefactor = 0.5,   0.6“ / pix   ----> 1.2"/pix:
    :return:
    """
    #This block can be modified to build your own file list------------------------
    if isinstance(target_file, list):
        file_list = target_file
        ctime = ut.read_time_from_header(file_list[0])
    else:
        file_list = [target_file]
        aia_kw_list = ['Z.131.', 'Z.171.', 'Z.211.', 'Z.94.', 'Z.335.', 'Z.193.']
        for ckw in aia_kw_list:
            if ckw in target_file: continue
            cfl = ut.makelist(tdir=tar_dir, keyword1=ckw,keyword2='fits')
            ctime = ut.read_time_from_header(target_file,exptime_correct=False)
            file_list.append(ut.find_closest_files(target_datetime=ctime, file_list=cfl))
    savefile_name += '_' + ctime.isot.replace('/','_')
    if fov is not None:
        savefile_name += '_fov'
        for c_coor in [item for sublist in fov for item in sublist]:
            savefile_name += '_'+str(c_coor)
    savefile_name += '.p'
    print('The input file list is: ', file_list)

    aia_tresp_en_file = ut.find_file_in_package('demregpy','tresp/aia_tresp_en.dat')
    savefile = os.path.join(save_path, savefile_name)
    if os.path.exists(savefile):
        print('file exists, exit now')
        return 1
    scalefactor = scalefactor #original resolution
    iter_st = time.time()
    amaps = sunpy.map.Map(file_list)
    # Get the wavelengths of the maps and get index of sort for this list of maps
    wvn0 = [m.meta['wavelnth'] for m in amaps]
    # print(wvn0)
    srt_id = sorted(range(len(wvn0)), key=wvn0.__getitem__)
    #print('sorted_list is : ', srt_id)
    amaps = [amaps[i] for i in srt_id]
    # print([m.meta['wavelnth'] for m in amaps])
    channels = [94, 131, 171, 193, 211, 335] * u.angstrom
    #cctime = atime.Time(time_string, scale='utc')
    nc = len(channels)
    degs = np.empty(nc)
    for i in np.arange(nc):
        degs[i] = degradation(channels[i], ctime, calibration_version=10)
    aprep = []
    for m in amaps:
        m_temp = update_pointing(m)
        aprep.append(register(m_temp))
    # Get the durations for the DN/px/s normalisation and
    # wavenlength to check the order - should already be sorted above
    wvn = [m.meta['wavelnth'] for m in aprep]
    durs = [m.meta['exptime'] for m in aprep]
    # Convert to numpy arrays as make things easier later
    durs = np.array(durs)
    #print(durs)
    wvn = np.array(wvn)
    worder = np.argsort(wvn)
    #print(worder)

    trin = io.readsav(aia_tresp_en_file)
    tresp_logt = np.array(trin['logt'])
    nt = len(tresp_logt)
    nf = len(trin['tr'][:])
    trmatrix = np.zeros((nt, nf))
    for i in range(0, nf):
        trmatrix[:, i] = trin['tr'][i]
    gains = np.array([18.3, 17.6, 17.7, 18.3, 18.3, 17.6])
    dn2ph = gains * np.array([94, 131, 171, 193, 211, 335]) / 3397.
    temps = np.logspace(5.7, 7.6, num=42)
    # Temperature bin mid-points for DEM plotting
    mlogt = ([np.mean([(np.log10(temps[i])), np.log10((temps[i + 1]))]) \
              for i in np.arange(0, len(temps) - 1)])

    # --------------------------------------------------------------------------------------------

    tmp_aprep = []
    for api, cap in enumerate(aprep):
        tmp_aprep.append(ut.make_sub_map(cur_map=cap, fov=fov))
    aprep = tmp_aprep
    #aprep[0].peek()

    cmap_shape = aprep[0].data.shape
    data_cube = np.zeros((cmap_shape[0], cmap_shape[1], 6))
    num_pix = (1 / scalefactor) ** 2
    #if scalefactor = 0.5,   0.6“ / pix   ----> 1.2"/pix
    #num_pix = 1.0
    rdnse = np.array([1.14, 1.18, 1.15, 1.20, 1.20, 1.18])*np.sqrt(num_pix)/num_pix
    for mi, m in enumerate(aprep):
        data_cube[:, :, mi] = m.data / degs[mi] / durs[mi]
        dn2ph_cube = np.broadcast_to(dn2ph, (cmap_shape[0], cmap_shape[1], 6))
        degs_cube = np.broadcast_to(degs, (cmap_shape[0], cmap_shape[1], 6))
        rdnse_cube = np.broadcast_to(rdnse, (cmap_shape[0], cmap_shape[1], 6))
        durs_cube = np.broadcast_to(durs, (cmap_shape[0], cmap_shape[1], 6))
    shotnoise = (dn2ph_cube * data_cube * num_pix) ** 0.5 / dn2ph_cube / num_pix / degs_cube
    edata_cube = (shotnoise ** 2 + rdnse_cube ** 2) ** 0.5 / durs_cube
    pe_time = time.time()
    print('it take to {} to prep in iter_time'.format(pe_time - iter_st))
    print('{} * {}, {} pixels to calculate in total'.format(data_cube.shape[0], data_cube.shape[1],
                                                            data_cube.shape[0] * data_cube.shape[1]))

    dem, edem, elogt, chisq, dn_reg = dn2dem(data_cube, edata_cube, trmatrix, tresp_logt, temps)
    if plot_res:
        fig = plt.figure(figsize=(8, 9))
        for j in range(10):
            fig = plt.subplot(3, 4, j + 1)
            plt.imshow(np.log10(dem[:, :, j * 4] + 1e-20), 'inferno', vmin=19, vmax=23, origin='lower')
            ax = plt.gca()
            ax.set_title('%.1f' % (5.6 + j * 2 * 0.1))
            fig = plt.subplot(3, 4, 12)
            plt.imshow(data_cube[:,:,0], origin='lower')
    res_tuple = (dem, edem, elogt, chisq, dn_reg, temps,fov)
    #pickle.dump(res_tuple, open('/Volumes/Data/20220511/dem/demreg/test_res.p', 'wb'))
    pickle.dump(res_tuple, open(savefile, 'wb'))
    return res_tuple

    #return 1


def main():
    aia_file = ''
    aia_dir = ''
    save_path = ''
    fov = [[905, -295], [965, -235]]
    cal_dem_map(save_path=save_path, target_file=aia_file, tar_dir=aia_dir,
                fov = fov, plot_res=False)

if __name__ == '__main__':
    main()