from astropy.time import Time
import os
from astropy.io import fits
import time
import re
from datetime import datetime, timedelta
import pkg_resources
from astropy.coordinates import SkyCoord
from astropy import units as u
import pickle

def time_it(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        time_spent = end_time - start_time
        print(f"Time taken by '{func.__name__}': {time_spent:.6f} seconds")
        return result
    return wrapper


def find_file_in_package(package, filename):
    try:
        # Get the path of the file
        file_path = pkg_resources.resource_filename(package, filename)
        return file_path
    except:
        return None

def read_time_from_header(fits_file, exptime_correct=False):
    """
    Date and time when observation of this image STARTED!
    :param fits_file:
    :return: if the date is not in normal form, e.g. EUI image which has 'May' in it, return the original string
    """
    format_list = ['DATE-OBS', 'date_obs', 'T_OBS']
    hdulist = fits.open(fits_file, mode='readonly')
    for chdu in hdulist:
        for cformat in format_list:
            try:
                result = chdu.header[cformat]
                cexpt = chdu.header['exptime']
                #print(result)
                try:
                    if exptime_correct:
                        return Time(result)+timedelta(seconds=cexpt/2.0)
                    else:
                        return Time(result)
                except:
                    return result
            except Exception as e:
                pass
    return False

def extract_datetime_from_filename(filename):
    # Define the patterns for the date and time format in the filename
    if os.path.isabs(filename):
        filename = os.path.basename(filename)
    patterns = [
        (r'\d{4}-\d{2}-\d{2}T\d{6}Z', '%Y-%m-%dT%H%M%SZ'),  # YYYY-MM-DDTHHMMSSZ
        (r'\d{8}T\d{9}', '%Y%m%dT%H%M%S%f')  # YYYYMMDDTHHMMSSZ
    ]
    for pattern, format in patterns:
        # Use regex to find the pattern in the filename
        match = re.search(pattern, filename)
        if match is not None:
            # Extract the date and time string
            datetime_str = match.group()
            # Convert the date and time string to a datetime object
            datetime_obj = datetime.strptime(datetime_str, format)
            return Time(datetime_obj)
    # If no patterns matched
    return None

def makelist(tdir='', keyword1='', keyword2='', exclude=None):
    li = []
    # for root, dirs, files in os.walk(tdir):
    root = os.getcwd()
    files = os.listdir(tdir)
    for file in files:
        if exclude is None:
            if keyword1 in file and keyword2 in file:
                li.append(os.path.join(tdir, file))
        else:
            if keyword1 in file and keyword2 in file and exclude not in file:
                li.append(os.path.join(tdir, file))

    return li

def find_closest_files(target_datetime, file_list, use_filename=True, number_of_results=3, verify_by_header=True, **kwargs):
    # Extract datetime and calculate time difference for each file
    time_diffs = []
    for filename in file_list:
        if use_filename:
            file_datetime = extract_datetime_from_filename(filename)
        else:
            if verify_by_header:
                verify_by_header = False
                #warnings('Its nonsence to verify the header time when you are read time from the header')
            file_datetime = read_time_from_header(filename, **kwargs)
        if file_datetime is not None:
            time_diff = abs(target_datetime - file_datetime)
            time_diffs.append((time_diff, filename))

    # Sort by time difference and select the three files with smallest difference
    closest_files = sorted(time_diffs, key=lambda x: x[0])[:number_of_results]
    if verify_by_header:
        time_diffs_header = []
        for diff, filename in closest_files:
            file_datetime = read_time_from_header(filename, **kwargs)
            if file_datetime is not None:
                time_diff = abs(target_datetime - file_datetime)
                time_diffs_header.append((time_diff, filename))
        closest_files_header = sorted(time_diffs_header, key=lambda x: x[0])
        #print(closest_files[0][1], closest_files_header[0][1])
        # if closest_files[0][1] == closest_files_header[0][1]:
        #     return closest_files[0][1]
        # else:
        #     raise ValueError('Time from the file name is not reliable!')
        return closest_files_header[0][1]

    # Return just the filenames
    return closest_files[0][1]

def make_sub_map(cur_map, fov):
    bole = SkyCoord(fov[0][0] * u.arcsec, fov[0][1] * u.arcsec, frame=cur_map.coordinate_frame)
    tori = SkyCoord(fov[1][0] * u.arcsec, fov[1][1] * u.arcsec, frame=cur_map.coordinate_frame)
    sub_map = cur_map.submap(bole, top_right=tori)
    return sub_map


def temperature_em_map(dem_res_save, t_range=None):
    import numpy as np
    if isinstance(dem_res_save, str):
        with open(dem_res_save, 'rb') as csf:
            dem_res = pickle.load(csf)
        csf.close()
    else:
        dem_res = dem_res_save
    mean_t = (dem_res[5][1:]+dem_res[5][0:-1])/2.0
    if t_range is not None:
        t_min_idx = np.nanargmin(abs(mean_t- t_range[0]))
        t_max_idx = np.nanargmin(abs(mean_t- t_range[1]))
    else:
        t_min_idx, t_max_idx = (0,len(mean_t))


    #weights_sum = np.sum(dem_res[0][:,:,t_min_idx:t_max_idx]/dem_res[1][:,:,t_min_idx:t_max_idx], axis=2)
    #weighted_indices_sum = np.sum(dem_res[0][:,:,t_min_idx:t_max_idx]/dem_res[1][:,:,t_min_idx:t_max_idx] * mean_t[np.newaxis, np.newaxis, t_min_idx:t_max_idx], axis=2)
    #error = dem_res[1][:,:,t_min_idx:t_max_idx]/dem_res[0][:,:,t_min_idx:t_max_idx]
    ##todo: uncertainty on temperature: https://mingjiejian.github.io/2019/12/06/error-weighted-mean/
    error = np.ones_like(dem_res[1][:,:,t_min_idx:t_max_idx])/dem_res[1][:,:,t_min_idx:t_max_idx]
    weights_sum = np.sum(dem_res[0][:,:,t_min_idx:t_max_idx]*error, axis=2)
    weighted_indices_sum = np.sum(dem_res[0][:,:,t_min_idx:t_max_idx]*error * mean_t[np.newaxis, np.newaxis, t_min_idx:t_max_idx], axis=2)
    weighted_indices_sum_all = np.sum(dem_res[0][:, :, :] * mean_t[np.newaxis, np.newaxis, :], axis=2)
    centroid_indices = weighted_indices_sum / weights_sum
    return centroid_indices, weighted_indices_sum_all
