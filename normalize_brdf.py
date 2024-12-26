__author__ = "Shi Qiu"

# This script is used to normalize the BRDF of Landsat-8 to the ideal solar angles.
# The solar zenith angle is modeled as a function of the day of the year and latitude for the sun-synchronous orbit, derived from an intermediate variable known as the local overpass time. 
# The local overpass time is calculated using Eq. (4) from Li et al. 2019, based on an astronomical model parameterized by the local overpass time at the Equator, satellite inclination, and latitude.
# We do not provide the brdf core model in python, but you can request it from Dr. Zhang (https://github.com/hankui).
# Users need to intsall the libraries needed in this script by self.

import os
from pathlib import Path
import numpy as np
import glob
import click
import pysolar
import model
import pandas
from datetime import datetime, timezone


def model_reflectance(vza, sza, raa, b):
    # parameters for the model
    # _pars_array = [[ 774, 372, 79], \
	# 		  [1306, 580,178], \sjobs
	# 		  [1690, 574,227], \
	# 		  [3093,1535,330], \
	# 		  [3430,1154,453], \
	# 		  [2658, 639,387], \
	# 		  ]/10000.0 # convert 0 to 1
    _pars_array = [
        [0.0774, 0.0372, 0.0079],
        [0.1306, 0.0580, 0.0178],
        [0.1690, 0.0574, 0.0227],
        [0.3093, 0.1535, 0.0330],
        [0.3430, 0.1154, 0.0453],
        [0.2658, 0.0639, 0.0387]
    ]
    # This is the core model of the BRDF model, which is not provided here.
    # ask Dr. Zhang for the core model
    kk = model.Kernels(vza, sza, raa, doIntegrals=False, RossHS=False, RossType='Thick', LiType='Sparse', normalise=1, RecipFlag=True, MODISSPARSE=True)
    K = np.ones([3,len(vza)])
    K[1,:] = kk.Ross[:]
    K[2,:] = kk.Li[:]
    # print("K:", np.transpose(K).shape, np.array(_pars_array[b]).shape, type(np.array(np.transpose(_pars_array[b]))))
    fwd = np.matmul(np.transpose(K), np.array(_pars_array[b]))
    # model.Kernels.printKernels(kk,header=True,reflectance=False,file="tempfile");
    return fwd

def getNadirLocalTime(lat, localtime = 10.5, inclination = 98.2):
    # mean latitudinal overpass time 10.5 hours for Sentinel-2
    # see eq. 4 in Li, Z., Zhang, H.K., Roy, D.P., Investigation of Sentinel-2 bidirectional reflectance hot-spot sensing conditions. IEEE Transactions on Geoscience and Remote Sensing. ​In press.	double Mlat_hours; /* mean latitudinal overpass time*/
    if lat > (180 - inclination): # high latitudes, we use the degree cloest to the inclination, suggested by Dr. Zhang
        _lat = (180 - inclination)
        print("Warning: the degree cloest to the inclination {%f} was used for the high altitude" % _lat)
    elif lat < -(180 - inclination):
        _lat = -(180 - inclination)
        print("Warning: the degree cloest to the inclination {%f} was used for the high altitude" % _lat)
    else:
        _lat = lat
    return localtime - np.rad2deg(np.arcsin(np.tan(np.deg2rad(_lat))/np.tan(np.deg2rad(inclination)))) / 15  # in hours for descending orbit

def convertLocalTime2GMT(dLongitude, doy, time):
    gmt_time = time - dLongitude / 15  # convert local time to GMT
    iDay_gmt = doy

    if gmt_time > 24:
        iDay_gmt += 1
        gmt_time -= 24
        if iDay_gmt > 365:
            iDay_gmt = 365
    elif gmt_time < 0:
        iDay_gmt -= 1
        gmt_time += 24
        if iDay_gmt < 1:
            iDay_gmt = 1  # assuming non-leap year for simplicity

    dHours_gmt = int(gmt_time)
    left_mins = (gmt_time - dHours_gmt) * 60
    dMinutes_gmt = int(left_mins)
    dSeconds_gmt = (left_mins - dMinutes_gmt) * 60


    return iDay_gmt, "{:d}:{:d}:{:f}".format(dHours_gmt, dMinutes_gmt, dSeconds_gmt)


@click.command()
@click.option("--ci",       "-i", type=int, help="The core's id",               default=4)
@click.option("--cn",       "-n", type=int, help="The number of cores",         default=200)
@click.option("--des",      "-s", type=str, help="The directory of working on", default='/gpfs/sharedfs1/zhulab/Shi/ProjectGlobalGreening/GlobalPathRow')
@click.option("--target",   "-t", type=str, help="The directory of working on", default='SurfaceReflectanceIdealAngle')
def main(ci, cn, des, target):
    print('#######################################################')
    print('Working on appending ideal solar angles on the table')
    print('ci: {}'.format(ci))
    print('cn: {}'.format(cn))
    print('des: {}'.format(des))
    print('target: {}'.format(target))
    
    targetout = target + "BRDF"
    
    # check out the mean local time of Landsat8 and Sentinel2, from https://github.com/hankui/landsat_sentinel_solar_angle/blob/main/Pro.Landsat.Sentinel-2.sz.r
    # for Landsat-8 only
    l8_localtime = 10.18333333333
    l8_inclination = 98.2

    # for Sentinel-2 only
    s2_inclination = 98.62
    s2_localtime = 10.5

    # for HLS
    hls_localtime = (s2_localtime+l8_localtime)/2
    hls_inclination = (s2_inclination+l8_inclination)/2
    
    
    # read global grid tiles
    tiles_list = glob.glob(os.path.join(des, 'p*'))
    tiles_list = sorted(tiles_list)

    # Process each task
    tasks_range = range(ci - 1, len(tiles_list), cn)
    for i in tasks_range:
        dir_tile = tiles_list[i]
        tile_name = Path(dir_tile).stem
        
        print("***********************************************************")
        print("Working for {}".format(tile_name))
        
        dir_data = os.path.join(dir_tile, 'LandsatARD', target)
        files_list = glob.glob(os.path.join(dir_data, 'p*.csv'))
        for file_csv in files_list:
            file_csv_out = os.path.join(dir_tile, 'LandsatARD', targetout, Path(file_csv).name)
            Path(file_csv_out).parent.mkdir(parents=True, exist_ok=True)
            if os.path.exists(file_csv_out):
                print("Skip the file (data exist): {}".format(file_csv_out))
                continue
            dt_all = pandas.read_csv(file_csv)
            if len(dt_all) == 0:
                print("Skip the file (data empty): {}".format(file_csv_out))
                continue
            # in case when we have produced the angles
            if 'isza' not in dt_all.columns:
                # convert to datatime and julian date
                dt_all['datetime'] = pandas.to_datetime(dt_all.year * 1000 + dt_all.doy, format='%Y%j')
                # dt_all['juliandate']= dt_all.datetime.apply(lambda x: sum(gcal2jd(x.year, x.month, x.day))) + 0.5 # plus + 0.5, for locating 12:00:00, and go through with int
                # only remain the year, month, day since we will merge landsat and modis according to juliandate, which are same year, month, and day.
                dt_all['month']    = dt_all.datetime.apply(lambda x: x.month)
                dt_all['day']      = dt_all.datetime.apply(lambda x: x.day)

                for index, row in dt_all.iterrows():
                    lat_localtime = getNadirLocalTime(row['lat'], l8_localtime, l8_inclination)
                    if np.isnan(lat_localtime):
                        print("Error: lat_localtime is nan")
                        print(row)
                        continue

                    [_, lat_gmt] = convertLocalTime2GMT(row['lon'], row['doy'], lat_localtime)
                    # UTC is effectively the new name for GMT. It has very minor differences, but none that will impact you in that scenario.
                    # local time to UTC
                    row_date = datetime.strptime('{:04d}-{:02d}-{:02d} {} UTC+0000'.format(row['year'], row['month'], row['day'], lat_gmt), '%Y-%m-%d %H:%M:%S.%f %Z%z')
                    row_date= row_date.replace(tzinfo=timezone.utc)
                    dt_all.loc[index,'isza'] = 90 - pysolar.solar.get_altitude(row['lat'] , row['lon'], row_date)
                    dt_all.loc[index,'isaa'] = pysolar.solar.get_azimuth(row['lat'] , row['lon'], row_date) - 360 # make it same as to Landsat meta

            # BRDF normalization
            # delete the rows with nan in the ideal solar angles
            dt_all = dt_all.dropna(subset=['isza'])
            solar_zenith_output = dt_all['isza']
            for band_id, band_name in enumerate(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']):
                
                srf1 = model_reflectance(dt_all['vza'], dt_all['sza'], dt_all['vaa'] - dt_all['saa'], band_id)
                srf0 = model_reflectance(np.ones(len(dt_all['vza']))*0, solar_zenith_output, np.ones(len(dt_all['vza']))*0, band_id)
                dt_all[band_name] = ((srf0/srf1*dt_all[band_name]))
            # save the file
            dt_all.to_csv(file_csv_out + ".part", index=False)
            # rename the file
            os.rename(file_csv_out + ".part", file_csv_out)
        
if __name__ == "__main__":
    main()
