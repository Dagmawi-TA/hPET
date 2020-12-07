#---------------------------------------------------------#
# This module contain the functions to aggregate the hourly
# PET values to daily PET values.


# Dagmawi Teklu Asfaw
# July, 2020
#----------------------------------------------------------#
import sys
import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt
import datetime as dt
import h5py
# ---------------------------------------------------------#

def wrapper(year, stage):
    datapath = '/bp1store/geog-tropical/data/ERA-Land/driving_data/' # data path
    daily_pet, lat, lon = read_hourly_pet(datapath, year, stage)

    # 3. write the output to a netcdf file
    filename=datapath+str(year)+'_'+str(stage)+'_dailyPET.nc'
    nc_write(daily_pet, lat, lon, filename)



def read_hourly_pet(datapath, year, stage):
    # read the offset
    offsetdata=Dataset(datapath + 'timezone_offset.nc')
     # # read hourly pet
    data = h5py.File(datapath + str(year)+'_cache.hdf5', 'r')
    # read previous and next year pet (only one day data)
    # previous year
    data_prev = h5py.File(datapath + str(year-1)+'_cache.hdf5', 'r')
       # next year
    data_next = h5py.File(datapath + str(year+1)+'_cache.hdf5', 'r')

    # before we make the daily sum for pet we need to account for 
    # the time offset. some part the time is offset is negative and
    # in some part it is positive while there are places with no time
    # offset from the UTC. Therefore, we need to add hours and subtract hours
    # from the begning and end of each pet timeseries to make all values 
    # similar on the day which is (00:00-23:59).

    # the data is too large hence it is divided into 2 stages by area
    # looping through the time is very complex hence this approach is taken.

    # loop through the days 168 at a time to overcome memory issues
    if year%4 == 0:
        dlen = 366 # data length leap year
    else:
        dlen = 365 # data length non leap year
    
    
    if stage == 1: # left side off the globe
        
        hourly_pet = data["eto_val"][:, :, :1800] 
        offset=offsetdata.variables['offset'][:,:1800]
        latitude = offsetdata.variables['latitude'][:]
        longitude = offsetdata.variables['longitude'][:1800]
        previous_petdata = data_prev["eto_val"][-128:, :, :1800]  # 128 hours of data to handle the offsets in the ocean (-127)
        next_petdata = data_next["eto_val"][:128, :, :1800]
        daily_pet=np.zeros((dlen,hourly_pet.shape[1],hourly_pet.shape[2]))
        
        # 1. for negative offset
        # negative and zero offset
        for j in range(0,hourly_pet.shape[1]):
            for k in range(0,hourly_pet.shape[2]):
                offsetval=offset[j,k]
                if offsetval<=0:
                    offsetval=abs(offsetval) # need the absolute value to index the array
                    x1=np.concatenate((hourly_pet[offsetval:,j,k], next_petdata[:offsetval,j,k]))
                    x1=np.reshape(x1,(dlen,24)) # reshape the array (365,24) for normal year and (366,24) for leap year
                    daily_pet[:,j,k] = np.sum(x1,axis=1)
                else:
                    # positive offset
                    x2=np.concatenate((previous_petdata[-offsetval:,j,k], hourly_pet[:(hourly_pet.shape[0]-offsetval),j,k]))
                    x2=np.reshape(x2,(dlen,24))
                    daily_pet[:,j,k] = np.sum(x2,axis=1)
                   

    elif stage == 2: # right side of the globe
        hourly_pet = data["eto_val"][:, :, 1800:] 
        offset=offsetdata.variables['offset'][:,1800:]
        latitude = offsetdata.variables['latitude'][:]
        longitude = offsetdata.variables['longitude'][1800:]
        previous_petdata = data_prev["eto_val"][-128:, :, 1800:]  # 128 hours of data to handle the offsets in the ocean (-127)
        next_petdata = data_next["eto_val"][:128, :, 1800:]
        daily_pet=np.zeros((dlen,hourly_pet.shape[1],hourly_pet.shape[2]))
        # 1. for negative offset
        # negative and zero offset
        for j in range(0,hourly_pet.shape[1]):
            for k in range(0,hourly_pet.shape[2]):
                offsetval=offset[j,k]
                if offsetval<=0:
                    offsetval=abs(offsetval) # need the absolute value to index the array
                    x1=np.concatenate((hourly_pet[offsetval:,j,k], next_petdata[:offsetval,j,k]))
                    x1=np.reshape(x1,(dlen,24))
                    daily_pet[:,j,k] = np.sum(x1,axis=1)
                else:
                    # positive offset
                    x2=np.concatenate((previous_petdata[-offsetval:,j,k], hourly_pet[:(hourly_pet.shape[0]-offsetval),j,k]))
                    x2=np.reshape(x2,(dlen,24)) # reshape the array (365,24) for normal year and (366,24) for leap year
                    daily_pet[:,j,k] = np.sum(x2,axis=1)
    
    else:
        print('Wrong stage input please put 1 or 2!')
    
    return daily_pet, latitude, longitude


def nc_write(pet, lat, lon, filename):
    """
    this function write the PET on a netCDF file.

    :param: pet: PET (time,lat,lon)
    :param: lat: latitude
    :param:lon: longitude
    :param:filename: the file name to write the values with .nc extension

    :return:  produce a netCDF file in the same directory.
    """

    ds = Dataset(filename, mode='w', format='NETCDF4_CLASSIC')

    time = ds.createDimension('time', None)
    latitude = ds.createDimension('latitude', len(lat))
    longitude = ds.createDimension('longitude', len(lon))
   
    time = ds.createVariable('time', np.float32, ('time',))
    latitude = ds.createVariable('latitude', np.float32, ('latitude',))  #, fill_value=-32767
    longitude = ds.createVariable('longitude', np.float32, ('longitude',))
    pet_val = ds.createVariable('pet', 'f4', ('time','latitude','longitude'),zlib=True)
    
    # units
    time.units='hours since 1981-01-01 00:00:00.0'
    time.calendar='proleptic_gregorian'
    latitude.units='degree_north'
    longitude.units='degree_east'
    pet_val.units='mm/day' 
    # values
    time[:] = np.arange(pet.shape[0])
    latitude[:] = lat
    longitude [:] = lon
    pet_val[:,:,:] = pet

    ds.close()

# -------------------------------------------------- #
if __name__=='__main__':
    wrapper(int(sys.argv[1]),int(sys.argv[2]))

