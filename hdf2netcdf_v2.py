import h5py
import numpy as np
from netCDF4 import Dataset
import sys

def wrapper(year, stage):
    datapath = '/bp1store/geog-tropical/data/ERA-Land/driving_data/'
    surfsolar = Dataset(datapath + str(year) + '_surface_net_solar_radiation.nc')
    lat = surfsolar.variables['latitude'][:]
    lon = surfsolar.variables['longitude'][:]
    hdf5_store = h5py.File(datapath + str(year)+'_cache.hdf5', 'r')

    if stage == 1:
        pet = hdf5_store["eto_val"][:4320, :, :] 
        filename = datapath + str(year)+ '_' + str(stage) +'_pet.nc'
        nc_write(pet, lat, lon, filename)
    else:
        pet = hdf5_store["eto_val"][4320:, :, :] 
        filename = datapath + str(year)+ '_'  + str(stage) +'_pet.nc'
        nc_write(pet, lat, lon, filename)
	

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
    pet_val.units='mm' 
    # values
    time[:] = np.arange(pet.shape[0])
    latitude[:] = lat
    longitude [:] = lon
    pet_val[:,:,:] = pet

    ds.close()
    
    return None    
if __name__=='__main__':
    wrapper(int(sys.argv[1]),int(sys.argv[2]))
