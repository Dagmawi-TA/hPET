import numpy as np
from netCDF4 import Dataset

def wrapper(startyear,endyear,latmin,latmax,lonmin,lonmax,regionname,t_resolution,output_path):
    """
    This is a wrapper function to run for downloading hPET and dPET data.
    All the arguments need to be given at the end of the script before running
    the script.

    :param startyear: the begging year to start the data download (min = 1981, max = 2019)
    :param endyear: the last year of data to be  downloaded (min = 1981, max = 2019)
    :param latmin: the minimum latitude value of the region (float)
    :param latmax: the maximum latitude value of the region (float)
    :param lonmin: the minimum longitude value of the region (float)
    :param lonmax: the maximum longitude value of the region (float)
    :param regionname: name of the region (it could be any name the user wants) (string)
    :param t_resolution: the time resolution to be downloaded (daily or hourly)
    :param output_path:  the file path to store the downloaded data (string)
    :return:
    """

    if t_resolution == 'daily':
        datapath = '/bp1store/geog-tropical/data/ERA-Land/driving_data/global_pet/daily_pet/'
    elif t_resolution == 'hourly':
        datapath = '/bp1store/geog-tropical/data/ERA-Land/driving_data/global_pet/hourly_pet/'
    else:
        raise ValueError("t_resolution is wrong please write 'daily' or 'hourly'")

    # set up the year array loop through each year to download the data
    years = np.arange(startyear,endyear+1)
    for y in range(0,len(years)):
        year=int(years[y])
        region_extract(datapath,year,latmin,latmax,lonmin,lonmax,regionname,t_resolution,output_path)
        print(year)



def region_extract(datapath,year,latmin,latmax,lonmin,lonmax,regionname,t_resolution,output_path):
    """
    This function extract the data from the global hPET and dPET file and write a new
    netCDF file with a file name <year>_<t_resolution>_pet_<regionname>.nc in the output_path
    provided.

    :param datapath: the file path where the hPET data is stored (url)
    :param year: the year for which data is going to be downloaded (integer)
    :param latmin: the minimum latitude value (float)
    :param latmax: the maximum latitude value (float)
    :param lonmin: the minimum longitude value (float)
    :param lonmax: the maximum longitude value (float)
    :param regionname: name of the region (it could be any name the user wants) (string)
    :param t_resolution: the time resolution to be downloaded (daily or hourly)
    :param output_path:  the file path to store the downloaded data (string)
    :return: hPET or dPET data in a netCDF file
    """

    if t_resolution == 'daily':
        fname = '_daily_pet.nc'
    elif t_resolution == 'hourly':
        fname = '_hourly_pet.nc'
    else:
        raise ValueError("t_resolution is wrong please write 'daily' or 'hourly'")

    pet_hr = Dataset(datapath + str(year) + fname)
    lats = pet_hr.variables['latitude'][:]
    lons = pet_hr.variables['longitude'][:]
    
    # extract the min and max index
    latminind, lonminind = nearest_point(latmin, lonmin, lats, lons)
    latmaxind, lonmaxind = nearest_point(latmax, lonmax, lats, lons)
 
    # read athe data pet
    reg_data=pet_hr.variables['pet'][:, latmaxind:latminind, lonminind:lonmaxind]  
    
    newlats=lats[latmaxind:latminind]
    newlons=lons[lonminind:lonmaxind]

    filename=output_path+str(year)+'_'+t_resolution+'_pet_'+regionname+'.nc'
    nc_write(reg_data, newlats, newlons, filename)

    return None
    

def nearest_point(lat_var, lon_var, lats, lons):
    """
    This function identify the nearest grid location index for a specific lat-lon
    point.
    :param lat_var: the latitude
    :param lon_var: the longitude
    :param lats: all available latitude locations in the data
    :param lons: all available longitude locations in the data
    :return: the lat_index and lon_index
    """
    # this part is to handle if lons are givn 0-360 or -180-180
    if any(lons > 180.0) and (lon_var < 0.0):
        lon_var = lon_var + 360.0
    else:
        lon_var = lon_var
        
    lat = lats
    lon = lons

    if lat.ndim == 2:
        lat = lat[:, 0]
    else:
        pass
    if lon.ndim == 2:
        lon = lon[0, :]
    else:
        pass

    index_a = np.where(lat >= lat_var)[0][-1]
    index_b = np.where(lat <= lat_var)[0][-1]

    if abs(lat[index_a] - lat_var) >= abs(lat[index_b] - lat_var):
        index_lat = index_b
    else:
        index_lat = index_a

    index_a = np.where(lon >= lon_var)[0][0]
    index_b = np.where(lon <= lon_var)[0][0]
    if abs(lon[index_a] - lon_var) >= abs(lon[index_b] - lon_var):
        index_lon = index_b
    else:
        index_lon = index_a

    return index_lat, index_lon


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
    latitude = ds.createVariable('latitude', np.float32, ('latitude',))
    longitude = ds.createVariable('longitude', np.float32, ('longitude',))
    pet_val = ds.createVariable('pet', 'f4', ('time','latitude','longitude'))
   
    time.units = 'days since 1981-01-01'
    time.calendar = 'proleptic_gregorian'
    time[:] = np.arange(pet.shape[0])
    latitude[:] = lat
    longitude [:] = lon
    pet_val[:,:,:] = pet
    ds.close()
    
    return None    
##### ---------------------------######
# to extract data please fill the arguments in the function and run the script.
if __name__ == '__main__':
    # example (please change these values to your specification)
    # input arguments
    startyear = 1981
    endyear = 1983
    latmin=-0.7
    latmax=1.7
    lonmin=35.8
    lonmax=38.4
    regionname='kenya'
    t_resolution ='daily'
    output_path = '/home/fp20123/andres_kenya/'
    # run script
    wrapper(startyear,endyear,latmin,latmax,lonmin,lonmax,regionname,t_resolution,output_path)
