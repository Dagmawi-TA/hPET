# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 13:17:35 2020

@author: fp20123
"""

#---------------------------------------------------------#
# This module contain all the functions to calculate
# the Potential Evapotranspiration (PET).
# The method used are explained in FAO reference mannual
# which can be acessed at this link
# http://www.fao.org/docrep/X0490E/x0490e05.htm
# The procedure uses the Penman Monteith method to calculate
# PET from climatic variables.
#
# Dagmawi Teklu Asfaw
# March, 2020
#----------------------------------------------------------#
import sys
import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt
import datetime as dt
import h5py
# ---------------------------------------------------------#

def wrapper(year):

    # 1. read all the netcdf files from ERA5 pre_erapev_m,  pre_erapev_m,
    datapath = '/bp1store/geog-tropical/data/ERA-Land/driving_data/' # data path
    
    latitude,longitude,pre_surface_net_solar_radiation_J_m2, pre_surface_net_thermal_radiation_J_m2, \
    surfsolar,surfthermal,tmean,tdew,surfpres,uwnd,vwnd,erapet=read_data(year, datapath)

    print('Data reading done!')

    # 2 extract all the data in array form for each week and calculate the PET
    PET_mm_hr = data_array(datapath,year,latitude,longitude,pre_surface_net_solar_radiation_J_m2,
        pre_surface_net_thermal_radiation_J_m2, surfsolar, surfthermal, tmean, tdew, surfpres, uwnd, vwnd, erapet)

    # 3. write the output to a netcdf file
    # data is too large to be writtenn on RAM hence write as hdf file on hard drive
    # and change to netcdf file later !!!!!
    #fname=datapath+str(year)+'_PET.nc'
    #nc_write(PET_mm_hr, latitude, longitude, fname)

    print('PET writting is done!')


    # NB this is for further work
    # Daily PET
    # 5. run  hourly2daily()
    # 6. run calc_pet() with pet_time='daily'
    # 7. write the output to a netcdf file

    return None


def read_data(year, datapath):
    """
    This is the function to read each year data for all the variables.
    :param: year: the year of the data to be read.
    :return: latitude, longitude and the other 7 variables
    """

    # previous year data reading  pre_erapev_m,
    # 1. previous year data , pre_erapev_m
    pre_surface_net_solar_radiation_J_m2, pre_surface_net_thermal_radiation_J_m2 = read_data_previous(year,datapath)
    
    # 2. read all the netcdf files from ERA5
    surfsolar = Dataset(datapath + str(year) + '_surface_net_solar_radiation.nc')
    surfthermal = Dataset(datapath + str(year) + '_surface_net_thermal_radiation.nc')
    tmean = Dataset(datapath + str(year) + '_2m_temperature.nc')
    tdew = Dataset(datapath + str(year) + '_2m_dewpoint_temperature.nc')
    surfpres = Dataset(datapath + str(year) + '_surface_pressure.nc')
    uwnd = Dataset(datapath + str(year) + '_10m_u_component_of_wind.nc')
    vwnd = Dataset(datapath + str(year) + '_10m_v_component_of_wind.nc')
    erapet = Dataset(datapath + str(year) + '_potential_evaporation.nc')

    latitude = surfsolar.variables['latitude'][:]
    longitude = surfsolar.variables['longitude'][:]

    return latitude,longitude,pre_surface_net_solar_radiation_J_m2, pre_surface_net_thermal_radiation_J_m2, \
    surfsolar,surfthermal,tmean,tdew,surfpres,uwnd,vwnd,erapet


def read_data_previous(year, datapath):
    """
    This is the function to read each year data for all the variables.
    :param: year: the year of the data to be read.
    :return: latitude, longitude and the other 7 variables
    """
    # channge the year (string) to intiger and subtruct 1 year
    print(year)
    lastyear = year - 1
    # 1. read all the netcdf files from ERA5
    surfsolar = Dataset(datapath + str(lastyear) + '_surface_net_solar_radiation.nc')
    surfthermal = Dataset(datapath + str(lastyear) + '_surface_net_thermal_radiation.nc')
    surface_net_solar_radiation_J_m2 = surfsolar.variables['ssr'][-2:, :, :]  # (time,latitude,longitude)
    surface_net_thermal_radiation_J_m2 = surfthermal.variables['str'][-2:, :, :]
    
    # var name change
    conv_surface_net_solar_radiation_J_m2 = surface_net_solar_radiation_J_m2
    conv_surface_net_thermal_radiation_J_m2 = surface_net_thermal_radiation_J_m2

    # change data type to float32
    conv_surface_net_solar_radiation_J_m2 = change_dtype(conv_surface_net_solar_radiation_J_m2, 'float32')
    conv_surface_net_thermal_radiation_J_m2 = change_dtype(conv_surface_net_thermal_radiation_J_m2, 'float32')
    
    del surface_net_solar_radiation_J_m2
    del surface_net_thermal_radiation_J_m2
    
    return conv_surface_net_solar_radiation_J_m2,conv_surface_net_thermal_radiation_J_m2


def data_array(datapath,year,latitude,longitude,pre_surface_net_solar_radiation_J_m2, pre_surface_net_thermal_radiation_J_m2, \
    surfsolar,surfthermal,tmean,tdew,surfpres,uwnd,vwnd,erapet): 

    # loop through the days 168 at a time to overcome memory issues
    if year%4 == 0:
        dlen = 8784 # data length leap year
    else:
        dlen = 8760 # data length non leap year
    
    # create h5py file to write the eto values
    hdf5_store = h5py.File(datapath+str(year)+'_cache.hdf5', "a")
    eto_val = hdf5_store.create_dataset("eto_val", (dlen,pre_surface_net_solar_radiation_J_m2.shape[1],pre_surface_net_solar_radiation_J_m2.shape[2]), compression="gzip")

    for i in range(0,dlen,168):
        surface_net_solar_radiation_J_m2 = surfsolar.variables['ssr'][i:i+168, :, :]  # (time,latitude,longitude)
        surface_net_thermal_radiation_J_m2 = surfthermal.variables['str'][i:i+168, :, :]
        temperature2m_K = tmean.variables['t2m'][i:i+168, :, :]
        dewpoint2m_K = tdew.variables['d2m'][i:i+168, :, :]
        surface_pressure_Pa = surfpres.variables['sp'][i:i+168, :, :]
        u10m_m_s = uwnd.variables['u10'][i:i+168, :, :]
        v10m_m_s = vwnd.variables['v10'][i:i+168, :, :]

        # change data type to float32
        surface_net_solar_radiation_J_m2 = change_dtype(surface_net_solar_radiation_J_m2,'float32')
        surface_net_thermal_radiation_J_m2 = change_dtype(surface_net_thermal_radiation_J_m2,'float32')
        temperature2m_K = change_dtype(temperature2m_K,'float32')
        dewpoint2m_K = change_dtype(dewpoint2m_K,'float32')
        surface_pressure_Pa = change_dtype(surface_pressure_Pa,'float32')
        u10m_m_s = change_dtype(u10m_m_s,'float32')
        v10m_m_s = change_dtype(v10m_m_s,'float32')
   
        print(i)
        if i==0:
            # net radiation conversion to hourly values from cummulative
            nsr=np.concatenate((pre_surface_net_solar_radiation_J_m2,surface_net_solar_radiation_J_m2[:-1,:,:]),axis=0)
            nsr=np.delete(nsr,0,axis=0)
            netsolarrad = surface_net_solar_radiation_J_m2 - nsr
            # list of index this will be used for the other two variables too (solar and thermal radiation)
            index=np.arange(0,netsolarrad.shape[0])
            ind=np.where(index%24 == 0)
            ind=np.array(ind) + 1 
            
            netsolarrad[ind,:,:]=surface_net_solar_radiation_J_m2[ind,:,:]

            # net thermal radiation
            ntr=np.concatenate((pre_surface_net_thermal_radiation_J_m2, surface_net_thermal_radiation_J_m2[:-1,:,:]),axis=0)
            ntr=np.delete(ntr,0,axis=0)
            netthermalrad_v = surface_net_thermal_radiation_J_m2 - ntr
            netthermalrad_v[ind,:,:] = surface_net_thermal_radiation_J_m2[ind,:,:] # day change is set to the actual value

            # var name change
            conv_surface_net_solar_radiation_J_m2 = netsolarrad
            conv_surface_net_thermal_radiation_J_m2 = netthermalrad_v
            
            # delete variables not required
            del netsolarrad
            del netthermalrad_v
            del nsr
            del ntr
            
        else:
            # previous week data 
            last_surface_net_solar_radiation_J_m2 = surfsolar.variables['ssr'][i-168:i, :, :]
            last_surface_net_thermal_radiation_J_m2 = surfthermal.variables['str'][i-168:i, :, :]

            surf_netsolarrad_last = last_surface_net_solar_radiation_J_m2[-2:,:,:]
            surf_netthermalrad_last = last_surface_net_thermal_radiation_J_m2[-2:,:,:] 

            surf_netsolarrad_last = change_dtype(surf_netsolarrad_last,'float32')
            surf_netthermalrad_last = change_dtype(surf_netthermalrad_last,'float32')

            del last_surface_net_solar_radiation_J_m2
            del last_surface_net_thermal_radiation_J_m2


            # net radiation conversion to hourly values from cummulative
            nsr=np.concatenate((surf_netsolarrad_last,surface_net_solar_radiation_J_m2[:-1,:,:]),axis=0)
            nsr=np.delete(nsr,0,axis=0)
            netsolarrad = surface_net_solar_radiation_J_m2 - nsr
            # list of index this will be used for the other two variables too (solar and thermal radiation)
            index=np.arange(0,netsolarrad.shape[0])
            ind=np.where(index%24 == 0)
            ind=np.array(ind) + 1 
            
            netsolarrad[ind,:,:]=surface_net_solar_radiation_J_m2[ind,:,:]

            # net thermal radiation
            ntr=np.concatenate((surf_netthermalrad_last, surface_net_thermal_radiation_J_m2[:-1,:,:]),axis=0)
            ntr=np.delete(ntr,0,axis=0)
            netthermalrad_v = surface_net_thermal_radiation_J_m2 - ntr
            netthermalrad_v[ind,:,:] = surface_net_thermal_radiation_J_m2[ind,:,:] # day change is set to the actual value

            # var name change
            conv_surface_net_solar_radiation_J_m2 = netsolarrad
            conv_surface_net_thermal_radiation_J_m2 = netthermalrad_v
            
            # delete variables not required
            del netsolarrad
            del netthermalrad_v
            del nsr
            del ntr

        # 2. run unit_conv() to get all the write variables in right units
        surface_pressure_KPa, temperature2m_C, dewpoint2m_C, net_radiation_MJ_m2, windspeed2m_m_s,soil_hf = unit_conv(temperature2m_K,
                                                                                                          dewpoint2m_K, u10m_m_s,
                                                                                                          v10m_m_s,
                                                                                                          conv_surface_net_solar_radiation_J_m2,
                                                                                                          conv_surface_net_thermal_radiation_J_m2,
                                                                                                          surface_pressure_Pa)
        print('unit_conv done!')
        del temperature2m_K
        del dewpoint2m_K 
        del u10m_m_s
        del v10m_m_s
        del conv_surface_net_solar_radiation_J_m2
        del conv_surface_net_thermal_radiation_J_m2
        del surface_pressure_Pa
        del surface_net_solar_radiation_J_m2
        del surface_net_thermal_radiation_J_m2
        
        # hourly PET
        # 3. Run calc_pet() with pet_time='hourly'
        ET0_mm_hr = calculate_pet(surface_pressure_KPa,  # Surface pressure KPa
                              temperature2m_C,           # Daily mean temperature at 2 m
                              dewpoint2m_C,              # Daily mean dewpoint temperature at 2 m
                              windspeed2m_m_s,           # Windspeed at 2 m
                              net_radiation_MJ_m2,       # Total daily net downward radiation MJ/m2/day
                              soil_hf,                   # factor used to get the hourly soil heat flux from net radiation
                              'hourly')

        
        eto_val[i:i+168,:,:]=ET0_mm_hr
        #eto_val = np.append(eto_val, ET0_mm_hr[:,:,:], axis=0)
        print('eto calc done!')
        del ET0_mm_hr
        # delete all arrays to make space
        del surface_pressure_KPa
        del temperature2m_C
        del dewpoint2m_C
        del windspeed2m_m_s
        del net_radiation_MJ_m2
        del soil_hf
    
    # delete all arrays to make space
    del surfsolar
    del surfthermal
    del tmean
    del tdew
    del surfpres
    del uwnd
    del vwnd 
    del erapet
    
    return eto_val


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


def change_dtype(var, dtype):
    """
    This function change array data type to the given  datatype

    :param: var: the variable 
    :param dtype: the new data type
    """
    var = var.astype(dtype)
    
    return var


def unit_conv(block_temperature2m_K, block_dewpoint2m_K, block_u10m_m_s,
              block_v10m_m_s,block_surface_net_solar_radiation_J_m2,
              block_surface_net_thermal_radiation_J_m2,surface_pressure_Pa):
    """
    This function is used to convert the units of the input variables to approprate units.
    """
    
    # Convert units, etc of hourly data.
    row_temperature2m_C = block_temperature2m_K - 273.15
    row_dewpoint2m_C = block_dewpoint2m_K - 273.15
    del block_temperature2m_K, block_dewpoint2m_K

    # Wind speed at 2 m (use wind profile to scale from 10 m)
    row_u10m_m_s = block_u10m_m_s
    row_v10m_m_s = block_v10m_m_s
    temp_windspeed10m_m_s = np.sqrt(row_u10m_m_s**2 + row_v10m_m_s**2)
    row_windspeed2m_m_s = temp_windspeed10m_m_s*(4.87/(np.log(67.8*10-5.42)))
    del block_u10m_m_s, block_v10m_m_s, temp_windspeed10m_m_s

    # Net downward radiation.
    row_surface_net_solar_radiation_J_m2 = block_surface_net_solar_radiation_J_m2
    row_surface_net_thermal_radiation_J_m2 = block_surface_net_thermal_radiation_J_m2
    row_net_radiation_MJ_m2 = (row_surface_net_solar_radiation_J_m2 +
                               row_surface_net_thermal_radiation_J_m2) / 1e6
    del block_surface_net_solar_radiation_J_m2, block_surface_net_thermal_radiation_J_m2
    
    # Surface pressure Pa to KPa.
    surface_pressure_KPa = surface_pressure_Pa / 1000.
    del surface_pressure_Pa

    # soil heat flux
    # change the dtype to 'float16' 
    row_surface_net_solar_radiation_J_m2 = row_surface_net_solar_radiation_J_m2.astype('float16') 
    soil_hf = np.copy(row_net_radiation_MJ_m2)
    # soil heat flux (condition, day, night) 
    soil_hf = np.where(row_surface_net_solar_radiation_J_m2 > 0.0, soil_hf*0.1, soil_hf*0.5)
      
    del row_surface_net_solar_radiation_J_m2
  
    return surface_pressure_KPa,row_temperature2m_C,row_dewpoint2m_C,row_net_radiation_MJ_m2,row_windspeed2m_m_s,soil_hf


def calculate_pet(surface_pressure_KPa,     # surface pressure KPa
                  temperature2m_C,          # Daily mean temperature at 2 m
                  dewpoint2m_C,             # Daily mean dewpoint temperature at 2 m
                  windspeed2m_m_s,          # Windspeed at 2 m
                  net_radiation_MJ_m2,      # Total daily net downward radiation MJ/m2/day
                  soil_hf,                  # factor used to get the soil heat flux
                  pet_time):                # 'daily' or 'hourly'  ETo value
    """
    This is the function that calculate the PET based on the PM method.
    """
    # Constants.
    lmbda = 2.45  # Latent heat of vaporization [MJ kg -1] (simplification in the FAO PenMon (latent heat of about 20°C)
    cp = 1.013e-3 # Specific heat at constant pressure [MJ kg-1 °C-1]
    eps = 0.622   # Ratio molecular weight of water vapour/dry air

    # Soil heat flux density [MJ m-2 day-1] - set to 0 following eq 42 in FAO
    G = soil_hf    
   
    # Atmospheric pressure [kPa] eq 7 in FAO.
    P_kPa = surface_pressure_KPa #101.3*((293.0-0.0065*height_m) / 293.0)**5.26

    # Psychrometric constant (gamma symbol in FAO) eq 8 in FAO.
    psychometric_kPa_c = cp*P_kPa / (eps*lmbda)

    # Saturation vapour pressure, eq 11 in FAO.
    svp_kPa = 0.6108*np.exp((17.27*temperature2m_C) / (temperature2m_C+237.3))

    # Delta (slope of saturation vapour pressure curve) eq 13 in FAO.
    delta_kPa_C = 4098.0*svp_kPa / (temperature2m_C+237.3)**2

    # Actual vapour pressure, eq 14 in FAO.
    avp_kPa = 0.6108*np.exp((17.27*dewpoint2m_C) / (dewpoint2m_C+237.3))

    # Saturation vapour pressure deficit.
    svpdeficit_kPa = svp_kPa - avp_kPa


    if pet_time == 'daily':
        # Calculate ET0, equation 6 in FAO
        numerator = 0.408*delta_kPa_C*(net_radiation_MJ_m2 - G) + \
            psychometric_kPa_c*(900/(temperature2m_C+273))*windspeed2m_m_s*svpdeficit_kPa
        denominator = delta_kPa_C + psychometric_kPa_c*(1 + 0.34*windspeed2m_m_s)
    
        ET0_mm_day = numerator / denominator
        return ET0_mm_day
    
    elif pet_time == 'hourly':
        # Calculate ET0, equation 53 in FAO
        numerator = 0.408*delta_kPa_C*(net_radiation_MJ_m2 - G) + \
            psychometric_kPa_c*(37/(temperature2m_C+273))*windspeed2m_m_s*svpdeficit_kPa
        denominator = delta_kPa_C + psychometric_kPa_c*(1 + 0.34*windspeed2m_m_s)
    
        ET0_mm_hr = numerator / denominator
        return ET0_mm_hr
    
    else:
        raise ValueError("time only takes 'daily' or 'hourly'")

# ------------------------------------------------------------#
if __name__=='__main__':
    start=dt.datetime.now()
    
    wrapper(int(sys.argv[1]))
    
    end=dt.datetime.now()
    diff=end-start
    print('Time it took to run: %s'%diff)
