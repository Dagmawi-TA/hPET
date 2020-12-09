

**************************************************************************************
Script name         	Description
**************************************************************************************
pet_valc_v3.3.py    	This is the script used to calculate the PET values
		    	using Penman Monteith method. 
                    	The versions pet_valc_v3.3_1981.py and pet_valc_v3.3_2019.py
                    	with year values are just a variation 
                    	of the the script for the first and last year.

hdf2netcdf_v2.py        This is the script used to change the .hdf5 files 
                    	of the PET values to a netCDF file (.nc). This is 
                    	done in two steps where the script produce two
                    	netcdf files dividing the time array in to two. 
 		    	Then a cdo tool is used to merge the netcdf files
                    	to get one netcdf file for each year.

hourly2daily_pet_v3.py  This script is used to convert the hourly PET values
                        to daily values. The daily values account for the 
                        time offset at each grid location. Similarly this 
                        conversion is done in two steps where the netcdf 
                        files produced are merged later using cdo tools.
                        The versions hourly2daily_pet_v3_1981.py and 
                    	hourly2daily_pet_v3_2019.py
                    	with year values are just a variation 
                    	of the the script for the first and last year.
download_hPET.py       This script contain code to extract selected regions from the 
                       globe to download hPET and dPET dataset fro the server. 
		       The user need to download the script and change the locations before
		       running the code to get the data.
*************************************************************************************                 
