
###########################################################################
#
# PYTHON 3 FUNCTIONS FOR read_datasets
#
# Built with the following conda environment:
#
# name: bhoover-obs_match_3d
# channels:
#   - conda-forge
#   - defaults
# dependencies:
#   - python=3
#   - numpy
#   - pandas
#   - pynio
#   - matplotlib
#   - cartopy
#   - jupyter
#   - netCDF4
#   - scikit-learn
#   - dask
#   - geopy
#   - pip
#   - pip:
#
###########################################################################
#
# Import modules
#
from os.path import exists
import sys
import numpy as np #....................................................... Array module
#import dask.array as da #.................................................. Dask module
#from dask.distributed import Client, LocalCluster #........................ Dask client modules (LocalCluster == runs on a local machine)
import datetime as dt #.................................................... Datetime module
#from geopy.distance import geodesic #...................................... Geodesic distance module
import time #.............................................................. Time module
from netCDF4 import Dataset #.............................................. netCDF module

from quality_controls import qc_aeolus
#from quality_controls import qc_aircraft
from quality_controls import qc_amv
#
###########################################################################

# Functions for reading input datasets
#	Read in all variables.

# -------------------------------------------------------------------------
# Check if pressure or height variable exists in NetCDF index file
#
#	INPUTS:
#		dset_file ........................... Index filename
#		Pvar ................................ Pressure differences variable name
#		Hvar ................................ Height differences variable name
#
#	OUTPUTS:
#		PHmatch ............................. Array of either pressure differences or height differences
#		PHstr ............................... String indicating the variable that exists (i.e., Pressure or Height)
#
def check_PHvar_exists(dset_file,Pvar,Hvar):

	# load dataset
  data_hdl = Dataset(dset_file)

  print("Pvar = "+str(Pvar))
  print("Hvar = "+str(Hvar))
  try:
  	# check if pressure diff array exists
    PHmatch	= np.asarray( data_hdl.variables[Pvar] )
    PHstr	= "Pressure"
  except:	
  	# if pressure diffs does NOT exist
    try:	
    	# check if height diffs array exists
      PHmatch	= np.asarray( data_hdl.variables[Hvar] )
      PHstr	= "Height"
    except: 
      print("ERROR in check_PHvar_exists: DP_match and HT_match do not exist! At least one must exist to continue.")
      sys.exit()

  data_hdl.close()

	# Return variable
  return PHmatch,PHstr

# -------------------------------------------------------------------------
# Read AEOLUS (for collocation)
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		dateB4 .............................. Day before "yyyymmddhh" in yyyymmdd format
#		driver_dset_type .................... Aeolus dataset type: orig=original (not reprocessed), for reprocessed files use baseline number (example: B10)
#		driver_wind_type .................... Aeolus wind data type: RayClear=Rayleigh-clear, MieCloud=Mie-cloudy
#		bool_drv_qc ......................... Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		indexesD ............................ Indices of input obs
#		qc_list ............................. List of QC applied (if applicable)
#
def read_aeolus(path_prefix,yyyymmddhh,dateB4,driver_dset_type,driver_wind_type,bool_drv_qc):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]
  hour = yyyymmddhh[8:10]
  
  yyyymmdd = yyyy+mm+dd
  
  yyB4 = dateB4[0:4]
  mmB4 = dateB4[4:6]
  ddB4 = dateB4[6:8]

	#-------------------------------------
        # Find hour limits ... for Aeolus datafiles

  if hour == "00":
    hB4 = "21"
    hA  = "03"
  elif hour == "06":
    hB4 = "03"
    hA  = "09"
  elif hour == "12":
    hB4 = "09"
    hA  = "15"
  elif hour == "18":
    hB4 = "15"
    hA  = "21"

	#-------------------------------------------------
    	# Define dataset

  if driver_dset_type=='orig':						# Aeolus original dataset (not reprocessed)
    driver_dset_type_str = 'original'
  elif driver_dset_type!='orig':					# Aeolus dataset reprocessed by ESA
    driver_dset_type_str = 'reprocessed/2'+str(driver_dset_type)	# 	Reprocessed with Baseline B10 processor

  tmp_driver_path = '/aeolus-dataset/netcdf/'+driver_dset_type_str+'/'+yyyy+'/'+mm+'/'
  tmp_driver_path_B4 = '/aeolus-dataset/netcdf/'+driver_dset_type_str+'/'+yyB4+'/'+mmB4+'/'

		# Paths/files
  driver_path     	= path_prefix+tmp_driver_path
  driver_path_B4  	= path_prefix+tmp_driver_path_B4

  driver_filename 	= 'Aeolus.L2B.'+driver_wind_type+'.1day.'+yyyymmdd+'.nc'
  driver_filename_B4 	= 'Aeolus.L2B.'+driver_wind_type+'.1day.'+dateB4+'.nc'

  driver_file     	= driver_path+driver_filename
  driver_file_B4  	= driver_path_B4+driver_filename_B4

  driver_exists = exists(driver_file)
  if driver_exists==False:
    print("ERROR: file "+driver_file+" does not exist!")
    sys.exit()
  driverB4_exists = exists(driver_file_B4)
  if hour=="00" and driverB4_exists==False:
    print("ERROR: file "+driver_file_B4+" does not exist!")
    sys.exit()
  print("AEOLUS file: "+str(driver_file))

		# Paths on FTP/web archive server (for output NetCDF only)
  str_driver_path 	= driver_path
  str_driver_path_B4 	= driver_path_B4

  if hour=="00":
    drv_src = str_driver_path_B4+driver_filename_B4+", "+driver_path+driver_filename
  else:
    drv_src = str_driver_path+driver_filename

		# Variable names
  drv_lat_var 		= 'latitude'
  drv_lon_var 		= 'longitude'
  drv_prs_var 		= 'pressure'
  drv_hgt_var 		= 'height_mid'
  drv_yr_var  		= 'year'
  drv_mm_var  		= 'month'
  drv_dy_var 		= 'day'
  drv_hr_var  		= 'hour'
  drv_mn_var  		= 'minute'

  drv_err_var 		= 'HLOS_error'
  drv_len_var 		= 'length'
  drv_hgtTop_var 	= 'height_top'
  drv_hgtBot_var 	= 'height_bot'

    	#-------------------------------------------------
  	# Load dataset
  
  	#`````````````````````````````````````
	# CURRENT DATE
  data_hdl = Dataset(driver_file)

  tdrvC_lat = np.asarray( data_hdl.variables[drv_lat_var] )
  tdrvC_lon = np.asarray( data_hdl.variables[drv_lon_var] )
  tdrvC_prs = np.asarray( data_hdl.variables[drv_prs_var] )
  tdrvC_hgt = np.asarray( data_hdl.variables[drv_hgt_var] )
  tdrvC_yr  = np.asarray( data_hdl.variables[drv_yr_var]  )
  tdrvC_mm  = np.asarray( data_hdl.variables[drv_mm_var]  )
  tdrvC_dy  = np.asarray( data_hdl.variables[drv_dy_var]  )
  tdrvC_hr  = np.asarray( data_hdl.variables[drv_hr_var]  )
  tdrvC_mn  = np.asarray( data_hdl.variables[drv_mn_var]  )

  tdrvC_err     = np.asarray( data_hdl.variables[drv_err_var]  )
  tdrvC_len     = np.asarray( data_hdl.variables[drv_len_var]  )
  tdrvC_hgtTop  = np.asarray( data_hdl.variables[drv_hgtTop_var]  )
  tdrvC_hgtBot  = np.asarray( data_hdl.variables[drv_hgtBot_var]  )

  data_hdl.close() 

	# check pressure units and convert to hPa
  if max(tdrvC_prs) > 10000.:
    tdrvC_prs = tdrvC_prs/100.

	# check height units and convert to km
  if max(tdrvC_hgt) > 1000.:
    tdrvC_hgt = tdrvC_hgt/1000.

	# check height Top units and convert to km
  if max(tdrvC_hgtTop) > 1000.:
    tdrvC_hgtTop = tdrvC_hgtTop/1000.

	# check height Bottom units and convert to km
  if max(tdrvC_hgtBot) > 1000.:
    tdrvC_hgtBot = tdrvC_hgtBot/1000.

  
  if hour == "00":  
  	#`````````````````````````````````````
	# DATE BEFORE CURRENT
    data_hdl = Dataset(driver_file_B4)

    tdrvB4_lat = np.asarray( data_hdl.variables[drv_lat_var] )
    tdrvB4_lon = np.asarray( data_hdl.variables[drv_lon_var] )
    tdrvB4_prs = np.asarray( data_hdl.variables[drv_prs_var] )
    tdrvB4_hgt = np.asarray( data_hdl.variables[drv_hgt_var] )
    tdrvB4_yr  = np.asarray( data_hdl.variables[drv_yr_var]  )
    tdrvB4_mm  = np.asarray( data_hdl.variables[drv_mm_var]  )
    tdrvB4_dy  = np.asarray( data_hdl.variables[drv_dy_var]  )
    tdrvB4_hr  = np.asarray( data_hdl.variables[drv_hr_var]  )
    tdrvB4_mn  = np.asarray( data_hdl.variables[drv_mn_var]  )

    tdrvB4_err	= np.asarray( data_hdl.variables[drv_err_var]  )
    tdrvB4_len	= np.asarray( data_hdl.variables[drv_len_var]  )
    tdrvB4_hgtTop = np.asarray( data_hdl.variables[drv_hgtTop_var]  )
    tdrvB4_hgtBot = np.asarray( data_hdl.variables[drv_hgtBot_var]  )

    data_hdl.close()

  	# check pressure units and convert to hPa
    if max(tdrvB4_prs) > 10000.:
      tdrvB4_prs = tdrvB4_prs/100.

  	# check height units and convert to km
    if max(tdrvB4_hgt) > 1000.:
      tdrvB4_hgt = tdrvB4_hgt/1000.

  	# check height Top units and convert to km
    if max(tdrvB4_hgtTop) > 1000.:
      tdrvB4_hgtTop = tdrvB4_hgtTop/1000.
  
  	# check height Bottom units and convert to km
    if max(tdrvB4_hgtBot) > 1000.:
      tdrvB4_hgtBot = tdrvB4_hgtBot/1000.

	#`````````````````````````````````````
	# Extract 6-hr range from Aeolus day arrays: +/- 3 hours around center analysis hour of current date (00, 06, 12, or 18 UTC)

    i_tdrvC_hr  = tdrvC_hr.astype("int")
    i_tdrvB4_hr = tdrvB4_hr.astype("int")
  
    tiHHC = np.where(((i_tdrvC_hr >= 0) * (i_tdrvC_hr < int(hA))))
    siHHC = np.asarray(tiHHC)
    iHHC  = siHHC.flatten()

    tiHHB4 = np.where(((i_tdrvB4_hr >= 21) * (i_tdrvB4_hr < 24)))  
    siHHB4 = np.asarray(tiHHB4)
    iHHB4  = siHHB4.flatten()

    hdrvC_lat = tdrvC_lat[iHHC]
    hdrvC_lon = tdrvC_lon[iHHC]
    hdrvC_prs = tdrvC_prs[iHHC]
    hdrvC_hgt = tdrvC_hgt[iHHC]
    hdrvC_yr  = tdrvC_yr [iHHC]
    hdrvC_mm  = tdrvC_mm [iHHC]
    hdrvC_dy  = tdrvC_dy [iHHC]
    hdrvC_hr  = tdrvC_hr [iHHC]
    hdrvC_mn  = tdrvC_mn [iHHC]
    hdrvC_err     = tdrvC_err [iHHC]
    hdrvC_len     = tdrvC_len [iHHC]
    hdrvC_hgtTop  = tdrvC_hgtTop [iHHC]
    hdrvC_hgtBot  = tdrvC_hgtBot [iHHC]

    hdrvB4_lat = tdrvB4_lat[iHHB4]
    hdrvB4_lon = tdrvB4_lon[iHHB4]
    hdrvB4_prs = tdrvB4_prs[iHHB4]
    hdrvB4_hgt = tdrvB4_hgt[iHHB4]
    hdrvB4_yr  = tdrvB4_yr [iHHB4]
    hdrvB4_mm  = tdrvB4_mm [iHHB4]
    hdrvB4_dy  = tdrvB4_dy [iHHB4]
    hdrvB4_hr  = tdrvB4_hr [iHHB4]
    hdrvB4_mn  = tdrvB4_mn [iHHB4]
    hdrvB4_err     = tdrvB4_err [iHHB4]
    hdrvB4_len     = tdrvB4_len [iHHB4]
    hdrvB4_hgtTop  = tdrvB4_hgtTop [iHHB4]
    hdrvB4_hgtBot  = tdrvB4_hgtBot [iHHB4]

	# Append current arrays to before date arrays
    tdrv_lat = np.append(hdrvB4_lat,hdrvC_lat,axis=0)
    tdrv_lon = np.append(hdrvB4_lon,hdrvC_lon,axis=0)
    tdrv_prs = np.append(hdrvB4_prs,hdrvC_prs,axis=0)
    tdrv_hgt = np.append(hdrvB4_hgt,hdrvC_hgt,axis=0)
    tdrv_yr  = np.append(hdrvB4_yr ,hdrvC_yr ,axis=0)
    tdrv_mm  = np.append(hdrvB4_mm ,hdrvC_mm ,axis=0)
    tdrv_dy  = np.append(hdrvB4_dy ,hdrvC_dy ,axis=0)
    tdrv_hr  = np.append(hdrvB4_hr ,hdrvC_hr ,axis=0)
    tdrv_mn  = np.append(hdrvB4_mn ,hdrvC_mn ,axis=0)

    tdrv_err     = np.append(hdrvB4_err    ,hdrvC_err ,axis=0)
    tdrv_len     = np.append(hdrvB4_len    ,hdrvC_len ,axis=0)
    tdrv_hgtTop  = np.append(hdrvB4_hgtTop ,hdrvC_hgtTop ,axis=0)
    tdrv_hgtBot  = np.append(hdrvB4_hgtBot ,hdrvC_hgtBot ,axis=0)
    
   	# apply QC if bool_drv_qc=True
    if bool_drv_qc:
      #tindexesDC = []		#indexes of appended array (not original Aeolus indexes)
      #qc_list    = []
      tindexesDC,qc_list = qc_aeolus(driver_wind_type,tdrv_prs,tdrv_err,tdrv_len,tdrv_hgtTop,tdrv_hgtBot)
      
      sindexesDC = np.asarray(tindexesDC)
      	#NOTE: indexesDC includes both current date (C) and before date (B) indices. This combo is saved in 
	#      the collocation index files (for 00 UTC files only). The plotting code knows this and applies
	#      the same approach as here (find 6-hr window, then append B and C indices) to find indices of 
	#      collocated pairs.
      indexesDC  = sindexesDC.flatten()

      d_lat = tdrv_lat[indexesDC]
      d_lon = tdrv_lon[indexesDC]
      d_prs = tdrv_prs[indexesDC]
      d_hgt = tdrv_hgt[indexesDC]
      d_yr  = tdrv_yr[indexesDC]
      d_mm  = tdrv_mm[indexesDC]
      d_dy  = tdrv_dy[indexesDC]
      d_hr  = tdrv_hr[indexesDC]
      d_mn  = tdrv_mn[indexesDC]
      
      d_err = tdrv_err[indexesDC]
      d_len = tdrv_len[indexesDC]

    elif not bool_drv_qc:
      print("Do not apply AEOLUS QC")

      sindexesDC = np.asarray(np.where(tdrvC_lat==tdrvC_lat))           #get all indices
      indexesDC  = sindexesDC.flatten()

      qc_list = "No QC applied"

      d_lat = tdrv_lat 
      d_lon = tdrv_lon 
      d_prs = tdrv_prs 
      d_hgt = tdrv_hgt 
      d_yr  = tdrv_yr  
      d_mm  = tdrv_mm  
      d_dy  = tdrv_dy  
      d_hr  = tdrv_hr  
      d_mn  = tdrv_mn
      
      d_err = tdrv_err
      d_len = tdrv_len

    del tdrv_lat
    del tdrv_lon
    del tdrv_prs
    del tdrv_hgt
    del tdrv_yr 
    del tdrv_mm 
    del tdrv_dy 
    del tdrv_hr 
    del tdrv_mn
    del tdrv_err
    del tdrv_len
    del tdrv_hgtTop
    del tdrv_hgtBot

    indexesD     = indexesDC

  else:	# HOUR = 06, 12, 18

	# apply QC if bool_drv_qc=True
    if bool_drv_qc:
      #tindexesDC = []
      #qc_list    = []
      tindexesDC,qc_list = qc_aeolus(driver_wind_type,tdrvC_prs,tdrvC_err,tdrvC_len,tdrvC_hgtTop,tdrvC_hgtBot)
    
      sindexesDC = np.asarray(tindexesDC)
      indexesDC  = sindexesDC.flatten()

      tdrv_lat = tdrvC_lat[indexesDC]
      tdrv_lon = tdrvC_lon[indexesDC]
      tdrv_prs = tdrvC_prs[indexesDC]
      tdrv_hgt = tdrvC_hgt[indexesDC]
      tdrv_yr  = tdrvC_yr[indexesDC]
      tdrv_mm  = tdrvC_mm[indexesDC]
      tdrv_dy  = tdrvC_dy[indexesDC]
      tdrv_hr  = tdrvC_hr[indexesDC]
      tdrv_mn  = tdrvC_mn[indexesDC]
      
      tdrv_err = tdrvC_err [indexesDC]
      tdrv_len = tdrvC_len [indexesDC]

    elif not bool_drv_qc:
      print("Do not apply AEOLUS QC")

      sindexesDC = np.asarray(np.where(tdrvC_lat==tdrvC_lat))		#get all indices
      indexesDC  = sindexesDC.flatten()

      qc_list = "No QC applied"

      tdrv_lat = tdrvC_lat 
      tdrv_lon = tdrvC_lon 
      tdrv_prs = tdrvC_prs 
      tdrv_hgt = tdrvC_hgt 
      tdrv_yr  = tdrvC_yr  
      tdrv_mm  = tdrvC_mm  
      tdrv_dy  = tdrvC_dy  
      tdrv_hr  = tdrvC_hr  
      tdrv_mn  = tdrvC_mn
      
      tdrv_err = tdrvC_err	 
      tdrv_len = tdrvC_len	 

    del tdrvC_lat
    del tdrvC_lon
    del tdrvC_prs
    del tdrvC_hgt
    del tdrvC_yr 
    del tdrvC_mm 
    del tdrvC_dy 
    del tdrvC_hr 
    del tdrvC_mn
    
    del tdrvC_err
    del tdrvC_len

	#-------------------------------------------------
  	# Variables for collocation

    i_tdrvC_hr = tdrv_hr.astype("int")

    tiHHC = np.where((i_tdrvC_hr >= int(hB4)) * (i_tdrvC_hr < int(hA)))
    siHHC = np.asarray(tiHHC)
    iHHC  = siHHC.flatten()

    d_lat = tdrv_lat[iHHC]
    d_lon = tdrv_lon[iHHC]
    d_prs = tdrv_prs[iHHC]
    d_hgt = tdrv_hgt[iHHC]
    d_yr  = tdrv_yr [iHHC]
    d_mm  = tdrv_mm [iHHC]
    d_dy  = tdrv_dy [iHHC]
    d_hr  = tdrv_hr [iHHC]
    d_mn  = tdrv_mn [iHHC]
    
    d_err = tdrv_err[iHHC]
    d_len = tdrv_len[iHHC]

    indexesD = indexesDC[iHHC]
   
    print("indexesD max idx = "+str(max(indexesD)))
 
  	# Return variables to MAIN
  return d_lat,d_lon,d_prs,d_hgt,d_yr,d_mm,d_dy,d_hr,d_mn,indexesD,qc_list,drv_src, d_err,d_len

# -------------------------------------------------------------------------
# Read AEOLUS (for plotting)
#	Don't QC or separate into 6-hr chunks.
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		dateB4 .............................. Day before "yyyymmddhh" in yyyymmdd format
#		driver_dset_type .................... Aeolus dataset type: orig=original (not reprocessed), for reprocessed files use baseline number (example: B10)
#		driver_wind_type .................... Aeolus wind data type: RayClear=Rayleigh-clear, MieCloud=Mie-cloudy
#		bool_drv_qc ......................... Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		indexesD ............................ Indices of input obs
#		qc_list ............................. List of QC applied (if applicable)
#
def read_aeolus_for_plotting(path_prefix,yyyymmddhh,dateB4,driver_dset_type,driver_wind_type,bool_drv_qc):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]
  hour = yyyymmddhh[8:10]
  
  yyyymmdd = yyyy+mm+dd
  
  yyB4 = dateB4[0:4]
  mmB4 = dateB4[4:6]
  ddB4 = dateB4[6:8]

	#-------------------------------------
        # Find hour limits ... for Aeolus datafiles

  if hour == "00":
    hB4 = "21"
    hA  = "03"
  elif hour == "06":
    hB4 = "03"
    hA  = "09"
  elif hour == "12":
    hB4 = "09"
    hA  = "15"
  elif hour == "18":
    hB4 = "15"
    hA  = "21"

	#-------------------------------------------------
    	# Define dataset

  if driver_dset_type=='orig':						# Aeolus original dataset (not reprocessed)
    driver_dset_type_str = 'original'
  elif driver_dset_type!='orig':					# Aeolus dataset reprocessed by ESA
    driver_dset_type_str = 'reprocessed/2'+str(driver_dset_type)	# 	Reprocessed with Baseline B10 processor

  tmp_driver_path = '/aeolus-dataset/netcdf/'+driver_dset_type_str+'/'+yyyy+'/'+mm+'/'
  tmp_driver_path_B4 = '/aeolus-dataset/netcdf/'+driver_dset_type_str+'/'+yyB4+'/'+mmB4+'/'

		# Paths/files
  driver_path     	= path_prefix+tmp_driver_path
  driver_path_B4  	= path_prefix+tmp_driver_path_B4

  driver_filename 	= 'Aeolus.L2B.'+driver_wind_type+'.1day.'+yyyymmdd+'.nc'
  driver_filename_B4 	= 'Aeolus.L2B.'+driver_wind_type+'.1day.'+dateB4+'.nc'

  driver_file     	= driver_path+driver_filename
  driver_file_B4  	= driver_path_B4+driver_filename_B4

  driver_exists = exists(driver_file)
  if driver_exists==False:
    print("ERROR: file "+driver_file+" does not exist!")
    sys.exit()
  driverB4_exists = exists(driver_file_B4)
  if hour=="00" and driverB4_exists==False:
    print("ERROR: file "+driver_file_B4+" does not exist!")
    sys.exit()

		# Paths on FTP/web archive server (for output NetCDF only)
  str_driver_path 	= driver_path
  str_driver_path_B4 	= driver_path_B4

  if hour=="00":
    drv_src = str_driver_path_B4+driver_filename_B4+", "+driver_path+driver_filename
  else:
    drv_src = str_driver_path+driver_filename

		# Variable names
  drv_lat_var 		= 'latitude'
  drv_lon_var 		= 'longitude'
  drv_prs_var 		= 'pressure'
  drv_hgt_var 		= 'height_mid'
  drv_yr_var  		= 'year'
  drv_mm_var  		= 'month'
  drv_dy_var 		= 'day'
  drv_hr_var  		= 'hour'
  drv_mn_var  		= 'minute'

  drv_err_var 		= 'HLOS_error'
  drv_len_var 		= 'length'
  drv_hgtTop_var 	= 'height_top'
  drv_hgtBot_var 	= 'height_bot'
  
  drv_spd_var		= 'HLOS_wind_velocity'
  drv_dir_var		= 'HLOS_azimuth_angle'

    	#-------------------------------------------------
  	# Load dataset
  
  	#`````````````````````````````````````
	# CURRENT DATE
  data_hdl = Dataset(driver_file)

  tdrvC_lat = np.asarray( data_hdl.variables[drv_lat_var] )
  tdrvC_lon = np.asarray( data_hdl.variables[drv_lon_var] )
  tdrvC_prs = np.asarray( data_hdl.variables[drv_prs_var] )
  tdrvC_hgt = np.asarray( data_hdl.variables[drv_hgt_var] )
  tdrvC_yr  = np.asarray( data_hdl.variables[drv_yr_var]  )
  tdrvC_mm  = np.asarray( data_hdl.variables[drv_mm_var]  )
  tdrvC_dy  = np.asarray( data_hdl.variables[drv_dy_var]  )
  tdrvC_hr  = np.asarray( data_hdl.variables[drv_hr_var]  )
  tdrvC_mn  = np.asarray( data_hdl.variables[drv_mn_var]  )

  tdrvC_err     = np.asarray( data_hdl.variables[drv_err_var]  )
  tdrvC_len     = np.asarray( data_hdl.variables[drv_len_var]  )
  tdrvC_hgtTop  = np.asarray( data_hdl.variables[drv_hgtTop_var]  )
  tdrvC_hgtBot  = np.asarray( data_hdl.variables[drv_hgtBot_var]  )
  
  tdrvC_spd	= np.asarray( data_hdl.variables[drv_spd_var] )
  tdrvC_dir	= np.asarray( data_hdl.variables[drv_dir_var] )

  data_hdl.close() 

	# check pressure units and convert to hPa
  if max(tdrvC_prs) > 10000.:
    tdrvC_prs = tdrvC_prs/100.

	# check height units and convert to km
  if max(tdrvC_hgt) > 1000.:
    tdrvC_hgt = tdrvC_hgt/1000.

	# check height Top units and convert to km
  if max(tdrvC_hgtTop) > 1000.:
    tdrvC_hgtTop = tdrvC_hgtTop/1000.

	# check height Bottom units and convert to km
  if max(tdrvC_hgtBot) > 1000.:
    tdrvC_hgtBot = tdrvC_hgtBot/1000.

  
  if hour == "00":  
  	#`````````````````````````````````````
	# DATE BEFORE CURRENT
    data_hdl = Dataset(driver_file_B4)

    tdrvB4_lat = np.asarray( data_hdl.variables[drv_lat_var] )
    tdrvB4_lon = np.asarray( data_hdl.variables[drv_lon_var] )
    tdrvB4_prs = np.asarray( data_hdl.variables[drv_prs_var] )
    tdrvB4_hgt = np.asarray( data_hdl.variables[drv_hgt_var] )
    tdrvB4_yr  = np.asarray( data_hdl.variables[drv_yr_var]  )
    tdrvB4_mm  = np.asarray( data_hdl.variables[drv_mm_var]  )
    tdrvB4_dy  = np.asarray( data_hdl.variables[drv_dy_var]  )
    tdrvB4_hr  = np.asarray( data_hdl.variables[drv_hr_var]  )
    tdrvB4_mn  = np.asarray( data_hdl.variables[drv_mn_var]  )

    tdrvB4_err	= np.asarray( data_hdl.variables[drv_err_var]  )
    tdrvB4_len	= np.asarray( data_hdl.variables[drv_len_var]  )
    tdrvB4_hgtTop = np.asarray( data_hdl.variables[drv_hgtTop_var]  )
    tdrvB4_hgtBot = np.asarray( data_hdl.variables[drv_hgtBot_var]  )
    
    tdrvB4_spd	= np.asarray( data_hdl.variables[drv_spd_var] )
    tdrvB4_dir	= np.asarray( data_hdl.variables[drv_dir_var] )

    data_hdl.close()

  	# check pressure units and convert to hPa
    if max(tdrvB4_prs) > 10000.:
      tdrvB4_prs = tdrvB4_prs/100.

  	# check height units and convert to km
    if max(tdrvB4_hgt) > 1000.:
      tdrvB4_hgt = tdrvB4_hgt/1000.

  	# check height Top units and convert to km
    if max(tdrvB4_hgtTop) > 1000.:
      tdrvB4_hgtTop = tdrvB4_hgtTop/1000.
  
  	# check height Bottom units and convert to km
    if max(tdrvB4_hgtBot) > 1000.:
      tdrvB4_hgtBot = tdrvB4_hgtBot/1000.

	#`````````````````````````````````````
	# Extract 6-hr range from Aeolus day arrays: +/- 3 hours around center analysis hour of current date (00, 06, 12, or 18 UTC)

    i_tdrvC_hr  = tdrvC_hr.astype("int")
    i_tdrvB4_hr = tdrvB4_hr.astype("int")
  
    tiHHC = np.where(((i_tdrvC_hr >= 0) * (i_tdrvC_hr < int(hA))))
    siHHC = np.asarray(tiHHC)
    iHHC  = siHHC.flatten()

    tiHHB4 = np.where(((i_tdrvB4_hr >= 21) * (i_tdrvB4_hr < 24)))  
    siHHB4 = np.asarray(tiHHB4)
    iHHB4  = siHHB4.flatten()

    hdrvC_lat = tdrvC_lat[iHHC]
    hdrvC_lon = tdrvC_lon[iHHC]
    hdrvC_prs = tdrvC_prs[iHHC]
    hdrvC_hgt = tdrvC_hgt[iHHC]
    hdrvC_yr  = tdrvC_yr [iHHC]
    hdrvC_mm  = tdrvC_mm [iHHC]
    hdrvC_dy  = tdrvC_dy [iHHC]
    hdrvC_hr  = tdrvC_hr [iHHC]
    hdrvC_mn  = tdrvC_mn [iHHC]
    hdrvC_err     = tdrvC_err [iHHC]
    hdrvC_len     = tdrvC_len [iHHC]
    hdrvC_hgtTop  = tdrvC_hgtTop [iHHC]
    hdrvC_hgtBot  = tdrvC_hgtBot [iHHC]
    hdrvC_spd	= tdrvC_spd [iHHC]
    hdrvC_dir	= tdrvC_dir [iHHC]

    hdrvB4_lat = tdrvB4_lat[iHHB4]
    hdrvB4_lon = tdrvB4_lon[iHHB4]
    hdrvB4_prs = tdrvB4_prs[iHHB4]
    hdrvB4_hgt = tdrvB4_hgt[iHHB4]
    hdrvB4_yr  = tdrvB4_yr [iHHB4]
    hdrvB4_mm  = tdrvB4_mm [iHHB4]
    hdrvB4_dy  = tdrvB4_dy [iHHB4]
    hdrvB4_hr  = tdrvB4_hr [iHHB4]
    hdrvB4_mn  = tdrvB4_mn [iHHB4]
    hdrvB4_err     = tdrvB4_err [iHHB4]
    hdrvB4_len     = tdrvB4_len [iHHB4]
    hdrvB4_hgtTop  = tdrvB4_hgtTop [iHHB4]
    hdrvB4_hgtBot  = tdrvB4_hgtBot [iHHB4]
    hdrvB4_spd	= tdrvB4_spd [iHHB4]
    hdrvB4_dir	= tdrvB4_dir [iHHB4]

	# Append current arrays to before date arrays
    tdrv_lat = np.append(hdrvB4_lat,hdrvC_lat,axis=0)
    tdrv_lon = np.append(hdrvB4_lon,hdrvC_lon,axis=0)
    tdrv_prs = np.append(hdrvB4_prs,hdrvC_prs,axis=0)
    tdrv_hgt = np.append(hdrvB4_hgt,hdrvC_hgt,axis=0)
    tdrv_yr  = np.append(hdrvB4_yr ,hdrvC_yr ,axis=0)
    tdrv_mm  = np.append(hdrvB4_mm ,hdrvC_mm ,axis=0)
    tdrv_dy  = np.append(hdrvB4_dy ,hdrvC_dy ,axis=0)
    tdrv_hr  = np.append(hdrvB4_hr ,hdrvC_hr ,axis=0)
    tdrv_mn  = np.append(hdrvB4_mn ,hdrvC_mn ,axis=0)

    tdrv_err     = np.append(hdrvB4_err    ,hdrvC_err ,axis=0)
    tdrv_len     = np.append(hdrvB4_len    ,hdrvC_len ,axis=0)
    tdrv_hgtTop  = np.append(hdrvB4_hgtTop ,hdrvC_hgtTop ,axis=0)
    tdrv_hgtBot  = np.append(hdrvB4_hgtBot ,hdrvC_hgtBot ,axis=0)
    
    tdrv_spd	= np.append(hdrvB4_spd, hdrvC_spd, axis=0)
    tdrv_dir	= np.append(hdrvB4_dir, hdrvC_dir, axis=0)
    
   	# apply QC if bool_drv_qc=True
#    if bool_drv_qc:
#      #tindexesDC = []		#indexes of appended array (not original Aeolus indexes)
#      #qc_list    = []
#      tindexesDC,qc_list = qc_aeolus(driver_wind_type,tdrv_prs,tdrv_err,tdrv_len,tdrv_hgtTop,tdrv_hgtBot)
#      
#      sindexesDC = np.asarray(tindexesDC)
#      	#NOTE: indexesDC includes both current date (C) and before date (B) indices. This combo is saved in 
#	#      the collocation index files (for 00 UTC files only). The plotting code knows this and applies
#	#      the same approach as here (find 6-hr window, then append B and C indices) to find indices of 
#	#      collocated pairs.
#      indexesDC  = sindexesDC.flatten()
#
#      d_lat = tdrv_lat[indexesDC]
#      d_lon = tdrv_lon[indexesDC]
#      d_prs = tdrv_prs[indexesDC]
#      d_hgt = tdrv_hgt[indexesDC]
#      d_yr  = tdrv_yr[indexesDC]
#      d_mm  = tdrv_mm[indexesDC]
#      d_dy  = tdrv_dy[indexesDC]
#      d_hr  = tdrv_hr[indexesDC]
#      d_mn  = tdrv_mn[indexesDC]
#      
#      d_err = tdrv_err[indexesDC]
#      d_len = tdrv_len[indexesDC]
#      d_spd = tdrv_spd[indexesDC]
#      d_dir = tdrv_dir[indexesDC]
#
#    elif not bool_drv_qc:
#      print("Do not apply AEOLUS QC")

    sindexesDC = np.asarray(np.where(tdrv_lat==tdrv_lat))           #get all indices
    indexesDC  = sindexesDC.flatten()

#      qc_list = "No QC applied"

    d_lat = tdrv_lat 
    d_lon = tdrv_lon 
    d_prs = tdrv_prs 
    d_hgt = tdrv_hgt 
    d_yr  = tdrv_yr  
    d_mm  = tdrv_mm  
    d_dy  = tdrv_dy  
    d_hr  = tdrv_hr  
    d_mn  = tdrv_mn
    
    d_err = tdrv_err
    d_len = tdrv_len
    d_spd = tdrv_spd
    d_dir = tdrv_dir

    del tdrv_lat
    del tdrv_lon
    del tdrv_prs
    del tdrv_hgt
    del tdrv_yr 
    del tdrv_mm 
    del tdrv_dy 
    del tdrv_hr 
    del tdrv_mn
    del tdrv_err
    del tdrv_len
    del tdrv_hgtTop
    del tdrv_hgtBot
    del tdrv_spd
    del tdrv_dir

    indexesD     = indexesDC

  else:	# HOUR = 06, 12, 18

	# apply QC if bool_drv_qc=True
#    if bool_drv_qc:
#      #tindexesDC = []
#      #qc_list    = []
#      tindexesDC,qc_list = qc_aeolus(driver_wind_type,tdrvC_prs,tdrvC_err,tdrvC_len,tdrvC_hgtTop,tdrvC_hgtBot)
#    
#      sindexesDC = np.asarray(tindexesDC)
#      indexesDC  = sindexesDC.flatten()
#
#      tdrv_lat = tdrvC_lat[indexesDC]
#      tdrv_lon = tdrvC_lon[indexesDC]
#      tdrv_prs = tdrvC_prs[indexesDC]
#      tdrv_hgt = tdrvC_hgt[indexesDC]
#      tdrv_yr  = tdrvC_yr[indexesDC]
#      tdrv_mm  = tdrvC_mm[indexesDC]
#      tdrv_dy  = tdrvC_dy[indexesDC]
#      tdrv_hr  = tdrvC_hr[indexesDC]
#      tdrv_mn  = tdrvC_mn[indexesDC]
#      
#      tdrv_err = tdrvC_err [indexesDC]
#      tdrv_len = tdrvC_len [indexesDC]
#      tdrv_spd = tdrvC_spd [indexesDC]
#      tdrv_dir = tdrvC_dir[indexesDC]
#
#    elif not bool_drv_qc:
#      print("Do not apply AEOLUS QC")

    sindexesDC = np.asarray(np.where(tdrvC_lat==tdrvC_lat))		#get all indices
    indexesDC  = sindexesDC.flatten()

#      qc_list = "No QC applied"

    tdrv_lat = tdrvC_lat 
    tdrv_lon = tdrvC_lon 
    tdrv_prs = tdrvC_prs 
    tdrv_hgt = tdrvC_hgt 
    tdrv_yr  = tdrvC_yr  
    tdrv_mm  = tdrvC_mm  
    tdrv_dy  = tdrvC_dy  
    tdrv_hr  = tdrvC_hr  
    tdrv_mn  = tdrvC_mn
    
    tdrv_err = tdrvC_err       
    tdrv_len = tdrvC_len       
    tdrv_spd = tdrvC_spd 
    tdrv_dir = tdrvC_dir

    del tdrvC_lat
    del tdrvC_lon
    del tdrvC_prs
    del tdrvC_hgt
    del tdrvC_yr 
    del tdrvC_mm 
    del tdrvC_dy 
    del tdrvC_hr 
    del tdrvC_mn
    
    del tdrvC_err
    del tdrvC_len
    del tdrvC_spd
    del tdrvC_dir

	#-------------------------------------------------
  	# Variables for collocation

    d_lat = tdrv_lat
    d_lon = tdrv_lon
    d_prs = tdrv_prs
    d_hgt = tdrv_hgt
    d_yr  = tdrv_yr 
    d_mm  = tdrv_mm 
    d_dy  = tdrv_dy 
    d_hr  = tdrv_hr 
    d_mn  = tdrv_mn 
    
    d_err = tdrv_err
    d_len = tdrv_len
    d_spd = tdrv_spd
    d_dir = tdrv_dir	

    indexesD = indexesDC
   
    print("indexesD max idx = "+str(max(indexesD)))
 
  	# Return variables to MAIN
  return d_lat,d_lon,d_prs,d_hgt,d_yr,d_mm,d_dy,d_hr,d_mn,indexesD,qc_list,drv_src, d_err,d_len,d_spd,d_dir

# -------------------------------------------------------------------------
# Read AIRCRAFT (from NCEP)
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		bool_qc ......................... Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		qc_list ............................. List of QC applied (if applicable)
#
def read_aircraft(path_prefix,yyyymmddhh,bool_qc):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

	#-------------------------------------------------
    	# Define dataset

  tmp_dset1_path = '/atmos-nc-dataset/aircraft/'+yyyy+'/'+mm+'/'+dd+'/'

		# Path/file
  dset1_path   		= path_prefix+tmp_dset1_path

  dset1_filename 	= 'gdas.'+yyyymmddhh+'.aircft.tm00.bufr_d.nc4'

  dset1_file   		= dset1_path+dset1_filename

  dset1_exists = exists(dset1_file)
  if dset1_exists==False:
    print("ERROR: file "+dset1_file+" does not exist!")
    sys.exit()							#exit script immediately

		# Path on FTP/web archive server (for output NetCDF only)
  str_dset1_path = dset1_path

  dset1_src = str_dset1_path+dset1_filename

		# Variable names
  dset1_lat_var = 'latitude'
  dset1_lon_var = 'longitude'
  dset1_prs_var = 'pressure'
  dset1_yr_var  = 'year'
  dset1_mm_var  = 'month'
  dset1_dy_var  = 'day'
  dset1_hr_var  = 'hour'
  dset1_mn_var  = 'minutes'
  
  dset1_ht1_var = 'flight_level'
  dset1_ht2_var = 'height'
  
  dset1_spd_var = 'wind_speed'
  dset1_dir_var = 'wind_direction'
		
    	#-------------------------------------------------
  	# Load dataset
	#	AIRCRAFT data is divided into Groups
  
  data_hdl = Dataset(dset1_file)

  grps = list(data_hdl.groups)

	# populate full arrays with Group1 (grps(0))
  td_lat = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lat_var] )
  td_lon = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lon_var] )
  td_yr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_yr_var] )
  td_mm  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mm_var] )
  td_dy  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dy_var] )
  td_hr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var] )
  td_mn  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mn_var] )

  if grps[0]=='NC004006' or grps[0]=='NC004009':
    td_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht2_var]  )
  else:
    td_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht1_var]  )
  
  td_spd = np.asarray( data_hdl.groups[grps[0]].variables[dset1_spd_var]  )
  td_dir = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dir_var]  )

	# append data from remaining Groups
  for x in range(len(grps)-1):
#      print("grp exists = "+grps[x+1])
      tdset1_lat = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lat_var] )
      tdset1_lon = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lon_var] )
      tdset1_yr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_yr_var]  )
      tdset1_mm  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mm_var]  )
      tdset1_dy  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dy_var]  )
      tdset1_hr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
      tdset1_mn  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mn_var]  )
      tdset1_spd = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_spd_var]  )
      tdset1_dir = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dir_var]  )

      td_lat = np.append(td_lat,tdset1_lat,axis=0)
      td_lon = np.append(td_lon,tdset1_lon,axis=0)
      td_yr  = np.append(td_yr, tdset1_yr, axis=0)
      td_mm  = np.append(td_mm, tdset1_mm, axis=0)
      td_dy  = np.append(td_dy, tdset1_dy, axis=0)
      td_hr  = np.append(td_hr, tdset1_hr, axis=0)
      td_mn  = np.append(td_mn, tdset1_mn, axis=0)
      td_spd = np.append(td_spd, tdset1_spd, axis=0)
      td_dir = np.append(td_dir, tdset1_dir, axis=0)

      if grps[x+1]=='NC004006' or grps[x+1]=='NC004009':
        z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht2_var]  )
      else:
        z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht1_var]  )

      td_hgt = np.append(td_hgt,z,axis=0)

  data_hdl.close()

	# apply QC if bool_dset1_qc=True
  if bool_qc:
    print("Apply AIRCRAFT QC: TBD")
#    tindexes1,qc_list = qc_aircraft(td_hgt)
#    
#    sindexes1 = np.asarray(tindexes1)
#    indexes1  = sindexes1.flatten()
# 
#    d_lat = td_lat[indexes1]
#    d_lon = td_lon[indexes1]
#    d_yr  = td_yr [indexes1]
#    d_mm  = td_mm [indexes1]
#    d_dy  = td_dy [indexes1]
#    d_hr  = td_hr [indexes1]
#    d_mn  = td_mn [indexes1]
#    d_hgt = td_hgt[indexes1]
#    d_spd = td_spd[indexes1]
#    d_dir = td_dir[indexes1]
    
  elif not bool_qc:
    print("Do not apply AIRCRAFT QC")

    sindexesDC = np.asarray(np.where(td_lat==td_lat))           #get all indices
    indexes1   = sindexesDC.flatten()

    qc_list = "No QC applied"

    d_lat = td_lat
    d_lon = td_lon
    d_yr  = td_yr
    d_mm  = td_mm
    d_dy  = td_dy
    d_hr  = td_hr
    d_mn  = td_mn
    d_hgt = td_hgt
    d_spd = td_spd
    d_dir = td_dir

	# create pressure array but fill with missing -999.
	#	AIRCRAFT pressures not available.
  d_prs = np.nan * np.ones_like(d_lat)

	# check height units and convert to km
  if max(d_hgt) > 1000.:
    d_hgt = d_hgt/1000.

	# Return variables to MAIN
  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes1,qc_list,dset1_src, d_spd,d_dir
  
# -------------------------------------------------------------------------
# Read AIRCRAFT (from NCEP)
#	Don't QC
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		bool_qc ......................... Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		qc_list ............................. List of QC applied (if applicable)
#
def read_aircraft_for_plotting(path_prefix,yyyymmddhh,bool_qc):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

	#-------------------------------------------------
    	# Define dataset

  tmp_dset1_path = '/atmos-nc-dataset/aircraft/'+yyyy+'/'+mm+'/'+dd+'/'

		# Path/file
  dset1_path   		= path_prefix+tmp_dset1_path

  dset1_filename 	= 'gdas.'+yyyymmddhh+'.aircft.tm00.bufr_d.nc4'

  dset1_file   		= dset1_path+dset1_filename

  dset1_exists = exists(dset1_file)
  if dset1_exists==False:
    print("ERROR: file "+dset1_file+" does not exist!")
    sys.exit()							#exit script immediately

		# Path on FTP/web archive server (for output NetCDF only)
  str_dset1_path = dset1_path

  dset1_src = str_dset1_path+dset1_filename

		# Variable names
  dset1_lat_var = 'latitude'
  dset1_lon_var = 'longitude'
  dset1_prs_var = 'pressure'
  dset1_yr_var  = 'year'
  dset1_mm_var  = 'month'
  dset1_dy_var  = 'day'
  dset1_hr_var  = 'hour'
  dset1_mn_var  = 'minutes'
  
  dset1_ht1_var = 'flight_level'
  dset1_ht2_var = 'height'
  
  dset1_spd_var = 'wind_speed'
  dset1_dir_var = 'wind_direction'
		
    	#-------------------------------------------------
  	# Load dataset
	#	AIRCRAFT data is divided into Groups
  
  data_hdl = Dataset(dset1_file)

  grps = list(data_hdl.groups)

	# populate full arrays with Group1 (grps(0))
  td_lat = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lat_var] )
  td_lon = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lon_var] )
  td_yr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_yr_var] )
  td_mm  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mm_var] )
  td_dy  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dy_var] )
  td_hr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var] )
  td_mn  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mn_var] )

  if grps[0]=='NC004006' or grps[0]=='NC004009':
    td_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht2_var]  )
  else:
    td_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht1_var]  )
  
  td_spd = np.asarray( data_hdl.groups[grps[0]].variables[dset1_spd_var]  )
  td_dir = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dir_var]  )

	# append data from remaining Groups
  for x in range(len(grps)-1):
#      print("grp exists = "+grps[x+1])
      tdset1_lat = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lat_var] )
      tdset1_lon = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lon_var] )
      tdset1_yr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_yr_var]  )
      tdset1_mm  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mm_var]  )
      tdset1_dy  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dy_var]  )
      tdset1_hr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
      tdset1_mn  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mn_var]  )
      tdset1_spd = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_spd_var]  )
      tdset1_dir = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dir_var]  )

      td_lat = np.append(td_lat,tdset1_lat,axis=0)
      td_lon = np.append(td_lon,tdset1_lon,axis=0)
      td_yr  = np.append(td_yr, tdset1_yr, axis=0)
      td_mm  = np.append(td_mm, tdset1_mm, axis=0)
      td_dy  = np.append(td_dy, tdset1_dy, axis=0)
      td_hr  = np.append(td_hr, tdset1_hr, axis=0)
      td_mn  = np.append(td_mn, tdset1_mn, axis=0)
      td_spd = np.append(td_spd, tdset1_spd, axis=0)
      td_dir = np.append(td_dir, tdset1_dir, axis=0)

      if grps[x+1]=='NC004006' or grps[x+1]=='NC004009':
        z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht2_var]  )
      else:
        z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht1_var]  )

      td_hgt = np.append(td_hgt,z,axis=0)

  data_hdl.close()

  sindexesDC = np.asarray(np.where(td_lat==td_lat))	      #get all indices
  indexes1   = sindexesDC.flatten()

  qc_list = "No QC applied"

  d_lat = td_lat
  d_lon = td_lon
  d_yr  = td_yr 
  d_mm  = td_mm 
  d_dy  = td_dy 
  d_hr  = td_hr 
  d_mn  = td_mn 
  d_hgt = td_hgt
  d_spd = td_spd
  d_dir = td_dir

	# create pressure array but fill with missing -999.
	#	AIRCRAFT pressures not available.
  d_prs = np.nan * np.ones_like(d_lat)

	# check height units and convert to km
  if max(d_hgt) > 1000.:
    d_hgt = d_hgt/1000.

  print("min/max aircraft: "+str(min(d_spd))+","+str(max(d_spd)))

	# Return variables to MAIN
  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes1,qc_list,dset1_src, d_spd,d_dir

# -------------------------------------------------------------------------
# Read AMV (for 4th International Comparison)
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		pct ................................. Minimum AMV quality indicator (QI) in % for QC
#		bool_qc ............................. Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#		pct ................................. AMV QI in %
#               dset_type_str ....................... Clear or Cloudy AMV
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		d_pccf .............................. AMV QI (percent confidence) in %
#               indexes2 ............................ Indices of input obs
#		qc_list ............................. List of QC applied (if applicable)
#
def read_amv_4th_int(path_prefix,yyyymmddhh,bool_qc,pct,dset_type_str,shortname):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

  shortabbrv = shortname
  print("shortabbrev = "+shortabbrv)

	#-------------------------------------------------
    	# Define dataset

#  tmp_dset2_path = 'atmos-nc-dataset/AMV/'+yyyy+'/'+mm+'/'+dd+'/'
  tmp_dset2_path = '/data/users/daves/intercomparison2021/'

		# Path/file
  path_prefix		= ''		#reset for this dataset only
  dset2_path   		= path_prefix+tmp_dset2_path

  if dset_type_str == 'Cloud':
    if shortabbrv=='BRZ':
      dset2_filename    = 'BRZ_4th_AMVIC_INPE_Test_2a_final.Cloud.nc'
      t_name = 'BRZ_'+str(dset_type_str)
    elif shortabbrv=='EUM':
      dset2_filename    = 'EUM_AMVIntm_Chan14_20191020113000Z_GOES_ASCII_Test21.Cloud.nc'
      t_name = 'EUM_'+str(dset_type_str)
    elif shortabbrv=='JMA':
      dset2_filename    = 'JMA_test21_1130.Cloud.nc'
      t_name = 'JMA_'+str(dset_type_str)
    elif shortabbrv=='KMA':
      dset2_filename    = 'KMA_test21_goes16_abi_ch14_amv_201910201130.Cloud.nc'
      t_name = 'KMA_'+str(dset_type_str)
    elif shortabbrv=='NOA':
      dset2_filename    = 'NOA_ASCII_AMV-4thInt_TEST2.GOES16.2019293.1130.CH_14.FD.Cloud.nc'
      t_name = 'NOA_'+str(dset_type_str)
    elif shortabbrv=='NWC':
      dset2_filename    = 'NWC_INTERCOMP2021_NWCSAFTEST21.Cloud.nc'
      t_name = 'NWC_'+str(dset_type_str)
  elif dset_type_str == 'Clear':
    if shortabbrv=='BRZ':
      dset2_filename    = 'BRZ_4th_AMVIC_INPE_Test_3_final.Clear.nc'
      t_name = 'BRZ_'+str(dset_type_str)
    elif shortabbrv=='EUM':
      dset2_filename    = 'EUM_AMVIntm_Chan08_20191020113000Z_GOES_ASCII_Test31.Clear.nc'
      t_name = 'EUM_'+str(dset_type_str)
    elif shortabbrv=='KMA':
      dset2_filename    = 'KMA_test31_goes16_abi_ch08_amv_201910201130.Clear.nc'
      t_name = 'KMA_'+str(dset_type_str)
    elif shortabbrv=='NOA':
      dset2_filename    = 'NOA_ASCII_AMV-4thInt_TEST3.GOES16.2019293.1130.CH_08.FD.Clear.nc'
      t_name = 'NOA_'+str(dset_type_str)
    elif shortabbrv=='NWC':
      dset2_filename    = 'NWC_INTERCOMP2021_NWCSAFTEST31.Clear.nc'
      t_name = 'NWC_'+str(dset_type_str)
    elif shortabbrv=='BRZ4':
      dset2_filename    = 'BRZ_4th_AMVIC_INPE_Test_4_final.Clear.nc'
      t_name = 'BRZ4_'+str(dset_type_str)
    elif shortabbrv=='EUM4':
      dset2_filename    = 'EUM_AMVIntm_Chan08_20191020113000Z_GOES_ASCII_Test41.Clear.nc'
      t_name = 'EUM4_'+str(dset_type_str)
    elif shortabbrv=='KMA4':
      dset2_filename    = 'KMA_test41_goes16_abi_ch08_amv_201910201130.Clear.nc'
      t_name = 'KMA4_'+str(dset_type_str)
    elif shortabbrv=='NOA4':
      dset2_filename    = 'NOA_ASCII_AMV-4thInt_TEST4.GOES16.2019293.1130.CH_08.FD.Clear.nc'
      t_name = 'NOA4_'+str(dset_type_str)
    elif shortabbrv=='NWC4':
      dset2_filename    = 'NWC_INTERCOMP2021_NWCSAFTEST41.Clear.nc'
      t_name = 'NWC4_'+str(dset_type_str)

  dset2_file   		= dset2_path+dset2_filename

  dset2_exists = exists(dset2_file)
  if dset2_exists==False:
    print("ERROR: file "+dset2_file+" does not exist!")
    sys.exit()

		# Path (for output NetCDF only)
  str_dset2_path = dset2_path

  dset2_src = str_dset2_path+dset2_filename

		# Variable names
  dset2_lat_var  = 'Ob_Lat'
  dset2_lon_var  = 'Ob_Lon'
  dset2_prs_var  = 'Ob_Prs'		# units = hPa
  dset2_yr_var   = 'Ob_Year'
  dset2_mm_var   = 'Ob_Month'
  dset2_dy_var   = 'Ob_Day'
  dset2_hr_var   = 'Ob_Hour'
  dset2_mn_var   = 'Ob_Min'
  dset2_pccf_var = 'Ob_QI'
  
  dset2_spd_var	 = 'Ob_Spd'
  dset2_dir_var  = 'Ob_Dir'

    	#-------------------------------------------------
  	# Load dataset

  data_hdl = Dataset(dset2_file)

  tdset2_lat = np.asarray( data_hdl.variables[dset2_lat_var] )
  tdset2_lon = np.asarray( data_hdl.variables[dset2_lon_var] )		#NOTE: lons are [-180,180]. This will be converted to [0,360] automatically before collocation.
  tdset2_prs = np.asarray( data_hdl.variables[dset2_prs_var] )
  tdset2_yr  = np.asarray( data_hdl.variables[dset2_yr_var]  )
  tdset2_mm  = np.asarray( data_hdl.variables[dset2_mm_var]  )
  tdset2_dy  = np.asarray( data_hdl.variables[dset2_dy_var]  )
  tdset2_hr  = np.asarray( data_hdl.variables[dset2_hr_var]  )
  tdset2_mn  = np.asarray( data_hdl.variables[dset2_mn_var]  )
  tdset2_pccf = np.asarray( data_hdl.variables[dset2_pccf_var]  )

  tdset2_spd = np.asarray( data_hdl.variables[dset2_spd_var] )
  tdset2_dir = np.asarray( data_hdl.variables[dset2_dir_var] )

  data_hdl.close()

	# check pressure units and convert to hPa
  if max(tdset2_prs) > 10000.:
    tdset2_prs = tdset2_prs/100.

	# apply QC if bool_dset2_qc=True
  if bool_qc:
    
    tindexes2,qc_list = qc_amv(tdset2_pccf,pct)
    
    sindexes2 = np.asarray(tindexes2)
    indexes2  = sindexes2.flatten()
 
    d_lat = tdset2_lat[indexes2]
    d_lon = tdset2_lon[indexes2]
    d_prs = tdset2_prs[indexes2]
    d_yr  = tdset2_yr [indexes2]
    d_mm  = tdset2_mm [indexes2]
    d_dy  = tdset2_dy [indexes2]
    d_hr  = tdset2_hr [indexes2]
    d_mn  = tdset2_mn [indexes2]
    
    d_spd = tdset2_spd[indexes2]
    d_dir = tdset2_dir[indexes2]

  elif not bool_qc:
    print("Do not apply AMV_4th_Int QC")

    sindexesDC = np.asarray(np.where(tdset2_lat==tdset2_lat))           #get all indices
    indexes2   = sindexesDC.flatten()

    qc_list = "No QC applied"

    d_lat = tdset2_lat
    d_lon = tdset2_lon
    d_yr  = tdset2_yr 
    d_mm  = tdset2_mm 
    d_dy  = tdset2_dy 
    d_hr  = tdset2_hr 
    d_mn  = tdset2_mn  
    d_prs = tdset2_prs
    
    d_spd = tdset2_spd
    d_dir = tdset2_dir

	# create height array but fill with missing -999.
	#	AMV heights not available.
  d_hgt = np.nan * np.ones_like(d_lat)

  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes2,qc_list,dset2_src, d_spd,d_dir,t_name
  
# -------------------------------------------------------------------------
# Read AMV (for 4th International Comparison)
#	Don't QC
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		pct ................................. Minimum AMV quality indicator (QI) in % for QC
#		bool_qc ............................. Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		d_pccf .............................. AMV QI (percent confidence) in %
#               indexes2 ............................ Indices of input obs
#		qc_list ............................. List of QC applied (if applicable)
#
def read_amv_4th_int_for_plotting(path_prefix,yyyymmddhh,bool_qc,pct,dset_type_str,shortname):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

  shortabbrv = shortname[0:4]
  print("shortabbrev = "+shortabbrv)

	#-------------------------------------------------
    	# Define dataset

#  tmp_dset2_path = 'atmos-nc-dataset/AMV/'+yyyy+'/'+mm+'/'+dd+'/'
  tmp_dset2_path = '/data/users/daves/intercomparison2021/'

		# Path/file
  path_prefix		= ''		#reset for this dataset only
  dset2_path   		= path_prefix+tmp_dset2_path

  if dset_type_str == 'Cloud':
    if shortabbrv=='BRZ_':
      dset2_filename	= 'BRZ_4th_AMVIC_INPE_Test_2a_final.Cloud.nc'
      t_name = 'BRZ_'+str(dset_type_str)
    elif shortabbrv=='EUM_':
      dset2_filename	= 'EUM_AMVIntm_Chan14_20191020113000Z_GOES_ASCII_Test21.Cloud.nc'
      t_name = 'EUM_'+str(dset_type_str)
    elif shortabbrv=='JMA_':      
      dset2_filename	= 'JMA_test21_1130.Cloud.nc'
      t_name = 'JMA_'+str(dset_type_str)
    elif shortabbrv=='KMA_':
      dset2_filename	= 'KMA_test21_goes16_abi_ch14_amv_201910201130.Cloud.nc'
      t_name = 'KMA_'+str(dset_type_str)
    elif shortabbrv=='NOA_':
      dset2_filename	= 'NOA_ASCII_AMV-4thInt_TEST2.GOES16.2019293.1130.CH_14.FD.Cloud.nc'
      t_name = 'NOA_'+str(dset_type_str)
    elif shortabbrv=='NWC_':
      dset2_filename	= 'NWC_INTERCOMP2021_NWCSAFTEST21.Cloud.nc'
      t_name = 'NWC_'+str(dset_type_str)
  elif dset_type_str == 'Clear':
    if shortabbrv=='BRZ_':
      dset2_filename	= 'BRZ_4th_AMVIC_INPE_Test_3_final.Clear.nc'
      t_name = 'BRZ_'+str(dset_type_str)
    elif shortabbrv=='EUM_':
      dset2_filename	= 'EUM_AMVIntm_Chan08_20191020113000Z_GOES_ASCII_Test31.Clear.nc'
      t_name = 'EUM_'+str(dset_type_str)
    elif shortabbrv=='KMA_':
      dset2_filename	= 'KMA_test31_goes16_abi_ch08_amv_201910201130.Clear.nc'
      t_name = 'KMA_'+str(dset_type_str)
    elif shortabbrv=='NOA_':
      dset2_filename	= 'NOA_ASCII_AMV-4thInt_TEST3.GOES16.2019293.1130.CH_08.FD.Clear.nc'
      t_name = 'NOA_'+str(dset_type_str)
    elif shortabbrv=='NWC_':
      dset2_filename	= 'NWC_INTERCOMP2021_NWCSAFTEST31.Clear.nc'
      t_name = 'NWC_'+str(dset_type_str)
    elif shortabbrv=='BRZ4':
      dset2_filename    = 'BRZ_4th_AMVIC_INPE_Test_4_final.Clear.nc'
      t_name = 'BRZ4_'+str(dset_type_str)
    elif shortabbrv=='EUM4':
      dset2_filename    = 'EUM_AMVIntm_Chan08_20191020113000Z_GOES_ASCII_Test41.Clear.nc'
      t_name = 'EUM4_'+str(dset_type_str)
    elif shortabbrv=='KMA4':
      dset2_filename    = 'KMA_test41_goes16_abi_ch08_amv_201910201130.Clear.nc'
      t_name = 'KMA4_'+str(dset_type_str)
    elif shortabbrv=='NOA4':
      dset2_filename    = 'NOA_ASCII_AMV-4thInt_TEST4.GOES16.2019293.1130.CH_08.FD.Clear.nc'
      t_name = 'NOA4_'+str(dset_type_str)
    elif shortabbrv=='NWC4':
      dset2_filename    = 'NWC_INTERCOMP2021_NWCSAFTEST41.Clear.nc'
      t_name = 'NWC4_'+str(dset_type_str)

  dset2_file   		= dset2_path+dset2_filename

  dset2_exists = exists(dset2_file)
  if dset2_exists==False:
    print("ERROR: file "+dset2_file+" does not exist!")
    sys.exit()

		# Path (for output NetCDF only)
  str_dset2_path = dset2_path

  dset2_src = str_dset2_path+dset2_filename

		# Variable names
  dset2_lat_var  = 'Ob_Lat'
  dset2_lon_var  = 'Ob_Lon'
  dset2_prs_var  = 'Ob_Prs'		# units = hPa
  dset2_yr_var   = 'Ob_Year'
  dset2_mm_var   = 'Ob_Month'
  dset2_dy_var   = 'Ob_Day'
  dset2_hr_var   = 'Ob_Hour'
  dset2_mn_var   = 'Ob_Min'
  dset2_pccf_var = 'Ob_QI'
  
  dset2_spd_var	 = 'Ob_Spd'
  dset2_dir_var  = 'Ob_Dir'

    	#-------------------------------------------------
  	# Load dataset

  data_hdl = Dataset(dset2_file)

  tdset2_lat = np.asarray( data_hdl.variables[dset2_lat_var] )
  tdset2_lon = np.asarray( data_hdl.variables[dset2_lon_var] )		#NOTE: lons are [-180,180]. This will be converted to [0,360] automatically before collocation.
  tdset2_prs = np.asarray( data_hdl.variables[dset2_prs_var] )
  tdset2_yr  = np.asarray( data_hdl.variables[dset2_yr_var]  )
  tdset2_mm  = np.asarray( data_hdl.variables[dset2_mm_var]  )
  tdset2_dy  = np.asarray( data_hdl.variables[dset2_dy_var]  )
  tdset2_hr  = np.asarray( data_hdl.variables[dset2_hr_var]  )
  tdset2_mn  = np.asarray( data_hdl.variables[dset2_mn_var]  )
  tdset2_pccf = np.asarray( data_hdl.variables[dset2_pccf_var]  )

  tdset2_spd = np.asarray( data_hdl.variables[dset2_spd_var] )
  tdset2_dir = np.asarray( data_hdl.variables[dset2_dir_var] )

  data_hdl.close()

	# check pressure units and convert to hPa
  if max(tdset2_prs) > 10000.:
    tdset2_prs = tdset2_prs/100.

  sindexesDC = np.asarray(np.where(tdset2_lat==tdset2_lat))	      #get all indices
  indexes2   = sindexesDC.flatten()

  qc_list = "No QC applied"

  d_lat = tdset2_lat
  d_lon = tdset2_lon
  d_yr  = tdset2_yr 
  d_mm  = tdset2_mm 
  d_dy  = tdset2_dy 
  d_hr  = tdset2_hr 
  d_mn  = tdset2_mn  
  d_prs = tdset2_prs
  
  d_spd = tdset2_spd
  d_dir = tdset2_dir

	# create height array but fill with missing -999.
	#	AMV heights not available.
  d_hgt = np.nan * np.ones_like(d_lat)

  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes2,qc_list,dset2_src, d_spd,d_dir

# -------------------------------------------------------------------------
# Read AMV (from NCEP)
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		pct ................................. Minimum AMV quality indicator (QI) in % for QC
#		bool_qc ............................. Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		d_pccf .............................. AMV QI (percent confidence) in %
#               indexes2 ............................ Indices of input obs
#		qc_list ............................. List of QC applied (if applicable)
#
def read_amv_ncep(path_prefix,yyyymmddhh,bool_qc,pct,qi_choice):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

	#-------------------------------------------------
    	# Define dataset

  tmp_dset2_path = '/atmos-nc-dataset/AMV/'+yyyy+'/'+mm+'/'+dd+'/'

		# Path/file
  dset2_path   		= path_prefix+tmp_dset2_path

  dset2_filename 	= 'gdas.'+yyyymmddhh+'.satwnd.tm00.bufr_d.nc4'

  dset2_file   		= dset2_path+dset2_filename

  dset2_exists = exists(dset2_file)
  if dset2_exists==False:
    print("ERROR: file "+dset2_file+" does not exist!")
    sys.exit()

		# Path on FTP/web archive server (for output NetCDF only)
  str_dset2_path = dset2_path

  dset2_src = str_dset2_path+dset2_filename

		# Variable names
  dset2_lat_var  = 'latitude'
  dset2_lon_var  = 'longitude'
  dset2_prs_var  = 'pressure'
  dset2_yr_var   = 'year'
  dset2_mm_var   = 'month'
  dset2_dy_var   = 'day'
  dset2_hr_var   = 'hour'
  dset2_mn_var   = 'minutes'
  if qi_choice == 'YES_FC':  
    dset2_pccf_var = 'percent_confidence_yes_forecast'
  else:
    dset2_pccf_var = 'percent_confidence_no_forecast'

  dset2_satname_var = 'satellite_name'
  dset2_wcm_var  = 'wind_calculation_method'
  dset2_ham_var  = 'height_assignment_method'
  dset2_spd_var  = 'wind_speed'
  dset2_dir_var  = 'wind_direction'  
		
    	#-------------------------------------------------
  	# Load dataset

  data_hdl = Dataset(dset2_file)

  tdset2_lat = np.asarray( data_hdl.variables[dset2_lat_var] )
  tdset2_lon = np.asarray( data_hdl.variables[dset2_lon_var] )
  tdset2_prs = np.asarray( data_hdl.variables[dset2_prs_var] )
  tdset2_yr  = np.asarray( data_hdl.variables[dset2_yr_var]  )
  tdset2_mm  = np.asarray( data_hdl.variables[dset2_mm_var]  )
  tdset2_dy  = np.asarray( data_hdl.variables[dset2_dy_var]  )
  tdset2_hr  = np.asarray( data_hdl.variables[dset2_hr_var]  )
  tdset2_mn  = np.asarray( data_hdl.variables[dset2_mn_var]  )
  tdset2_pccf = np.asarray( data_hdl.variables[dset2_pccf_var]  )
  
  tdset2_satname = np.asarray( data_hdl.variables[dset2_satname_var]  )
  tdset2_wcm = np.asarray( data_hdl.variables[dset2_wcm_var]  )
  tdset2_ham = np.asarray( data_hdl.variables[dset2_ham_var]  )
  tdset2_spd = np.asarray( data_hdl.variables[dset2_spd_var]  )
  tdset2_dir = np.asarray( data_hdl.variables[dset2_dir_var]  )

  data_hdl.close()

	# check pressure units and convert to hPa
  if max(tdset2_prs) > 10000.:
    tdset2_prs = tdset2_prs/100.

	# apply QC if bool_dset2_qc=True
  if bool_qc:

    tindexes2,qc_list = qc_amv(tdset2_pccf,pct)
    
    sindexes2 = np.asarray(tindexes2)
    indexes2  = sindexes2.flatten()
 
    d_lat = tdset2_lat[indexes2]
    d_lon = tdset2_lon[indexes2]
    d_prs = tdset2_prs[indexes2]
    d_yr  = tdset2_yr [indexes2]
    d_mm  = tdset2_mm [indexes2]
    d_dy  = tdset2_dy [indexes2]
    d_hr  = tdset2_hr [indexes2]
    d_mn  = tdset2_mn [indexes2]
    
    d_satname = tdset2_satname[indexes2]
    d_wcm = tdset2_wcm[indexes2]
    d_ham = tdset2_ham[indexes2]
    d_spd = tdset2_spd[indexes2]
    d_dir = tdset2_dir[indexes2]

  elif not bool_qc:
    print("Do not apply AMV_NCEP QC")

    sindexesDC = np.asarray(np.where(tdset2_lat==tdset2_lat))           #get all indices
    indexes2   = sindexesDC.flatten()

    qc_list = "No QC applied"

    d_lat = tdset2_lat
    d_lon = tdset2_lon
    d_yr  = tdset2_yr 
    d_mm  = tdset2_mm 
    d_dy  = tdset2_dy 
    d_hr  = tdset2_hr 
    d_mn  = tdset2_mn  
    d_prs = tdset2_prs
    
    d_satname = tdset2_satname
    d_wcm = tdset2_wcm
    d_ham = tdset2_ham
    d_spd = tdset2_spd
    d_dir = tdset2_dir

	# create height array but fill with missing -999.
	#	AMV heights not available.
  d_hgt = np.nan * np.ones_like(d_lat)

  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes2,qc_list,dset2_src, d_satname,d_wcm,d_ham,d_spd,d_dir
  
# -------------------------------------------------------------------------
# Read AMV (from NCEP)
#	Don't QC
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		pct ................................. Minimum AMV quality indicator (QI) in % for QC
#		bool_qc ............................. Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		d_pccf .............................. AMV QI (percent confidence) in %
#               indexes2 ............................ Indices of input obs
#		qc_list ............................. List of QC applied (if applicable)
#
def read_amv_ncep_for_plotting(path_prefix,yyyymmddhh,bool_qc,pct,qi_choice):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

	#-------------------------------------------------
    	# Define dataset

  tmp_dset2_path = '/atmos-nc-dataset/AMV/'+yyyy+'/'+mm+'/'+dd+'/'

		# Path/file
  dset2_path   		= path_prefix+tmp_dset2_path

  dset2_filename 	= 'gdas.'+yyyymmddhh+'.satwnd.tm00.bufr_d.nc4'

  dset2_file   		= dset2_path+dset2_filename

  dset2_exists = exists(dset2_file)
  if dset2_exists==False:
    print("ERROR: file "+dset2_file+" does not exist!")
    sys.exit()

		# Path on FTP/web archive server (for output NetCDF only)
  str_dset2_path = dset2_path

  dset2_src = str_dset2_path+dset2_filename

		# Variable names
  dset2_lat_var  = 'latitude'
  dset2_lon_var  = 'longitude'
  dset2_prs_var  = 'pressure'
  dset2_yr_var   = 'year'
  dset2_mm_var   = 'month'
  dset2_dy_var   = 'day'
  dset2_hr_var   = 'hour'
  dset2_mn_var   = 'minutes'
  if qi_choice == 'YES_FC':  
    dset2_pccf_var = 'percent_confidence_yes_forecast'
  else:
    dset2_pccf_var = 'percent_confidence_no_forecast'

  dset2_satname_var = 'satellite_name'
  dset2_wcm_var  = 'wind_calculation_method'
  dset2_ham_var  = 'height_assignment_method'
  dset2_spd_var  = 'wind_speed'
  dset2_dir_var  = 'wind_direction'  
		
    	#-------------------------------------------------
  	# Load dataset

  data_hdl = Dataset(dset2_file)

  tdset2_lat = np.asarray( data_hdl.variables[dset2_lat_var] )
  tdset2_lon = np.asarray( data_hdl.variables[dset2_lon_var] )
  tdset2_prs = np.asarray( data_hdl.variables[dset2_prs_var] )
  tdset2_yr  = np.asarray( data_hdl.variables[dset2_yr_var]  )
  tdset2_mm  = np.asarray( data_hdl.variables[dset2_mm_var]  )
  tdset2_dy  = np.asarray( data_hdl.variables[dset2_dy_var]  )
  tdset2_hr  = np.asarray( data_hdl.variables[dset2_hr_var]  )
  tdset2_mn  = np.asarray( data_hdl.variables[dset2_mn_var]  )
  tdset2_pccf = np.asarray( data_hdl.variables[dset2_pccf_var]  )
  
  tdset2_satname = np.asarray( data_hdl.variables[dset2_satname_var]  )
  tdset2_wcm = np.asarray( data_hdl.variables[dset2_wcm_var]  )
  tdset2_ham = np.asarray( data_hdl.variables[dset2_ham_var]  )
  tdset2_spd = np.asarray( data_hdl.variables[dset2_spd_var]  )
  tdset2_dir = np.asarray( data_hdl.variables[dset2_dir_var]  )

  data_hdl.close()

	# check pressure units and convert to hPa
  if max(tdset2_prs) > 10000.:
    tdset2_prs = tdset2_prs/100.


  sindexesDC = np.asarray(np.where(tdset2_lat==tdset2_lat))	      #get all indices
  indexes2   = sindexesDC.flatten()

  qc_list = "No QC applied"

  d_lat = tdset2_lat
  d_lon = tdset2_lon
  d_yr  = tdset2_yr 
  d_mm  = tdset2_mm 
  d_dy  = tdset2_dy 
  d_hr  = tdset2_hr 
  d_mn  = tdset2_mn  
  d_prs = tdset2_prs
  
  d_satname = tdset2_satname
  d_wcm = tdset2_wcm
  d_ham = tdset2_ham
  d_spd = tdset2_spd
  d_dir = tdset2_dir

	# create height array but fill with missing -999.
	#	AMV heights not available.
  d_hgt = np.nan * np.ones_like(d_lat)

  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes2,qc_list,dset2_src, d_satname,d_wcm,d_ham,d_spd,d_dir

# -------------------------------------------------------------------------
# Read Loon stratospheric balloon data
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		bool_qc ......................... Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		qc_list ............................. List of QC applied (if applicable)
#
def read_loon(path_prefix,yyyymmddhh,bool_qc):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

	#-------------------------------------------------
    	# Define dataset

  tmp_dset1_path = '/atmos-nc-dataset/Loon/'+yyyy+'/'+mm+'/'+dd+'/'

		# Path/file
  dset1_path   		= path_prefix+tmp_dset1_path

  dset1_filename 	= 'Loon_'+yyyymmddhh+'.nc4'

  dset1_file   		= dset1_path+dset1_filename

  dset1_exists = exists(dset1_file)
  if dset1_exists==False:
    print("ERROR: file "+dset1_file+" does not exist!")
    sys.exit()							#exit script immediately

		# Path on FTP/web archive server (for output NetCDF only)
  str_dset1_path = dset1_path

  dset1_src = str_dset1_path+dset1_filename

		# Variable names
  dset1_lat_var = 'latitude'
  dset1_lon_var = 'longitude'
  dset1_prs_var = 'pressure'
  dset1_yr_var  = 'year'
  dset1_mm_var  = 'month'
  dset1_dy_var  = 'day'
  dset1_hr_var  = 'hour'
  dset1_mn_var  = 'minute'
  dset1_hgt_var = 'altitude'
  
  dset1_azm_var = 'solar_azimuth_angle'
  dset1_elv_var = 'solar_elevation_angle'
  dset1_spd_var = 'wind_speed'
  dset1_dir_var = 'wind_direction'
		
    	#-------------------------------------------------
  	# Load dataset
	#	AIRCRAFT data is divided into Groups
  
  data_hdl = Dataset(dset1_file)

  ttdset1_lat = np.asarray( data_hdl.variables[dset1_lat_var] )
  ttdset1_lon = np.asarray( data_hdl.variables[dset1_lon_var] )
  ttdset1_prs = np.asarray( data_hdl.variables[dset1_prs_var] )
  ttdset1_yr  = np.asarray( data_hdl.variables[dset1_yr_var]  )
  ttdset1_mm  = np.asarray( data_hdl.variables[dset1_mm_var]  )
  ttdset1_dy  = np.asarray( data_hdl.variables[dset1_dy_var]  )
  ttdset1_hr  = np.asarray( data_hdl.variables[dset1_hr_var]  )
  ttdset1_mn  = np.asarray( data_hdl.variables[dset1_mn_var]  )
  ttdset1_hgt = np.asarray( data_hdl.variables[dset1_hgt_var]  )
  
  ttdset1_azm = np.asarray( data_hdl.variables[dset1_azm_var]  )
  ttdset1_elv = np.asarray( data_hdl.variables[dset1_elv_var]  )
  ttdset1_spd = np.asarray( data_hdl.variables[dset1_spd_var]  )
  ttdset1_dir = np.asarray( data_hdl.variables[dset1_dir_var]  )

  data_hdl.close()

	# check if arrays are missing and omit them
  idx = np.where(ttdset1_yr > -999)
  tdset1_lat = ttdset1_lat[idx]
  tdset1_lon = ttdset1_lon[idx]
  tdset1_yr  = ttdset1_yr [idx]
  tdset1_mm  = ttdset1_mm [idx]
  tdset1_dy  = ttdset1_dy [idx]
  tdset1_hr  = ttdset1_hr [idx]
  tdset1_mn  = ttdset1_mn [idx]
  tdset1_hgt = ttdset1_hgt[idx]
  tdset1_prs = ttdset1_prs[idx]
  tdset1_azm = ttdset1_azm[idx]
  tdset1_elv = ttdset1_elv[idx]
  tdset1_spd = ttdset1_spd[idx]
  tdset1_dir = ttdset1_dir[idx]
  del idx

	# check pressure units and convert to hPa
  if max(tdset1_prs) > 10000.:
    tdset1_prs = tdset1_prs/100.
    
    	# check height units and convert to km
  print("tdset1_hgt before = "+str(tdset1_hgt))
  if max(tdset1_hgt) > 1000.:
    tdset1_hgt = tdset1_hgt/1000.

	# apply QC if bool_dset2_qc=True
  if bool_qc:
    print("LOON QC NOT APPLIED: TBD")
#    tindexes1,qc_list = qc_amv(tdset1_pccf,pct)
#    
#    sindexes1 = np.asarray(tindexes1)
#    indexes1  = sindexes2.flatten()
# 
#    d_lat = tdset1_lat[indexes1]
#    d_lon = tdset1_lon[indexes1]
#    d_prs = tdset1_prs[indexes1]
#    d_yr  = tdset1_yr [indexes1]
#    d_mm  = tdset1_mm [indexes1]
#    d_dy  = tdset1_dy [indexes1]
#    d_hr  = tdset1_hr [indexes1]
#    d_mn  = tdset1_mn [indexes1]
#    d_hgt = tdset1_hgt[indexes1]
#    d_azm = tdset1_azm[indexes1]
#    d_elv = tdset1_elv[indexes1]
#    d_spd = tdset1_spd[indexes1]
#    d_dir = tdset1_dir[indexes1]

  elif not bool_qc:
    print("Do not apply LOON QC")

    sindexesDC = np.asarray(np.where(tdset1_lat==tdset1_lat))           #get all indices
    indexes1   = sindexesDC.flatten()

    qc_list = "No QC applied"

    d_lat = tdset1_lat
    d_lon = tdset1_lon
    d_yr  = tdset1_yr 
    d_mm  = tdset1_mm 
    d_dy  = tdset1_dy 
    d_hr  = tdset1_hr 
    d_mn  = tdset1_mn  
    d_prs = tdset1_prs
    d_hgt = tdset1_hgt
    d_azm = tdset1_azm
    d_elv = tdset1_elv
    d_spd = tdset1_spd
    d_dir = tdset1_dir

  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes1,qc_list,dset1_src,d_azm,d_elv,d_spd,d_dir
  
# -------------------------------------------------------------------------
# Read Loon stratospheric balloon data
#	Don't QC
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		bool_qc ......................... Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		qc_list ............................. List of QC applied (if applicable)
#
def read_loon_for_plotting(path_prefix,yyyymmddhh,bool_qc):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

	#-------------------------------------------------
    	# Define dataset

  tmp_dset1_path = '/atmos-nc-dataset/Loon/'+yyyy+'/'+mm+'/'+dd+'/'

		# Path/file
  dset1_path   		= path_prefix+tmp_dset1_path

  dset1_filename 	= 'Loon_'+yyyymmddhh+'.nc4'

  dset1_file   		= dset1_path+dset1_filename

  dset1_exists = exists(dset1_file)
  if dset1_exists==False:
    print("ERROR: file "+dset1_file+" does not exist!")
    sys.exit()							#exit script immediately

		# Path on FTP/web archive server (for output NetCDF only)
  str_dset1_path = dset1_path

  dset1_src = str_dset1_path+dset1_filename

		# Variable names
  dset1_lat_var = 'latitude'
  dset1_lon_var = 'longitude'
  dset1_prs_var = 'pressure'
  dset1_yr_var  = 'year'
  dset1_mm_var  = 'month'
  dset1_dy_var  = 'day'
  dset1_hr_var  = 'hour'
  dset1_mn_var  = 'minute'
  dset1_hgt_var = 'altitude'
  
  dset1_azm_var = 'solar_azimuth_angle'
  dset1_elv_var = 'solar_elevation_angle'
  dset1_spd_var = 'wind_speed'
  dset1_dir_var = 'wind_direction'
		
    	#-------------------------------------------------
  	# Load dataset
	#	AIRCRAFT data is divided into Groups
  
  data_hdl = Dataset(dset1_file)

  ttdset1_lat = np.asarray( data_hdl.variables[dset1_lat_var] )
  ttdset1_lon = np.asarray( data_hdl.variables[dset1_lon_var] )
  ttdset1_prs = np.asarray( data_hdl.variables[dset1_prs_var] )
  ttdset1_yr  = np.asarray( data_hdl.variables[dset1_yr_var]  )
  ttdset1_mm  = np.asarray( data_hdl.variables[dset1_mm_var]  )
  ttdset1_dy  = np.asarray( data_hdl.variables[dset1_dy_var]  )
  ttdset1_hr  = np.asarray( data_hdl.variables[dset1_hr_var]  )
  ttdset1_mn  = np.asarray( data_hdl.variables[dset1_mn_var]  )
  ttdset1_hgt = np.asarray( data_hdl.variables[dset1_hgt_var]  )
  
  ttdset1_azm = np.asarray( data_hdl.variables[dset1_azm_var]  )
  ttdset1_elv = np.asarray( data_hdl.variables[dset1_elv_var]  )
  ttdset1_spd = np.asarray( data_hdl.variables[dset1_spd_var]  )
  ttdset1_dir = np.asarray( data_hdl.variables[dset1_dir_var]  )

  data_hdl.close()

	# check if arrays are missing and omit them
  idx = np.where(ttdset1_yr > -999)
  tdset1_lat = ttdset1_lat[idx]
  tdset1_lon = ttdset1_lon[idx]
  tdset1_yr  = ttdset1_yr [idx]
  tdset1_mm  = ttdset1_mm [idx]
  tdset1_dy  = ttdset1_dy [idx]
  tdset1_hr  = ttdset1_hr [idx]
  tdset1_mn  = ttdset1_mn [idx]
  tdset1_hgt = ttdset1_hgt[idx]
  tdset1_prs = ttdset1_prs[idx]
  tdset1_azm = ttdset1_azm[idx]
  tdset1_elv = ttdset1_elv[idx]
  tdset1_spd = ttdset1_spd[idx]
  tdset1_dir = ttdset1_dir[idx]
  del idx

	# check pressure units and convert to hPa
  if max(tdset1_prs) > 10000.:
    tdset1_prs = tdset1_prs/100.
    
    	# check height units and convert to km
  print("tdset1_hgt before = "+str(tdset1_hgt))
  if max(tdset1_hgt) > 1000.:
    tdset1_hgt = tdset1_hgt/1000.

  sindexesDC = np.asarray(np.where(tdset1_lat==tdset1_lat))	      #get all indices
  indexes1   = sindexesDC.flatten()

  qc_list = "No QC applied"

  d_lat = tdset1_lat
  d_lon = tdset1_lon
  d_yr  = tdset1_yr 
  d_mm  = tdset1_mm 
  d_dy  = tdset1_dy 
  d_hr  = tdset1_hr 
  d_mn  = tdset1_mn  
  d_prs = tdset1_prs
  d_hgt = tdset1_hgt
  d_azm = tdset1_azm
  d_elv = tdset1_elv
  d_spd = tdset1_spd
  d_dir = tdset1_dir
 
  print("min/max loon: "+str(min(d_spd))+","+str(max(d_spd)))

  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes1,qc_list,dset1_src,d_azm,d_elv,d_spd,d_dir

# -------------------------------------------------------------------------
# Read RADIOSONDE (from NCEP)
#
#	INPUTS:
#		path_prefix ......................... Full path up to directories /aeolus... or /atmos... : where this archive is located on local machine
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		bool_qc ............................. Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		qc_list ............................. List of QC applied (if applicable)
#
def read_raobs(path_prefix,yyyymmddhh,bool_qc):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

	#-------------------------------------------------
    	# Define dataset

  tmp_dset1_path = '/atmos-nc-dataset/radiosonde/'+yyyy+'/'+mm+'/'+dd+'/'

		# Path/file
  dset1_path   		= path_prefix+tmp_dset1_path

  dset1_filename 	= 'gdas.'+yyyymmddhh+'.adpupa.tm00.bufr_d.nc4'

  dset1_file   		= dset1_path+dset1_filename

  dset1_exists = exists(dset1_file)
  if dset1_exists==False:
    print("ERROR: file "+dset1_file+" does not exist!")
    sys.exit()							#exit script immediately

		# Path on FTP/web archive server (for output NetCDF only)
  str_dset1_path = dset1_path

  dset1_src = str_dset1_path+dset1_filename

		# Variable names
  dset1_lat_var = 'latitude'
  dset1_lon_var = 'longitude'
  dset1_yr_var  = 'year'
  dset1_mm_var  = 'month'
  dset1_dy_var  = 'day'
  dset1_hr_var  = 'hour'
  dset1_mn_var  = 'minutes'
  dset1_prs_var = 'pressure'
  dset1_ht_var  = 'height'
  
  dset1_spd_var = 'wind_speed'
  dset1_dir_var = 'wind_direction'
		
    	#-------------------------------------------------
  	# Load dataset
	#	RADIOSONDE data is divided into Groups
 
  nsondes = []
  nlevels = []
 
  data_hdl = Dataset(dset1_file)

  grps = list(data_hdl.groups)
  ngrps = np.size(grps)
  print("RAOBS groups list = "+str(grps))
  print("RAOBS ngroups = "+str(ngrps))

	# populate full arrays with Group1 (grps(0))
  td_lat = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lat_var] )
  td_lon = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lon_var] )
  td_yr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_yr_var] )
  td_mm  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mm_var] )
  td_dy  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dy_var] )
  td_hr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var] )
  
  td_spd = np.asarray( data_hdl.groups[grps[0]].variables[dset1_spd_var] )
  td_dir = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dir_var] )

  print("RAOBS lat size, shape: "+str(np.size(td_lat))+" | "+str(np.shape(td_lat)))

  if grps[0]=='NC002001' or grps[0]=='NC002002' or grps[0]=='NC002003' or grps[0]=='NC002004' or grps[0]=='NC002005' or grps[0]=='NC002009':
    td_mn = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var]  )
    td_mn[:] = 0.0
  else:
    td_mn = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mn_var] )
  
  if grps[0]=='NC002001' or grps[0]=='NC002002' or grps[0]=='NC002003' or grps[0]=='NC002004' or grps[0]=='NC002005' or grps[0]=='NC002006' or grps[0]=='NC002009' or grps[0]=='NC002015':
    	# groups with pressure variable but not height
    td_prs = np.asarray( data_hdl.groups[grps[0]].variables[dset1_prs_var]  )
    td_hgt = np.nan * np.ones_like(td_prs)
  else:
    	# groups with height variable but not pressure
    td_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht_var]  )
    td_prs = np.nan * np.ones_like(td_hgt)

	# find shape of array (nsondes, nlevels)
  grp_shape = np.shape(td_prs)
  nsondes.append(grp_shape[0])
  nlevels.append(grp_shape[1])

	# assign same lat,lon,yr,mm,dy,hr,mn to all levels per sonde (so that all vars have same dimensions)
  ttd_lat = np.nan * np.ones_like(td_prs)
  ttd_lon = np.nan * np.ones_like(td_prs)
  ttd_yr  = np.nan * np.ones_like(td_prs)
  ttd_mm  = np.nan * np.ones_like(td_prs)
  ttd_dy  = np.nan * np.ones_like(td_prs)
  ttd_hr  = np.nan * np.ones_like(td_prs)
  ttd_mn  = np.nan * np.ones_like(td_prs)
  ttd_prs = np.nan * np.ones_like(td_prs)
  ttd_hgt = np.nan * np.ones_like(td_prs)
  
  ttd_spd = np.nan * np.ones_like(td_prs)
  ttd_dir = np.nan * np.ones_like(td_prs)

  for isonde in range(grp_shape[0]):
    ttd_lat[isonde,:] = td_lat[isonde]
    ttd_lon[isonde,:] = td_lon[isonde]
    ttd_yr [isonde,:] = td_yr [isonde]
    ttd_mm [isonde,:] = td_mm [isonde]
    ttd_dy [isonde,:] = td_dy [isonde]
    ttd_hr [isonde,:] = td_hr [isonde]
    ttd_mn [isonde,:] = td_mn [isonde]
    ttd_prs[isonde,:] = td_prs[isonde]
    ttd_hgt[isonde,:] = td_hgt[isonde]
    
    ttd_spd[isonde,:] = td_spd[isonde]
    ttd_dir[isonde,:] = td_dir[isonde]
    

	# conform arrays to 1D
	#	.flatten() COPIES the original array and conforms its dimensions to 1D
	#	Ex) a = [[1,2,3], [4,5,6]]
	#	    b = a.flatten()
	#	    print(b) --> prints b = [1,2,3,4,5,6]
  fd_lat = ttd_lat.flatten()
  fd_lon = ttd_lon.flatten()
  fd_yr  = ttd_yr.flatten()
  fd_mm  = ttd_mm.flatten()
  fd_dy  = ttd_dy.flatten()
  fd_hr  = ttd_hr.flatten()
  fd_mn  = ttd_mn.flatten()
  fd_prs = ttd_prs.flatten()
  fd_hgt = ttd_hgt.flatten()
  
  fd_spd = ttd_spd.flatten()
  fd_dir = ttd_dir.flatten()

	# append data from remaining Groups
  for x in range(len(grps)-1):
#      print("grp exists = "+grps[x+1])
      tdset1_lat = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lat_var] )
      tdset1_lon = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lon_var] )
      tdset1_yr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_yr_var]  )
      tdset1_mm  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mm_var]  )
      tdset1_dy  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dy_var]  )
      tdset1_hr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
      
      tdset1_spd = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_spd_var]  )
      tdset1_dir = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dir_var]  )

      if grps[x+1]=='NC002001' or grps[x+1]=='NC002002' or grps[x+1]=='NC002003' or grps[x+1]=='NC002004' or grps[x+1]=='NC002005' or grps[x+1]=='NC002009':
      	# groups without minutes variable
        tdset1_mn = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
        tdset1_mn[:] = 0.0
      else:
      	# groups with minutes variable
        tdset1_mn = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mn_var]  )

      if grps[x+1]=='NC002001' or grps[x+1]=='NC002002' or grps[x+1]=='NC002003' or grps[x+1]=='NC002004' or grps[x+1]=='NC002005' or grps[x+1]=='NC002006' or grps[x+1]=='NC002009' or grps[x+1]=='NC002015':
      	# groups with pressure variable but not height
        p = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_prs_var]  )
        z = np.nan * np.ones_like(p)
        tgrp_shape = np.shape(p)
      else:
      	# groups with height variable but not pressure
        z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht_var]  )
        p = np.nan * np.ones_like(z)
        tgrp_shape = np.shape(z)

	# find shape of array (nsondes, nlevels)
      tnsondes   = tgrp_shape[0]
      tnlevels   = tgrp_shape[1]

      nsondes.append(tnsondes)
      nlevels.append(tnlevels)

	# assign same lat,lon,yr,mm,dy,hr,mn to all levels per sonde (so that all vars have same dimensions)
      ttdset1_lat = np.nan * np.ones_like(p)
      ttdset1_lon = np.nan * np.ones_like(p)
      ttdset1_yr  = np.nan * np.ones_like(p)
      ttdset1_mm  = np.nan * np.ones_like(p)
      ttdset1_dy  = np.nan * np.ones_like(p)
      ttdset1_hr  = np.nan * np.ones_like(p)
      ttdset1_mn  = np.nan * np.ones_like(p)
      ttdset1_prs = np.nan * np.ones_like(p)
      ttdset1_hgt = np.nan * np.ones_like(p)
      
      ttdset1_spd = np.nan * np.ones_like(p)
      ttdset1_dir = np.nan * np.ones_like(p)

      for isonde in range(tnsondes):
        ttdset1_lat[isonde,:] = tdset1_lat[isonde]
        ttdset1_lon[isonde,:] = tdset1_lon[isonde]
        ttdset1_yr [isonde,:] = tdset1_yr [isonde]
        ttdset1_mm [isonde,:] = tdset1_mm [isonde]
        ttdset1_dy [isonde,:] = tdset1_dy [isonde]
        ttdset1_hr [isonde,:] = tdset1_hr [isonde]
        ttdset1_mn [isonde,:] = tdset1_mn [isonde]
        ttdset1_prs[isonde,:] = p[isonde]
        ttdset1_hgt[isonde,:] = z[isonde]
	
        ttdset1_spd[isonde,:] = tdset1_spd[isonde]
        ttdset1_dir[isonde,:] = tdset1_dir[isonde]
 
	# conform arrays to 1D
        #       .flatten() COPIES the original array and conforms its dimensions to 1D
        #       Ex) a = [[1,2,3], [4,5,6]]
        #           b = a.flatten()
        #           print(b) --> prints b = [1,2,3,4,5,6]
      fdset1_lat = ttdset1_lat.flatten()
      fdset1_lon = ttdset1_lon.flatten()
      fdset1_yr  = ttdset1_yr.flatten()
      fdset1_mm  = ttdset1_mm.flatten()
      fdset1_dy  = ttdset1_dy.flatten()
      fdset1_hr  = ttdset1_hr.flatten()
      fdset1_mn  = ttdset1_mn.flatten()
      fdset1_prs = ttdset1_prs.flatten()
      fdset1_hgt = ttdset1_hgt.flatten()
      
      fdset1_spd = ttdset1_spd.flatten()
      fdset1_dir = ttdset1_dir.flatten()

	# append 1D arrays to first group's 1D arrays
      fd_lat = np.append(fd_lat,fdset1_lat,axis=0)
      fd_lon = np.append(fd_lon,fdset1_lon,axis=0)
      fd_yr  = np.append(fd_yr, fdset1_yr, axis=0)
      fd_mm  = np.append(fd_mm, fdset1_mm, axis=0)
      fd_dy  = np.append(fd_dy, fdset1_dy, axis=0)
      fd_hr  = np.append(fd_hr, fdset1_hr, axis=0)
      fd_mn  = np.append(fd_mn, fdset1_mn, axis=0)
      fd_prs = np.append(fd_prs,fdset1_prs,axis=0)
      fd_hgt = np.append(fd_hgt,fdset1_hgt,axis=0)
      
      fd_spd = np.append(fd_spd,fdset1_spd,axis=0)
      fd_dir = np.append(fd_dir,fdset1_dir,axis=0)

  data_hdl.close()

	# apply QC if bool_dset1_qc=True
  if bool_qc:
    print("RADIOSONDE bool_qc variable TBA")	
  elif not bool_qc:
    print("Do not apply RADIOSONDE QC")

    sindexesDC = np.asarray(np.where(fd_prs==fd_prs))           #get all indices
    indexes1   = sindexesDC.flatten()

    qc_list = "No QC applied"

    d_lat = fd_lat
    d_lon = fd_lon
    d_prs = fd_prs
    d_hgt = fd_hgt
    d_yr  = fd_yr
    d_mm  = fd_mm
    d_dy  = fd_dy
    d_hr  = fd_hr
    d_mn  = fd_mn
    
    d_spd = fd_spd
    d_dir = fd_dir

	# check pressure units and convert to hPa
  if max(d_prs) > 10000.:
    d_prs = d_prs/100.

	# check height units and convert to km
  if max(d_hgt) > 1000.:
    d_hgt = d_hgt/1000.

  print("RAOBS FINAL SHAPE: "+str(np.shape(d_prs)))
  print("min/max raobs: "+str(min(d_spd))+","+str(max(d_spd)))

	# Return variables to MAIN
  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes1,qc_list,dset1_src,nsondes,nlevels,ngrps, d_spd,d_dir
  
# -------------------------------------------------------------------------
# Read RADIOSONDE (from NCEP)
#	Don't QC
#
#	INPUTS:
#		path_prefix ......................... Full path up to directories /aeolus... or /atmos... : where this archive is located on local machine
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		bool_qc ............................. Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		qc_list ............................. List of QC applied (if applicable)
#
def read_raobs_for_plotting(path_prefix,yyyymmddhh,bool_qc):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

	#-------------------------------------------------
    	# Define dataset

  tmp_dset1_path = '/atmos-nc-dataset/radiosonde/'+yyyy+'/'+mm+'/'+dd+'/'

		# Path/file
  dset1_path   		= path_prefix+tmp_dset1_path

  dset1_filename 	= 'gdas.'+yyyymmddhh+'.adpupa.tm00.bufr_d.nc4'

  dset1_file   		= dset1_path+dset1_filename

  dset1_exists = exists(dset1_file)
  if dset1_exists==False:
    print("ERROR: file "+dset1_file+" does not exist!")
    sys.exit()							#exit script immediately

		# Path on FTP/web archive server (for output NetCDF only)
  str_dset1_path = dset1_path

  dset1_src = str_dset1_path+dset1_filename

		# Variable names
  dset1_lat_var = 'latitude'
  dset1_lon_var = 'longitude'
  dset1_yr_var  = 'year'
  dset1_mm_var  = 'month'
  dset1_dy_var  = 'day'
  dset1_hr_var  = 'hour'
  dset1_mn_var  = 'minutes'
  dset1_prs_var = 'pressure'
  dset1_ht_var  = 'height'
  
  dset1_spd_var = 'wind_speed'
  dset1_dir_var = 'wind_direction'
		
    	#-------------------------------------------------
  	# Load dataset
	#	RADIOSONDE data is divided into Groups
 
  nsondes = []
  nlevels = []
 
  data_hdl = Dataset(dset1_file)

  grps = list(data_hdl.groups)
  ngrps = np.size(grps)
  print("RAOBS groups list = "+str(grps))
  print("RAOBS ngroups = "+str(ngrps))

	# populate full arrays with Group1 (grps(0))
  td_lat = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lat_var] )
  td_lon = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lon_var] )
  td_yr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_yr_var] )
  td_mm  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mm_var] )
  td_dy  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dy_var] )
  td_hr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var] )
  
  td_spd = np.asarray( data_hdl.groups[grps[0]].variables[dset1_spd_var] )
  td_dir = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dir_var] )

  print("RAOBS lat size, shape: "+str(np.size(td_lat))+" | "+str(np.shape(td_lat)))

  if grps[0]=='NC002001' or grps[0]=='NC002002' or grps[0]=='NC002003' or grps[0]=='NC002004' or grps[0]=='NC002005' or grps[0]=='NC002009':
    td_mn = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var]  )
    td_mn[:] = 0.0
  else:
    td_mn = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mn_var] )
  
  if grps[0]=='NC002001' or grps[0]=='NC002002' or grps[0]=='NC002003' or grps[0]=='NC002004' or grps[0]=='NC002005' or grps[0]=='NC002006' or grps[0]=='NC002009' or grps[0]=='NC002015':
    	# groups with pressure variable but not height
    td_prs = np.asarray( data_hdl.groups[grps[0]].variables[dset1_prs_var]  )
    td_hgt = np.nan * np.ones_like(td_prs)
  else:
    	# groups with height variable but not pressure
    td_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht_var]  )
    td_prs = np.nan * np.ones_like(td_hgt)

	# find shape of array (nsondes, nlevels)
  grp_shape = np.shape(td_prs)
  nsondes.append(grp_shape[0])
  nlevels.append(grp_shape[1])

	# assign same lat,lon,yr,mm,dy,hr,mn to all levels per sonde (so that all vars have same dimensions)
  ttd_lat = np.nan * np.ones_like(td_prs)
  ttd_lon = np.nan * np.ones_like(td_prs)
  ttd_yr  = np.nan * np.ones_like(td_prs)
  ttd_mm  = np.nan * np.ones_like(td_prs)
  ttd_dy  = np.nan * np.ones_like(td_prs)
  ttd_hr  = np.nan * np.ones_like(td_prs)
  ttd_mn  = np.nan * np.ones_like(td_prs)
  ttd_prs = np.nan * np.ones_like(td_prs)
  ttd_hgt = np.nan * np.ones_like(td_prs)
  
  ttd_spd = np.nan * np.ones_like(td_prs)
  ttd_dir = np.nan * np.ones_like(td_prs)

  for isonde in range(grp_shape[0]):
    ttd_lat[isonde,:] = td_lat[isonde]
    ttd_lon[isonde,:] = td_lon[isonde]
    ttd_yr [isonde,:] = td_yr [isonde]
    ttd_mm [isonde,:] = td_mm [isonde]
    ttd_dy [isonde,:] = td_dy [isonde]
    ttd_hr [isonde,:] = td_hr [isonde]
    ttd_mn [isonde,:] = td_mn [isonde]
    ttd_prs[isonde,:] = td_prs[isonde]
    ttd_hgt[isonde,:] = td_hgt[isonde]
    
    ttd_spd[isonde,:] = td_spd[isonde]
    ttd_dir[isonde,:] = td_dir[isonde]
    #print("RAOBS spd: "+str(td_spd[isonde]))

	# conform arrays to 1D
	#	.flatten() COPIES the original array and conforms its dimensions to 1D
	#	Ex) a = [[1,2,3], [4,5,6]]
	#	    b = a.flatten()
	#	    print(b) --> prints b = [1,2,3,4,5,6]
  fd_lat = ttd_lat.flatten()
  fd_lon = ttd_lon.flatten()
  fd_yr  = ttd_yr.flatten()
  fd_mm  = ttd_mm.flatten()
  fd_dy  = ttd_dy.flatten()
  fd_hr  = ttd_hr.flatten()
  fd_mn  = ttd_mn.flatten()
  fd_prs = ttd_prs.flatten()
  fd_hgt = ttd_hgt.flatten()
  
  fd_spd = ttd_spd.flatten()
  fd_dir = ttd_dir.flatten()

	# append data from remaining Groups
  for x in range(len(grps)-1):
#      print("grp exists = "+grps[x+1])
      tdset1_lat = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lat_var] )
      tdset1_lon = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lon_var] )
      tdset1_yr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_yr_var]  )
      tdset1_mm  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mm_var]  )
      tdset1_dy  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dy_var]  )
      tdset1_hr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
      
      tdset1_spd = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_spd_var]  )
      tdset1_dir = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dir_var]  )

      if grps[x+1]=='NC002001' or grps[x+1]=='NC002002' or grps[x+1]=='NC002003' or grps[x+1]=='NC002004' or grps[x+1]=='NC002005' or grps[x+1]=='NC002009':
      	# groups without minutes variable
        tdset1_mn = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
        tdset1_mn[:] = 0.0
      else:
      	# groups with minutes variable
        tdset1_mn = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mn_var]  )

      if grps[x+1]=='NC002001' or grps[x+1]=='NC002002' or grps[x+1]=='NC002003' or grps[x+1]=='NC002004' or grps[x+1]=='NC002005' or grps[x+1]=='NC002006' or grps[x+1]=='NC002009' or grps[x+1]=='NC002015':
      	# groups with pressure variable but not height
        p = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_prs_var]  )
        z = np.nan * np.ones_like(p)
        tgrp_shape = np.shape(p)
      else:
      	# groups with height variable but not pressure
        z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht_var]  )
        p = np.nan * np.ones_like(z)
        tgrp_shape = np.shape(z)

	# find shape of array (nsondes, nlevels)
      tnsondes   = tgrp_shape[0]
      tnlevels   = tgrp_shape[1]

      nsondes.append(tnsondes)
      nlevels.append(tnlevels)

	# assign same lat,lon,yr,mm,dy,hr,mn to all levels per sonde (so that all vars have same dimensions)
      ttdset1_lat = np.nan * np.ones_like(p)
      ttdset1_lon = np.nan * np.ones_like(p)
      ttdset1_yr  = np.nan * np.ones_like(p)
      ttdset1_mm  = np.nan * np.ones_like(p)
      ttdset1_dy  = np.nan * np.ones_like(p)
      ttdset1_hr  = np.nan * np.ones_like(p)
      ttdset1_mn  = np.nan * np.ones_like(p)
      ttdset1_prs = np.nan * np.ones_like(p)
      ttdset1_hgt = np.nan * np.ones_like(p)
      
      ttdset1_spd = np.nan * np.ones_like(p)
      ttdset1_dir = np.nan * np.ones_like(p)

      for isonde in range(tnsondes):
        ttdset1_lat[isonde,:] = tdset1_lat[isonde]
        ttdset1_lon[isonde,:] = tdset1_lon[isonde]
        ttdset1_yr [isonde,:] = tdset1_yr [isonde]
        ttdset1_mm [isonde,:] = tdset1_mm [isonde]
        ttdset1_dy [isonde,:] = tdset1_dy [isonde]
        ttdset1_hr [isonde,:] = tdset1_hr [isonde]
        ttdset1_mn [isonde,:] = tdset1_mn [isonde]
        ttdset1_prs[isonde,:] = p[isonde]
        ttdset1_hgt[isonde,:] = z[isonde]
	
        ttdset1_spd[isonde,:] = tdset1_spd[isonde]
        ttdset1_dir[isonde,:] = tdset1_dir[isonde]
 
	# conform arrays to 1D
        #       .flatten() COPIES the original array and conforms its dimensions to 1D
        #       Ex) a = [[1,2,3], [4,5,6]]
        #           b = a.flatten()
        #           print(b) --> prints b = [1,2,3,4,5,6]
      fdset1_lat = ttdset1_lat.flatten()
      fdset1_lon = ttdset1_lon.flatten()
      fdset1_yr  = ttdset1_yr.flatten()
      fdset1_mm  = ttdset1_mm.flatten()
      fdset1_dy  = ttdset1_dy.flatten()
      fdset1_hr  = ttdset1_hr.flatten()
      fdset1_mn  = ttdset1_mn.flatten()
      fdset1_prs = ttdset1_prs.flatten()
      fdset1_hgt = ttdset1_hgt.flatten()
      
      fdset1_spd = ttdset1_spd.flatten()
      fdset1_dir = ttdset1_dir.flatten()

	# append 1D arrays to first group's 1D arrays
      fd_lat = np.append(fd_lat,fdset1_lat,axis=0)
      fd_lon = np.append(fd_lon,fdset1_lon,axis=0)
      fd_yr  = np.append(fd_yr, fdset1_yr, axis=0)
      fd_mm  = np.append(fd_mm, fdset1_mm, axis=0)
      fd_dy  = np.append(fd_dy, fdset1_dy, axis=0)
      fd_hr  = np.append(fd_hr, fdset1_hr, axis=0)
      fd_mn  = np.append(fd_mn, fdset1_mn, axis=0)
      fd_prs = np.append(fd_prs,fdset1_prs,axis=0)
      fd_hgt = np.append(fd_hgt,fdset1_hgt,axis=0)
      
      fd_spd = np.append(fd_spd,fdset1_spd,axis=0)
      fd_dir = np.append(fd_dir,fdset1_dir,axis=0)

  data_hdl.close()

  sindexesDC = np.asarray(np.where(fd_prs==fd_prs))	      #get all indices
  indexes1   = sindexesDC.flatten()

  qc_list = "No QC applied"

  d_lat = fd_lat
  d_lon = fd_lon
  d_prs = fd_prs
  d_hgt = fd_hgt
  d_yr  = fd_yr
  d_mm  = fd_mm
  d_dy  = fd_dy
  d_hr  = fd_hr
  d_mn  = fd_mn
  
  d_spd = fd_spd
  d_dir = fd_dir

	# check pressure units and convert to hPa
  if max(d_prs) > 10000.:
    d_prs = d_prs/100.

	# check height units and convert to km
  if max(d_hgt) > 1000.:
    d_hgt = d_hgt/1000.

  print("RAOBS FINAL SHAPE: "+str(np.shape(d_prs)))

  print("min/max raobs: "+str(min(d_spd))+","+str(max(d_spd)))

	# Return variables to MAIN
  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes1,qc_list,dset1_src,nsondes,nlevels,ngrps, d_spd,d_dir

# -------------------------------------------------------------------------
# Read Collocation Index File
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		pct ................................. Minimum AMV quality indicator (QI) in % for QC
#		bool_qc ............................. Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		d_pccf .............................. AMV QI (percent confidence) in %
#               indexes2 ............................ Indices of input obs
#		qc_list ............................. List of QC applied (if applicable)
#
def read_index_file(inpath,yyyymmddhh,idx_file_str,ndset):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

	#-------------------------------------------------
    	# Define dataset

#  tmp_dset2_path = '/atmos-nc-dataset/collocation/'+yyyy+'/'+mm+'/'+dd+'/'
#  tmp_dset_path = '/data/users/klukens/collocation/longterm_anl/'

	# Paths 
  #dset_path   		= path_prefix+tmp_dset_path
  dset_path   		= inpath

  dset_filename 	= 'index.'+yyyymmddhh+idx_file_str

  dset_file   		= dset_path+dset_filename

  dset_exists = exists(dset_file)
  if dset_exists==False:
    print("ERROR: index file "+dset_file+" does not exist!")
    sys.exit()
    
    	# Path on FTP/web archive server (for output NetCDF only)
  str_dset_path = dset_path

  dset_src = str_dset_path+dset_filename

	# Variables names
		# DEPENDENT dataset name
  dset_name_var = 'dset'+str(ndset)

    		# indices of matched obs
			# DRIVER indices
  idx_drv_dset_var = 'idx_drv_'
  			# DEPENDENT indices
  idx_dset_var     = 'idx_'

		# collocation differences per matched pair
  DT_match_var  = 'DT_match_drv_'
  DP_match_var  = 'DP_match_drv_'
  DPlog_match_var  = 'DPlog_match_drv_'
  HT_match_var  = 'HT_match_drv_'
  GCD_match_var = 'GCD_match_drv_'
  
    	#-------------------------------------------------
  	# Load file and extract data

	# load dataset
  data_hdl = Dataset(dset_file)

  	# Extract data 
		# DEPENDENT dataset names
  dset_name      = np.asarray( data_hdl.variables[dset_name_var] )
  for attname in data_hdl.variables[dset_name_var].ncattrs():
    if attname == 'short_name':
      dset_shortname = getattr(data_hdl.variables[dset_name_var],attname)
      print("dset_shortname:")
      print(dset_shortname)
		# indices of matched obs
  idx_drv_dset  = np.asarray( data_hdl.variables[idx_drv_dset_var+dset_name_var] )
  idx_dset      = np.asarray( data_hdl.variables[idx_dset_var+dset_name_var] )

		# collocation differences per matched pair
  if np.size(idx_drv_dset) > 0:
    DT_match		= np.asarray( data_hdl.variables[DT_match_var+dset_name_var] )
    GCD_match		= np.asarray( data_hdl.variables[GCD_match_var+dset_name_var] )  
    Vert_match,PHstr	= check_PHvar_exists(dset_file,DP_match_var+dset_name_var,HT_match_var+dset_name_var)
  else:
    DT_match 	= 0
    GCD_match 	= 0
    Vert_match 	= 0
    PHstr 	= "NO_MATCHES"

  data_hdl.close()

  return dset_shortname,idx_drv_dset,idx_dset,DT_match,GCD_match,Vert_match,PHstr

# -------------------------------------------------------------------------

