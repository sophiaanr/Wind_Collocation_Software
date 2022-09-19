#########################################################################################################
# PLOT the Collocation of Winds from Multiple Datasets
# 	-- This script has been specifically modified to collocate Aircraft winds (dataset 1) and AMVs (dataset 2) with Aeolus winds (driver dataset)
#
# Developed by: Brett Hoover, UW-Madison/CIMSS
#	        David Santek, UW-Madison/CIMSS
#
# Modified by:  Katherine E. Lukens, UMD/CISESS, NOAA/NESDIS/STAR
#
# History: 
#	2021		B. Hoover		Created/developed program.
#	2021-10-25	K.E. Lukens		Added aircraft, AMV input configs. Added quality controls (QC).
#	2021-11-22	K.E. Lukens		Added date before consideration.
#	2021-11-23	K.E. Lukens		Created QC module and added relevant calls.
#	2022-03-08	K.E. Lukens		Added metadata to output files.
#
# Usage: python3 ThisScriptName.py $driver_wind_type $TF_drv_qc $TF_dset2_qc $pct $dateIN $driver_dset_type_str
#
#	driver_wind_type ..................... Aeolus wind type: RayClear = Rayleigh-clear winds, MieCloud = Mie-cloudy winds
#	TF_drv_qc ............................ Choose to apply Aeolus QC as listed in ./quality_controls.py: 0 = no QC, 1 = apply QC
#	TF_dset2_qc .......................... Choose to apply AMV QC as listed in ./quality_controls.py: 0 = no QC, 1 = apply QC
#	pct .................................. Percent confidence (quality indicator in %) for AMV data. If TF_dset2_qc=0, this value will be ignored
#	dateIN ............................... Date over which to perform collocation. Format is yyyymmddhh (yyyy = year, mm = month, dd = day, hh = hour (can only be 00,06,12, or 18))
#	driver_dset_type_str.................. Choose if Aeolus reprocessed data will be used: orig = not reprocessed (original data), for reprocessed files use baseline number (example: B10)
#
# Input: NetCDF
#	Driver dataset ....................... Daily Aeolus winds archived on S4, in yyyymmdd format
#	  -- Driver is independent dataset: All other datasets are collocated to Driver obs.
#	Dataset 1 ............................ 6-hourly Aircraft winds archived on S4 (converted from NCEP prepBUFR: aircft), in yyyymmddhh format
#	  -- Dataset 1 is dependent dataset: Dataset 1 is collocated to Driver obs.
#	Dataset 2 ............................ 6-hourly Atmospheric Motion Vectors (AMVs) archived on S4 (converted from NCEP prepBUFR: satwnd), in yyyymmddhh format
#	  -- Dataset 2 is dependent dataset: Dataset 2 is collocated to Driver obs.
#	Etc.
#
# Output: NetCDF
#	Files contain:
#		1. Indices of collocated observations. Indices pertain to arrays in each source data file.
#		2. Metadata listing collocation criteria, QC criteria (if applicable), paths to source data
#
#########################################################################################################
# Import python modules
#########################################################################################################

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
#os.system("module load miniconda/3.8-s4")

from os.path import exists

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
#from cartopy.util import add_cyclic_point

from netCDF4 import Dataset

from read_data import read_index_file

from read_data import read_aeolus_for_plotting
from read_data import read_aircraft_for_plotting
from read_data import read_amv_4th_int_for_plotting
from read_data import read_amv_ncep_for_plotting
from read_data import read_loon_for_plotting
from read_data import read_raobs_for_plotting
from read_data import read_raobs

from tools_analysis import compute_hlos
from tools_analysis import get_unique_values
from tools_analysis import count_unique_values
#from tools_analysis import superob_matches
from tools_analysis import matched_vars
from tools_analysis import prs_to_hgt

from tools_plotting import density_scatter
from tools_plotting import hist_diffs
from tools_plotting import map_locations_ce
from tools_plotting import map_locations_ortho
from tools_plotting import map_locations_ortho_rotate
from tools_plotting import map_3d_prof
from tools_plotting import scatter_matches

#########################################################################################################
print("***** BEGIN MAIN PROGRAM *****")

fill = -999.0

#=============================================
# Read raw user input from command line

aeolus_wind_type        = sys.argv[1]           # Aeolus wind type
dateIN                  = sys.argv[2]           # Date in question: YYYYMMDDHH
input_path		= sys.argv[3]		# Input directory for index files
input_file_suffix	= sys.argv[4]		# Input index file suffix
output_path             = sys.argv[5]           # Output (archive) directory
archive_parent		= sys.argv[6]		# Archive parent path: path where home archive directory is located
avgthin_choice		= int(sys.argv[7])		# Choice to super-ob, thin, or plot all matches

# temporary assignment
qi_choice = "NO_FC"

	# set AMV_4th_Int AMV type: Clear or Cloudy
if aeolus_wind_type=='RayClear':
  amv4th_type_str = 'Clear'
elif aeolus_wind_type=='MieCloud':
  amv4th_type_str = 'Cloud'
else:		# choose manually
  #amv4th_type_str = 'Clear'
  amv4th_type_str = 'Cloud'

print("CHECK USER INPUT ARGUMENTS:")
print("... Aeolus wind type = "+str(aeolus_wind_type))
print("... date IN = "+str(dateIN))
print("... input path = "+str(input_path))
print("... input index file suffix = "+str(input_file_suffix))
print("... output path = "+str(output_path))
print("... archive parent path = "+str(archive_parent))
print("... avg/thin choice = "+str(avgthin_choice))


#=============================================
# Full date arrays

mmARR 		= [ "01","02","03","04","05","06","07","08","09","10","11","12" ]
ddARRend 	= [ "31","28","31","30","31","30","31","31","30","31","30","31" ]
ddARR 		= [ "01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31" ]
hhARR		= [ "00","06","12","18" ]

#=============================================
# Find date before (B4) ... used for Aeolus input (driver dataset)

print("----- DATE = "+dateIN)

yyyy = dateIN[0:4]
mm   = dateIN[4:6]
dd   = dateIN[6:8]
hour = dateIN[8:10]

print("yyyy="+yyyy+",mm="+mm+",dd="+dd+",hour="+hour)

iyy = int(yyyy)
imm = int(mm)
idd = int(dd)

if dd == "01":
  if (mm == "01") | (mm == "02") | (mm == "04") | (mm == "06") | (mm == "08") | (mm == "09") | (mm == "11"):
    ddB4 = "31"
    if mm == "01":
      yyB4 = str(iyy-1)
      mmB4 = "12"
    else:
      yyB4 = yyyy
      immB4 = mmARR.index(mm)-1
      mmB4  = mmARR[immB4] 
  elif mm == "03":
    yyB4  = yyyy
    immB4 = mmARR.index(mm)-1
    mmB4  = mmARR[immB4]
    ddB4  = "28"
    if (iyy % 4) == 0:
      ddB4 = "29"
  elif (mm == "05") | (mm == "07") | (mm == "10") | (mm == "12"):
    yyB4  = yyyy
    immB4 = mmARR.index(mm)-1
    mmB4  = mmARR[immB4]
    ddB4  = "30"    
else:
  yyB4  = yyyy
  mmB4  = mm
  iddB4	= ddARR.index(dd)-1
  ddB4	= ddARR[iddB4]
 
 	#`````````````````````````````````````
yyyymmdd = yyyy+mm+dd		#Current date
dateB4   = yyB4+mmB4+ddB4	#Day before current date
 
print("Date current = "+yyyymmdd)
print("Date before  = "+dateB4) 

yyyymmddhh = yyyymmdd+hour

#=============================================
# Define Datasets

	#---------------------------------------------
	# Define and load INDEX FILE
	#---------------------------------------------
print("Define INDEX FILE")

		# read file
#tmp_dset2_path = '/atmos-nc-dataset/collocation/'+yyyy+'/'+mm+'/'+dd+'/'
#tmp_dset_path = '/data/users/klukens/collocation/longterm_anl/'
 
dset_path	      = input_path
dset_filename	      = 'index.'+yyyymmddhh+input_file_suffix
dset_file	      = dset_path+dset_filename

dset_exists = exists(dset_file)
if dset_exists==False:
  print("ERROR: index file "+dset_file+" does not exist!")
  sys.exit()
  
	      	# load dataset
data_hdl = Dataset(dset_file)

	   	# extract data
			# number of dependent datasets
ndset    = int(np.asarray( data_hdl.variables['ndset'] ))
			# collocation criteria
tim_max = np.asarray( data_hdl.variables['time_max'] )
prs_max = np.asarray( data_hdl.variables['pres_max'] )
hgt_max	= np.asarray( data_hdl.variables['hgt_max' ] )
dst_max = np.asarray( data_hdl.variables['dist_max'] )
			# DRIVER dataset attributes
drv_key = 'drv'
for attname in data_hdl.variables[drv_key].ncattrs():
  #att = getattr(data_hdl.variables['drv'],attname)
  #repr(nc_fid.variables[drv_key].getncattr(attname))
  att = data_hdl.variables[drv_key].getncattr(attname)
  print("attname = "+str(attname))
  print("... "+str(att))
  if attname == 'short_name':
    driver_name = att
  if attname == 'wind_type':
    aeolus_wind_type = att
  if attname == 'dataset_type':
    tatt = att
    if tatt.find("Reproc") != -1:
      tatt = tatt.strip("Reproc")
    else:
      tatt = "orig"
    aeolus_dset_type_str = tatt
    del tatt
  if attname == 'QC_flag':
    bool_drv_qc = False				# don't apply any QC. Need to extract the FULL dataset.
  if attname == 'QC_list':
    if ((driver_name.find('AMV') != -1) or (driver_name.find('amv') != -1) or (driver_name.find('BRZ') != -1) or (driver_name.find('EUM') != -1) or (driver_name.find('JMA') != -1) or (driver_name.find('KMA') != -1) or (driver_name.find('NOA') != -1) or (driver_name.find('NWC') != -1)):
      pct_drv = att
      pct_drv = pct_drv.strip(' %')
      #pct_drv = 0.0
  del att

print("driver_name = BRZ_"+amv4th_type_str)
del dset_path,dset_filename,dset_file,dset_exists 

idx_file_str = input_file_suffix
idx_file_vertstr = []
for i in range(ndset):
  print("i="+str(i))
  if i == 0:
    dset1_name,idx_drv_dset1,idx_dset1,DT_match1,GCD_match1,Vert_match1,Vert_str1 = read_index_file(input_path,yyyymmddhh,idx_file_str,int(i+1))
    idx_file_vertstr.append(Vert_str1)
  elif i == 1:
    dset2_name,idx_drv_dset2,idx_dset2,DT_match2,GCD_match2,Vert_match2,Vert_str2 = read_index_file(input_path,yyyymmddhh,idx_file_str,int(i+1))
    idx_file_vertstr.append(Vert_str2)
  elif i == 2:
    dset3_name,idx_drv_dset3,idx_dset3,DT_match3,GCD_match3,Vert_match3,Vert_str3 = read_index_file(input_path,yyyymmddhh,idx_file_str,int(i+1))
    idx_file_vertstr.append(Vert_str3)
  elif i == 3:
    dset4_name,idx_drv_dset4,idx_dset4,DT_match4,GCD_match4,Vert_match4,Vert_str4 = read_index_file(input_path,yyyymmddhh,idx_file_str,int(i+1))
    idx_file_vertstr.append(Vert_str4)
  elif i == 4:
    dset5_name,idx_drv_dset5,idx_dset5,DT_match5,GCD_match5,Vert_match5,Vert_str5 = read_index_file(input_path,yyyymmddhh,idx_file_str,int(i+1))
    idx_file_vertstr.append(Vert_str5)

	#---------------------------------------------
	# Define and load DRIVER dataset
	#---------------------------------------------
print("Define DRIVER dataset")

if driver_name=='Aeolus':
	# t_spd = HLOS wind velocity
	# t_dir = Aeolus azimuth angle
  drv_lat,drv_lon,drv_prs,drv_hgt,drv_yr,drv_mm,drv_dy,drv_hr,drv_mn,indexesD,qc_drv_list,drv_src,drv_err,drv_len,drv_spd,drv_dir = read_aeolus_for_plotting(archive_parent,yyyymmddhh,dateB4,aeolus_dset_type_str,aeolus_wind_type,bool_drv_qc)
elif driver_name=='Aircraft':
  drv_lat,drv_lon,drv_yr,drv_mm,drv_dy,drv_hr,drv_mn,drv_hgt,drv_prs,indexesD,qc_drv_list,drv_src,drv_spd,drv_dir = read_aircraft_for_plotting(archive_parent,yyyymmddhh,bool_drv_qc)
elif driver_name=='BRZ_'+amv4th_type_str or driver_name=='EUM_'+amv4th_type_str or driver_name=='JMA_'+amv4th_type_str or driver_name=='KMA_'+amv4th_type_str or driver_name=='NOA_'+amv4th_type_str or driver_name=='NWC_'+amv4th_type_str:
  drv_lat,drv_lon,drv_yr,drv_mm,drv_dy,drv_hr,drv_mn,drv_hgt,drv_prs,indexesD,qc_drv_list,drv_src,drv_spd,drv_dir = read_amv_4th_int_for_plotting(archive_parent,yyyymmddhh,bool_drv_qc,pct_drv,amv4th_type_str,driver_name)
elif driver_name=='AMV_NCEP':
  drv_lat,drv_lon,drv_yr,drv_mm,drv_dy,drv_hr,drv_mn,drv_hgt,drv_prs,indexesD,qc_drv_list,drv_src,drv_satname,drv_wcm,drv_ham,drv_spd,drv_dir = read_amv_ncep_for_plotting(archive_parent,yyyymmddhh,bool_drv_qc,pct_drv,qi_choice)
elif driver_name=='Loon':
  drv_lat,drv_lon,drv_yr,drv_mm,drv_dy,drv_hr,drv_mn,drv_hgt,drv_prs,indexesD,qc_drv_list,drv_src,drv_hgt,drv_azm,drv_elv,drv_spd,drv_dir = read_loon_for_plotting(archive_parent,yyyymmddhh,bool_drv_qc)
elif driver_name=='Radiosonde':
  drv_lat,drv_lon,drv_yr,drv_mm,drv_dy,drv_hr,drv_mn,drv_hgt,drv_prs,indexesD,qc_drv_list,drv_src,nsondes,nlevels,ngroups,drv_spd,drv_dir = read_raobs_for_plotting(archive_parent,yyyymmddhh,bool_drv_qc)
  #print("RAOBS spd PLOT: min/max = "+str(min(drv_spd))+","+str(max(drv_spd)))
else:
  print("ERROR: Driver dataset not defined!")
  sys.exit()

	# longitude check: set range to 0=360
for iP in range(np.size(drv_lon)):
  if drv_lon[iP] < 0.0:
    drv_lon[iP] = drv_lon[iP] + 360.

        # pressure check: set pressure to NaN if p<=0
for iP in range(np.size(drv_prs)):
  if drv_prs[iP] <= 0.0:
    drv_prs[iP] = np.nan

        # height check: set height to NaN if hgt<0
for iP in range(np.size(drv_hgt)):
  if drv_hgt[iP] < 0.0:
    drv_hgt[iP] = np.nan

print("... driver name: "+driver_name+" | drv size = "+str(np.size(drv_prs)))

	#---------------------------------------------
	# Define dependent datasets and append to each other to create total matching array for collocation
	#---------------------------------------------
print("Define DEPENDENT dataset(s)")

dependent_names = []
qc_dset_list = []
dset_src     = []

for i in range(ndset):
  print("i="+str(i))

  dset_key = "dset"+str(i+1)
  print("dset_key = "+dset_key)

	# get attributes for each dataset
  for attname in data_hdl.variables[dset_key].ncattrs():
    #att = getattr(data_hdl.variables[dset_key],attname)
    att = data_hdl.variables[dset_key].getncattr(attname)
    if attname == 'short_name':
      dependent_names.append(att)
      print("... dependent names = "+str(att))
    if attname == 'wind_type':
      aeolus_wind_type = att
    if attname == 'dataset_type':
      tatt = att
      if tatt.find("Reproc") != -1:
        tatt = tatt.strip("Reproc")
      else:
        tatt = "orig"
      aeolus_dset_type_str = tatt
      del tatt
    if attname == 'QC_flag':
      bool_dset_qc = False
    if attname == 'QC_list':
      if ((dependent_names[i].find('AMV') != -1) or (dependent_names[i].find('amv') != -1)):
        pct = att
        pct = pct.strip(' %')
    del att

  print("dependent_names = "+str(dependent_names))

  if idx_file_vertstr[i] != "NO_MATCHES":			# check if DEPENDENT dataset has any matches to DRIVER. If yes, proceed.
	# read dependent dataset
    if dependent_names[i]=='Aeolus':
  	# t_spd = HLOS wind velocity
	# t_dir = Aeolus azimuth angle
      t_lat,t_lon,t_prs,t_hgt,t_yr,t_mm,t_dy,t_hr,t_mn,t_indexes,t_qc_dset_list,t_src,t_err,t_len,t_spd,t_dir = read_aeolus_for_plotting(archive_parent,yyyymmddhh,dateB4,aeolus_dset_type_str,aeolus_wind_type,bool_dset_qc)
    elif dependent_names[i]=='Aircraft':
      t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,t_spd,t_dir = read_aircraft_for_plotting(archive_parent,yyyymmddhh,bool_dset_qc)
    elif dependent_names[i]=='AMV_4th_Int':
      t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,t_spd,t_dir = read_amv_4th_int_for_plotting(archive_parent,yyyymmddhh,bool_dset_qc,pct,amv4th_type_str)
    elif dependent_names[i]=='AMV_NCEP':
      t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,t_satname,t_wcm,t_ham,t_spd,t_dir = read_amv_ncep_for_plotting(archive_parent,yyyymmddhh,bool_dset_qc,pct,qi_choice)
    elif dependent_names[i]=='Loon':
      t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,t_azm,t_elv,t_spd,t_dir = read_loon_for_plotting(archive_parent,yyyymmddhh,bool_dset_qc)
    elif dependent_names[i]=='Radiosonde':
	# RADIOSONDE data are 2D (except for lat, lon, dates/times): [nsondes, nlevels]
	# 	nsondes = number of sondes
	# 	nlevels = number of vertical levels per sonde
	# nsondes, nlevels are used to find collocated indices that pertain to original data file
      t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,nsondes,nlevels,ngroups,t_spd,t_dir = read_raobs_for_plotting(archive_parent,yyyymmddhh,bool_dset_qc)
    else:
      print("ERROR: Dependent dataset "+str(i)+" ("+dependent_names[i]+") cannot be read by this program! Please add function to read_input_datasets module and try again.")
      sys.exit()

	# longitude check: set range to 0=360
    for iP in range(np.size(t_lon)):
      if t_lon[iP] < 0.0:
        t_lon[iP] = t_lon[iP] + 360.

	# pressure check: set pressure to NaN if p<=0
    for iP in range(np.size(t_prs)):
      if t_prs[iP] <= 0.0:
        t_prs[iP] = np.nan

	# height check: set height to NaN if hgt<0
    for iP in range(np.size(t_hgt)):
      if t_hgt[iP] < 0.0:
        t_hgt[iP] = np.nan
      
    print("size drv: "+str(np.size(drv_prs))+" | size dset: "+str(np.size(t_prs)))

	# extract matched obs for plotting
    if i==0:
      drv1_lat,drv1_lon,drv1_prs,drv1_hgt,drv1_spd,drv1_dir,dset1_lat,dset1_lon,dset1_prs,dset1_hgt,dset1_spd,dset1_dir = matched_vars(avgthin_choice,driver_name,dependent_names[i],idx_drv_dset1,idx_dset1,drv_lat,drv_lon,drv_prs,drv_hgt,drv_spd,drv_dir,t_lat,t_lon,t_prs,t_hgt,t_spd,t_dir)      

    if i==1:
      drv2_lat,drv2_lon,drv2_prs,drv2_hgt,drv2_spd,drv2_dir,dset2_lat,dset2_lon,dset2_prs,dset2_hgt,dset2_spd,dset2_dir = matched_vars(avgthin_choice,driver_name,dependent_names[i],idx_drv_dset2,idx_dset2,drv_lat,drv_lon,drv_prs,drv_hgt,drv_spd,drv_dir,t_lat,t_lon,t_prs,t_hgt,t_spd,t_dir)	 

    if i==2:
      drv3_lat,drv3_lon,drv3_prs,drv3_hgt,drv3_spd,drv3_dir,dset3_lat,dset3_lon,dset3_prs,dset3_hgt,dset3_spd,dset3_dir = matched_vars(avgthin_choice,driver_name,dependent_names[i],idx_drv_dset3,idx_dset3,drv_lat,drv_lon,drv_prs,drv_hgt,drv_spd,drv_dir,t_lat,t_lon,t_prs,t_hgt,t_spd,t_dir)	 

    if i==3:
      drv4_lat,drv4_lon,drv4_prs,drv4_hgt,drv4_spd,drv4_dir,dset4_lat,dset4_lon,dset4_prs,dset4_hgt,dset4_spd,dset4_dir = matched_vars(avgthin_choice,driver_name,dependent_names[i],idx_drv_dset4,idx_dset4,drv_lat,drv_lon,drv_prs,drv_hgt,drv_spd,drv_dir,t_lat,t_lon,t_prs,t_hgt,t_spd,t_dir)	 

    if i==4:
      drv5_lat,drv5_lon,drv5_prs,drv5_hgt,drv5_spd,drv5_dir,dset5_lat,dset5_lon,dset5_prs,dset5_hgt,dset5_spd,dset5_dir = matched_vars(avgthin_choice,driver_name,dependent_names[i],idx_drv_dset5,idx_dset5,drv_lat,drv_lon,drv_prs,drv_hgt,drv_spd,drv_dir,t_lat,t_lon,t_prs,t_hgt,t_spd,t_dir)

	# substring for plot filenames
    if avgthin_choice==-1:
      print("PLOT ALL MATCHES")
      avgthin_str = ".AllMatches"
    elif avgthin_choice==0:
      print("SUPER-OB (AVERAGE) MATCHES")
      avgthin_str = ".SuperOb"
    elif avgthin_choice==1:
      print("THIN MATCHES: TBA")
      avgthin_str = ".Thin"
    print("avgthin_str = "+avgthin_str)

data_hdl.close()

#==================================================================
# PLOTS
#==================================================================

dcolors = ["red","limegreen","blue","magenta","orange"]			# colors per DEPENDENT dataset
imark   = ["^","s","p","d","*"]						# markers per DEPENDENT dataset

	#----------------------------------------------------------
	# Maps: Locations of matched obs

alphaval = 0.25		# transparency factor (between 0 and 1, with 1=opaque)

		# select region to plot (CE map)
		#	0 = global
		#	1 = region automatically limited to the matched locations
region_flag = 0
#region_flag = 1

x_name = driver_name
y_name = ""	
if driver_name == 'Aeolus':
  x_name += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

jDs = []
jss = []
jDx = []
jxx = []
jDy = []
jyy = []
ja_tname = []
jmcolors = []
jmmarks  = []

for i in range(ndset):
  if idx_file_vertstr[i] != "NO_MATCHES":			# check if DEPENDENT dataset has any matches to DRIVER. If yes, proceed.
    if i == 0:
      Dts   = drv1_spd
      ts    = dset1_spd
      Dtx   = drv1_lon
      Dty   = drv1_lat
      tx    = dset1_lon
      ty    = dset1_lat
      tname = dset1_name
    if i == 1:
      Dts   = drv2_spd
      ts    = dset2_spd
      Dtx   = drv2_lon
      Dty   = drv2_lat
      tx    = dset2_lon
      ty    = dset2_lat
      tname = dset2_name
    if i == 2:
      Dts   = drv3_spd
      ts    = dset3_spd
      Dtx   = drv3_lon
      Dty   = drv3_lat
      tx    = dset3_lon
      ty    = dset3_lat
      tname = dset3_name
    if i == 3:
      Dts   = drv4_spd
      ts    = dset4_spd
      Dtx   = drv4_lon
      Dty   = drv4_lat
      tx    = dset4_lon
      ty    = dset4_lat
      tname = dset4_name
    if i == 4:
      Dts   = drv5_spd
      ts    = dset5_spd
      Dtx   = drv5_lon
      Dty   = drv5_lat
      tx    = dset5_lon
      ty    = dset5_lat
      tname = dset5_name
    
	# append
    jDs.append(Dts)
    jss.append(ts)
    jDx.append(Dtx)
    jxx.append(tx)
    jDy.append(Dty)
    jyy.append(ty)
    ja_tname.append(tname)
    jmcolors.append(dcolors[i])
    jmmarks.append(imark[i])

    	# y_name string for plot filename
    y_name += "_"+tname
    if tname == 'Aeolus':
      y_name += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

    del Dts,ts,Dtx,tx,Dty,ty,tname
    
shape_D = np.shape(jss)
    
    	# sort datasets by size. Plot largest first ... to smallest last.
aDs	 = np.asarray(jDs     )
ass	 = np.asarray(jss     )
aDx	 = np.asarray(jDx     )
axx	 = np.asarray(jxx     )
aDy	 = np.asarray(jDy     )
ayy	 = np.asarray(jyy     )
aa_tname = np.asarray(ja_tname)
amcolors = np.asarray(jmcolors)
ammarks  = np.asarray(jmmarks )

sizes = np.nan * np.ones(shape_D[0],dtype=int)
for i in range(shape_D[0]):
  sizes[i] = np.size(aDs[i])

print("MAP sizes = "+str(sizes))

sort_size = np.sort(sizes)[::-1]	      	# [::-1] = sort in descending order. For ascending, comment out [::-1]
idx_sort = []
for i in range(shape_D[0]):  	      		# loop to go thru sort_size
  for j in range(shape_D[0]):	      		# loop to go thru size(aDs)
    if np.size(aDs[j])==sort_size[i]:
      idx_sort.append(int(j))
print(type(idx_sort))
print(idx_sort)

Ds      = aDs[idx_sort]
ss      = ass[idx_sort]
Dx      = aDx[idx_sort]
xx      = axx[idx_sort]
Dy      = aDy[idx_sort]
yy      = ayy[idx_sort]
a_tname = aa_tname[idx_sort]
mcolors = amcolors[idx_sort]
mmarks  = ammarks[idx_sort]

print("MAP sort_size = "+str(sort_size)+" | idx_sort = "+str(idx_sort))
#print("aDs = "+str(aDs))
#print("Ds  = "+str(Ds))
#sys.exit()

del sort_size,idx_sort,sizes,shape_D
del jDs,jss,jDx,jxx,jDy,jyy,ja_tname,jmcolors,jmmarks
del aDs,ass,aDx,axx,aDy,ayy,aa_tname,amcolors,ammarks

	#```````````````````````````````````````````
	# Cylindrical Equidistant Map
	
marksize=15

map_locations_ce(Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,region_flag,marksize,mmarks,alphaval,mcolors,dateIN)

		# plot filename
outname = output_path+"MAP_CE.Match_Locations."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
		# save plot
plt.savefig(outname)

del outname

	#```````````````````````````````````````````
	# Orthgraphic Maps

#marksize=15

	# North Pole
		# lat/lon representing center point of projection (ORTHO maps)
                #       lon = 0.0   --> Polar Stereographic Projection
                #       lat = 90.0  --> North Pole
                #       lat = 0.0   --> Equator
                #       lat = -90.0 --> South Pole
central_lon = 0.0
central_lat = 90.0

map_locations_ortho(Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,marksize,mmarks,alphaval,mcolors,dateIN,central_lon,central_lat)

		# plot filename
outname = output_path+"MAP_Ortho.Match_Locations_NorthPole."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
		# save plot
plt.savefig(outname)

del outname

	# South Pole
                # lat/lon representing center point of projection (ORTHO maps)
                #       lon = 0.0   --> Polar Stereographic Projection
                #       lat = 90.0  --> North Pole
                #       lat = 0.0   --> Equator
                #       lat = -90.0 --> South Pole
central_lon = 0.0
central_lat = -90.0

map_locations_ortho(Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,marksize,mmarks,alphaval,mcolors,dateIN,central_lon,central_lat)

                # plot filename
outname = output_path+"MAP_Ortho.Match_Locations_SouthPole."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
                # save plot
plt.savefig(outname)

del outname

	# Prime Meridian center
                # lat/lon representing center point of projection (ORTHO maps)
                #       lon = 0.0   --> Polar Stereographic Projection
                #       lat = 90.0  --> North Pole
                #       lat = 0.0   --> Equator
                #       lat = -90.0 --> South Pole
#central_lon = 0.0
#central_lat = 0.0
#
#map_locations_ortho(Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,marksize,mmarks,alphaval,mcolors,dateIN,central_lon,central_lat)
#
#                # plot filename
#outname = output_path+"MAP_Ortho.Match_Locations_Lon0center."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
#                # save plot
#plt.savefig(outname)
#
#del outname

	# 180 Longitude center
                # lat/lon representing center point of projection (ORTHO maps)
                #       lon = 0.0   --> Polar Stereographic Projection
                #       lat = 90.0  --> North Pole
                #       lat = 0.0   --> Equator
                #       lat = -90.0 --> South Pole
#central_lon = 180.0
#central_lat = 0.0
#
#map_locations_ortho(Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,marksize,mmarks,alphaval,mcolors,dateIN,central_lon,central_lat)
#
#                # plot filename
#outname = output_path+"MAP_Ortho.Match_Locations_Lon180center."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
#                # save plot
#plt.savefig(outname)
#
#del outname

	# Eastern Hemisphere
                # lat/lon representing center point of projection (ORTHO maps)
                #       lon = 0.0   --> Polar Stereographic Projection
                #       lat = 90.0  --> North Pole
                #       lat = 0.0   --> Equator
                #       lat = -90.0 --> South Pole
#central_lon = 90.0
#central_lat = 0.0
#
#map_locations_ortho(Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,marksize,mmarks,alphaval,mcolors,dateIN,central_lon,central_lat)
#
#                # plot filename
#outname = output_path+"MAP_Ortho.Match_Locations_EastHem."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
#                # save plot
#plt.savefig(outname)
#
#del outname

	# Western Hemisphere
                # lat/lon representing center point of projection (ORTHO maps)
                #       lon = 0.0   --> Polar Stereographic Projection
                #       lat = 90.0  --> North Pole
                #       lat = 0.0   --> Equator
                #       lat = -90.0 --> South Pole
central_lon = 270.0
central_lat = 0.0

map_locations_ortho(Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,marksize,mmarks,alphaval,mcolors,dateIN,central_lon,central_lat)

                # plot filename
outname = output_path+"MAP_Ortho.Match_Locations_WestHem."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
                # save plot
plt.savefig(outname)

del outname

	#```````````````````````````````````````````
	# Rotating Orthgraphic Maps

	# Rotate around Equator
		# lat representing latitude around which projection rotates (ORTHO maps)
central_lat = 0

		# plot filename
outname = output_path+"MAP_Rotate.Match_Locations_centerlat"+str(central_lat)+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str

map_locations_ortho_rotate(output_path,outname,Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,marksize,mmarks,alphaval,mcolors,dateIN,central_lat)

del outname

	#```````````````````````````````````````````

del x_name,y_name,Ds,ss,Dx,xx,Dy,yy,a_tname,mcolors,mmarks

	#----------------------------------------------------------
	# Histograms of Collocation Difference
	
alphaval = 1.0		# transparency factor (between 0 and 1, with 1=opaque)

x_name  = driver_name
y_name  = ""
py_name = ""
zy_name = ""
if driver_name == 'Aeolus':
  x_name += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

latD = []
lonD = []

a_tmatch    = []
a_pmatch    = []
a_zmatch    = []
a_distmatch = []

a_tname     = []
a_ptname    = []
a_ztname    = []

a_tcolors    = []
a_pcolors    = []
a_zcolors    = []
a_distcolors = []

for i in range(ndset):
  pmatch = []
  zmatch = []
  if idx_file_vertstr[i] != "NO_MATCHES":			# check if DEPENDENT dataset has any matches to DRIVER. If yes, proceed.
    if i == 0:
      tmatch    = DT_match1
      distmatch = GCD_match1
      tlat	= drv1_lat
      tlon	= drv1_lon
      if Vert_str1=="Pressure":
        pmatch  = Vert_match1
      else:
        zmatch  = Vert_match1
      tname     = dset1_name
    if i == 1:
      tmatch    = DT_match2
      tlat	= drv2_lat
      tlon	= drv2_lon
      if Vert_str2=="Pressure":
        pmatch  = Vert_match2
      else:
        zmatch  = Vert_match2
      distmatch = GCD_match2
      tname     = dset2_name
    if i == 2:
      tmatch    = DT_match3
      tlat	= drv3_lat
      tlon	= drv3_lon
      if Vert_str3=="Pressure":
        pmatch  = Vert_match3
      else:
        zmatch  = Vert_match3
      distmatch = GCD_match3
      tname     = dset3_name
    if i == 3:
      tmatch    = DT_match4
      tlat	= drv4_lat
      tlon	= drv4_lon
      if Vert_str4=="Pressure":
        pmatch  = Vert_match4
      else:
        zmatch  = Vert_match4
      distmatch = GCD_match4
      tname     = dset4_name
    if i == 4:
      tmatch    = DT_match5
      tlat	= drv5_lat
      tlon	= drv5_lon
      if Vert_str5=="Pressure":
        pmatch  = Vert_match5
      else:
        zmatch  = Vert_match5
      distmatch = GCD_match5
      tname     = dset5_name
 
  	# y_name string for plot filename
    y_name += "_"+tname
    if tname == 'Aeolus':
      y_name += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

    	# arrays for histogram
    a_tname.append(tname)

    latD.append(tlat)
    lonD.append(tlon)
    
    a_tmatch.append(tmatch)
    a_distmatch.append(distmatch)

    a_tcolors.append(dcolors[i])
    a_distcolors.append(dcolors[i])

    if np.size(pmatch)>0:
      a_pmatch.append(pmatch)
      a_pcolors.append(dcolors[i])
      a_ptname.append(tname)
      py_name += "_"+tname
      del pmatch
    if np.size(zmatch)>0:
      zzmatch = zmatch*1000.0		# convert km to m
      a_zmatch.append(zzmatch)
      a_zcolors.append(dcolors[i])
      a_ztname.append(tname)
      zy_name += "_"+tname
      del zmatch,zzmatch
  
    del tmatch,distmatch,tname

	# call hist_diffs from tools_plotting.py
		# TIME
match_str = "Time"
units = "min"
hist_diffs(latD, lonD, a_tmatch, x_name, a_tname, alphaval, match_str, a_tcolors, units)
outname = output_path+"HIST.Match_Time_Diff."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
plt.savefig(outname)
del outname,a_tmatch,units
		# DISTANCE
match_str = "Distance"
units = "km"
hist_diffs(latD, lonD, a_distmatch, x_name, a_tname, alphaval, match_str, a_distcolors, units)
outname = output_path+"HIST.Match_Distance."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
plt.savefig(outname)
del outname,a_distmatch,units
		# VERTICAL
			# Pressure
if np.size(a_pmatch)>0:
  match_str = "Pressure"
  units = "hPa"
  hist_diffs(latD, lonD, a_pmatch, x_name, a_ptname, alphaval, match_str, a_pcolors, units)
  outname = output_path+"HIST.Match_Pressure_Diff."+str(dateIN)+".x_"+str(x_name)+".y"+str(py_name)+avgthin_str+".png"
  plt.savefig(outname)
  del outname,a_pmatch,units
			# Height
if np.size(a_zmatch)>0:
  match_str = "Height"
  units = "m"
  hist_diffs(latD, lonD, a_zmatch, x_name, a_ztname, alphaval, match_str, a_zcolors, units)
  outname = output_path+"HIST.Match_Height_Diff."+str(dateIN)+".x_"+str(x_name)+".y"+str(zy_name)+avgthin_str+".png"
  plt.savefig(outname)
  del outname,a_zmatch,units

del y_name

	#----------------------------------------------------------
	# Scatterplots of Matched Obs
	
marksize = 50
alphaval = 0.25

x_name = driver_name
y_name_hlos = ""
y_name_wspd = ""
y_name_pres = ""
y_name_hgt  = ""
if driver_name == 'Aeolus':
  x_name += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

latD  = []
lonD  = []

hlosD = []
hlos  = []
wspdD = []
wspd  = []
presD = []
pres  = []
hgtD  = []
hgt   = []

a_tname_hlos = []
a_tname_wspd = []
a_tname_pres = []
a_tname_hgt  = []
a_color_hlos = []
a_color_wspd = []
a_color_pres = []
a_color_hgt  = []
a_mark_hlos  = []
a_mark_wspd  = []
a_mark_pres  = []
a_mark_hgt   = []

Vert_str = []

for i in range(ndset):
  if idx_file_vertstr[i] != "NO_MATCHES":			# check if DEPENDENT dataset has any matches to DRIVER. If yes, proceed.
    if i == 0:
      tname = dset1_name
      txspd = drv1_spd
      txlat = drv1_lat
      txlon = drv1_lon
      tyspd = dset1_spd
      Vert_str = Vert_str1
      if Vert_str1=="Pressure":
        txp = drv1_prs
        typ = dset1_prs
      else:
        if np.isnan(drv1_hgt).all():
          txz = prs_to_hgt(drv1_prs)
        else:
          txz = drv1_hgt
        tyz = dset1_hgt
    if i == 1:
      tname = dset2_name
      txspd = drv2_spd
      txlat = drv2_lat
      txlon = drv2_lon
      tyspd = dset2_spd
      Vert_str = Vert_str2
      if Vert_str2=="Pressure":
        txp = drv2_prs
        typ = dset2_prs
      else:
        if np.isnan(drv2_hgt).all():
          txz = prs_to_hgt(drv2_prs)
        else:
          txz = drv2_hgt
        tyz = dset2_hgt
    if i == 2:
      tname = dset3_name
      txspd = drv3_spd
      txlat = drv3_lat
      txlon = drv3_lon
      tyspd = dset3_spd
      Vert_str = Vert_str3
      if Vert_str3=="Pressure":
        txp = drv3_prs
        typ = dset3_prs
      else:
        if np.isnan(drv3_hgt).all():
          txz = prs_to_hgt(drv3_prs)
        else:
          txz = drv3_hgt
        tyz = dset3_hgt
    if i == 3:
      tname = dset4_name
      txspd = drv4_spd
      txlat = drv4_lat
      txlon = drv4_lon
      tyspd = dset4_spd
      Vert_str = Vert_str4
      if Vert_str4=="Pressure":
        txp = drv4_prs
        typ = dset4_prs
      else:
        if np.isnan(drv4_hgt).all():
          txz = prs_to_hgt(drv4_prs)
        else:
          txz = drv4_hgt
        tyz = dset4_hgt
    if i == 4:
      tname = dset5_name
      txspd = drv5_spd
      txlat = drv5_lat
      txlon = drv5_lon
      tyspd = dset5_spd
      Vert_str = Vert_str5
      if Vert_str5=="Pressure":
        txp = drv5_prs
        typ = dset5_prs
      else:
        if np.isnan(drv5_hgt).all():
          txz = prs_to_hgt(drv5_prs)
        else:
          txz = drv5_hgt
        tyz = dset5_hgt

    if tname == 'Aeolus':
      tname += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

    latD.append(txlat)
    lonD.append(txlon)

	# append wind speeds
    if driver_name == 'Aeolus' or tname.find('Aeolus') != -1:	# find() != -1 --> string contains substring!
    	# HLOS wind
      hlosD.append(txspd)
      hlos.append(tyspd)
      a_tname_hlos.append(tname)
      a_color_hlos.append(dcolors[i])
      a_mark_hlos.append(imark[i])
      print("hlos color = "+str(a_color_hlos))
        # y_name string for plot filename
      y_name_hlos += "_"+tname
    if driver_name != 'Aeolus' and tname.find('Aeolus') == -1:	# find() == -1 --> string does NOT contain substring!
  	# non-HLOS wind
      wspdD.append(txspd)
      wspd.append(tyspd)
      a_tname_wspd.append(tname)
      a_color_wspd.append(dcolors[i])
      a_mark_wspd.append(imark[i])
      print("wspd color = "+str(a_color_wspd))
  	# y_name string for plot filename
      y_name_wspd += "_"+tname
	
	# append pressures or heights
    if Vert_str=="Pressure":
      presD.append(txp)
      pres.append(typ)
      a_tname_pres.append(tname)
      a_color_pres.append(dcolors[i])
      a_mark_pres.append(imark[i])
      #print("wspd color = "+str(a_color_wspd))
  	# y_name string for plot filename
      y_name_pres += "_"+tname
      del txp,typ		
    if Vert_str!="Pressure":
      hgtD.append(txz)
      hgt.append(tyz)
      a_tname_hgt.append(tname)
      a_color_hgt.append(dcolors[i])
      a_mark_hgt.append(imark[i])
      #print("wspd color = "+str(a_color_wspd))
  	# y_name string for plot filename
      y_name_hgt += "_"+tname
      del txz,tyz		

    print("SCATTER: tname = "+tname+" | color = "+str(dcolors[i])+" | mark = "+str(imark[i]))

    del txspd,tyspd,tname

shape_hlos = np.shape(hlosD)
shape_wspd = np.shape(wspdD)
shape_pres = np.shape(presD)
shape_hgt  = np.shape(hgtD )

	#```````````````````````````````````````````
	# plot scatterplots
		#```````````````````````````````````````````
		# HLOS WIND
if shape_hlos[0]>0:
  #print("shape D: "+str(np.shape(hlosD)))

  aDs      = np.asarray(hlosD       )
  ass	   = np.asarray(hlos        )
  aa_tname = np.asarray(a_tname_hlos)
  amcolors = np.asarray(a_color_hlos)
  ammarks  = np.asarray(a_mark_hlos )

		# sort datasets by size. Plot largest first ... to smallest last.
  sizes = np.nan * np.ones(shape_hlos[0],dtype=int)
  for i in range(shape_hlos[0]):
    sizes[i] = np.size(aDs[i])

  sort_size = np.sort(sizes)[::-1]		# [::-1] = sort in descending order. For ascending, comment out [::-1]
  idx_sort = []
  for i in range(shape_hlos[0]):		# loop to go thru sort_size
    for j in range(shape_hlos[0]):		# loop to go thru size(aDs)
      if np.size(aDs[j])==sort_size[i]:
        idx_sort.append(int(j))

  Ds      = aDs[idx_sort]
  ss      = ass[idx_sort]
  a_tname = aa_tname[idx_sort]
  mcolors = amcolors[idx_sort]
  mmarks  = ammarks[idx_sort]
  del sort_size,idx_sort,sizes
  del hlosD,hlos,a_tname_hlos,a_color_hlos,a_mark_hlos
  del aDs,ass,aa_tname,amcolors,ammarks

		# plot
  match_str = "HLOS"
  units = "m/s"
  scatter_matches(latD,lonD,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks)
  outname = output_path+"SCATTER.HLOS_Velocity."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hlos)+avgthin_str+".png"
  plt.savefig(outname)	
  del outname,units,Ds,ss,a_tname,mcolors,mmarks
  
  		#WIND SPEED (not HLOS)
if shape_wspd[0]>0:
  #print("shape D: "+str(np.shape(wspdD)))
		
  aDs      = np.asarray(wspdD       )
  ass	   = np.asarray(wspd        )
  aa_tname = np.asarray(a_tname_wspd)
  amcolors = np.asarray(a_color_wspd)
  ammarks  = np.asarray(a_mark_wspd )

		# sort datasets by size. Plot largest first ... to smallest last.
  sizes = np.nan * np.ones(shape_wspd[0],dtype=int)
  for i in range(shape_wspd[0]):
    sizes[i] = np.size(aDs[i])

  sort_size = np.sort(sizes)[::-1]		# [::-1] = sort in descending order. For ascending, comment out [::-1]
  idx_sort = []
  for i in range(shape_wspd[0]):		# loop to go thru sort_size
    for j in range(shape_wspd[0]):		# loop to go thru size(aDs)
      if np.size(aDs[j])==sort_size[i]:
        idx_sort.append(int(j))

  Ds      = aDs[idx_sort]
  ss      = ass[idx_sort]
  a_tname = aa_tname[idx_sort]
  mcolors = amcolors[idx_sort]
  mmarks  = ammarks[idx_sort]
  del sort_size,idx_sort,sizes
  del wspdD,wspd,a_tname_wspd,a_color_wspd,a_mark_wspd
  del aDs,ass,aa_tname,amcolors,ammarks

		# plot
  match_str = "Wind Speed"
  units = "m/s"
  scatter_matches(latD,lonD,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks)
  outname = output_path+"SCATTER.WindSpeed."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_wspd)+avgthin_str+".png"
  plt.savefig(outname)	
  del outname,units,Ds,ss,a_tname,mcolors,mmarks
  
  		#```````````````````````````````````````````
		# PRESSURE
if shape_pres[0]>0:
  #print("shape D: "+str(np.shape(presD)))
		
  aDs      = np.asarray(presD       )
  ass	   = np.asarray(pres        )
  aa_tname = np.asarray(a_tname_pres)
  amcolors = np.asarray(a_color_pres)
  ammarks  = np.asarray(a_mark_pres )

		# sort datasets by size. Plot largest first ... to smallest last.
  sizes = np.nan * np.ones(shape_pres[0],dtype=int)
  for i in range(shape_pres[0]):
    sizes[i] = np.size(aDs[i])

  sort_size = np.sort(sizes)[::-1]		# [::-1] = sort in descending order. For ascending, comment out [::-1]
  idx_sort = []
  for i in range(shape_pres[0]):		# loop to go thru sort_size
    for j in range(shape_pres[0]):		# loop to go thru size(aDs)
      if np.size(aDs[j])==sort_size[i]:
        idx_sort.append(int(j))

  Ds      = aDs[idx_sort]
  ss      = ass[idx_sort]
  a_tname = aa_tname[idx_sort]
  mcolors = amcolors[idx_sort]
  mmarks  = ammarks[idx_sort]
  del sort_size,idx_sort,sizes
  del presD,pres,a_tname_pres,a_color_pres,a_mark_pres
  del aDs,ass,aa_tname,amcolors,ammarks

		# plot
  match_str = "Pressure"
  units = "hPa"
  scatter_matches(latD,lonD,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks)
  outname = output_path+"SCATTER.Pressure."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_pres)+avgthin_str+".png"
  plt.savefig(outname)	
  del outname,units,Ds,ss,a_tname,mcolors,mmarks
  
  		#```````````````````````````````````````````
		# HEIGHT
if shape_hgt[0]>0:
  print("shape hgt D: "+str(np.shape(hgtD)))
  print("size hgt: D="+str(np.size(hgtD))+" | "+str(np.size(hgt)))
		
  aDs      = np.asarray(hgtD       )
  ass	   = np.asarray(hgt        )
  aa_tname = np.asarray(a_tname_hgt)
  amcolors = np.asarray(a_color_hgt)
  ammarks  = np.asarray(a_mark_hgt )

		# sort datasets by size. Plot largest first ... to smallest last.
  sizes = np.nan * np.ones(shape_hgt[0],dtype=int)
  for i in range(shape_hgt[0]):
    sizes[i] = np.size(aDs[i])

  sort_size = np.sort(sizes)[::-1]		# [::-1] = sort in descending order. For ascending, comment out [::-1]
  idx_sort = []
  for i in range(shape_hgt[0]):		# loop to go thru sort_size
    for j in range(shape_hgt[0]):		# loop to go thru size(aDs)
      if np.size(aDs[j])==sort_size[i]:
        idx_sort.append(int(j))

  Ds      = aDs[idx_sort]
  ss      = ass[idx_sort]
  a_tname = aa_tname[idx_sort]
  mcolors = amcolors[idx_sort]
  mmarks  = ammarks[idx_sort]
  del sort_size,idx_sort,sizes
  del hgtD,hgt,a_tname_hgt,a_color_hgt,a_mark_hgt
  del aDs,ass,aa_tname,amcolors,ammarks

		# plot
  match_str = "Height"
  units = "km"
  scatter_matches(latD,lonD,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks)
  outname = output_path+"SCATTER.Height."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hgt)+avgthin_str+".png"
  plt.savefig(outname)	
  del outname,units,Ds,ss,a_tname,mcolors,mmarks
  
  	#```````````````````````````````````````````

del x_name,y_name_hlos,y_name_wspd,y_name_pres,y_name_hgt,shape_hlos,shape_wspd,shape_pres,shape_hgt,latD,lonD

	#----------------------------------------------------------
	# Density Scatterplots: Wind Speed

x_name = driver_name
if driver_name == 'Aeolus':
  x_name += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

	#DRIVER - WIND SPEED
zzlabel  = Vert_str
tlat = []
tlon = []
tspd = []
tp   = []
tz   = []
thloslat = []
thloslon = []
thlosspd = []
thlosp   = []
thlosz   = []
for i in range(ndset):
  if i==0:
    Vert_str = Vert_str1
    try: drv1_lat
    except NameError: drv1_lat=None
    if drv1_lat is not None:
      if driver_name=="Aeolus" or dependent_names[i]=="Aeolus":
        thloslat = np.append(thloslat,drv1_lat,axis=0)
        thloslon = np.append(thloslon,drv1_lon,axis=0)
        thlosspd = np.append(thlosspd,drv1_spd,axis=0)
        if Vert_str1=="Pressure":
          thlosp = np.append(thlosp,drv1_prs,axis=0)
        else:
          thlosz = np.append(thlosz,drv1_hgt,axis=0)
      else:
        tlat = np.append(tlat,drv1_lat,axis=0)
        tlon = np.append(tlon,drv1_lon,axis=0)
        tspd = np.append(tspd,drv1_spd,axis=0)
        if Vert_str1=="Pressure":
          tp = np.append(tp,drv1_prs,axis=0)
        else:
          tz = np.append(tz,drv1_hgt,axis=0)
  if i==1:
    Vert_str = Vert_str2
    try: drv2_lat
    except NameError: drv2_lat=None
    if drv2_lat is not None:
      if driver_name=="Aeolus" or dependent_names[i]=="Aeolus":
        thloslat = np.append(thloslat,drv2_lat,axis=0)
        thloslon = np.append(thloslon,drv2_lon,axis=0)
        thlosspd = np.append(thlosspd,drv2_spd,axis=0)
        if Vert_str2=="Pressure":
          thlosp = np.append(thlosp,drv2_prs,axis=0)
        else:
          thlosz = np.append(thlosz,drv2_hgt,axis=0)
      else:
        tlat = np.append(tlat,drv2_lat,axis=0)
        tlon = np.append(tlon,drv2_lon,axis=0)
        tspd = np.append(tspd,drv2_spd,axis=0)
        if Vert_str2=="Pressure":
          tp = np.append(tp,drv2_prs,axis=0)
        else:
          tz = np.append(tz,drv2_hgt,axis=0)
  if i==2:
    Vert_str = Vert_str3
    try: drv3_lat
    except NameError: drv3_lat=None
    if drv3_lat is not None:
      if driver_name=="Aeolus" or dependent_names[i]=="Aeolus":
        thloslat = np.append(thloslat,drv3_lat,axis=0)
        thloslon = np.append(thloslon,drv3_lon,axis=0)
        thlosspd = np.append(thlosspd,drv3_spd,axis=0)
        if Vert_str3=="Pressure":
          thlosp = np.append(thlosp,drv3_prs,axis=0)
        else:
          thlosz = np.append(thlosz,drv3_hgt,axis=0)
      else:
        tlat = np.append(tlat,drv3_lat,axis=0)
        tlon = np.append(tlon,drv3_lon,axis=0)
        tspd = np.append(tspd,drv3_spd,axis=0)
        if Vert_str3=="Pressure":
          tp = np.append(tp,drv3_prs,axis=0)
        else:
          tz = np.append(tz,drv3_hgt,axis=0)
  if i==3:
    Vert_str = Vert_str4
    try: drv4_lat
    except NameError: drv4_lat=None
    if drv4_lat is not None:
      if driver_name=="Aeolus" or dependent_names[i]=="Aeolus":
        thloslat = np.append(thloslat,drv4_lat,axis=0)
        thloslon = np.append(thloslon,drv4_lon,axis=0)
        thlosspd = np.append(thlosspd,drv4_spd,axis=0)
        if Vert_str4=="Pressure":
          thlosp = np.append(thlosp,drv4_prs,axis=0)
        else:
          thlosz = np.append(thlosz,drv4_hgt,axis=0)
      else:
        tlat = np.append(tlat,drv4_lat,axis=0)
        tlon = np.append(tlon,drv4_lon,axis=0)
        tspd = np.append(tspd,drv4_spd,axis=0)
        if Vert_str4=="Pressure":
          tp = np.append(tp,drv4_prs,axis=0)
        else:
          tz = np.append(tz,drv4_hgt,axis=0)
  if i==4:
    Vert_str = Vert_str5
    try: drv5_lat
    except NameError: drv5_lat=None
    if drv5_lat is not None:
      if driver_name=="Aeolus" or dependent_names[i]=="Aeolus":
        thloslat = np.append(thloslat,drv5_lat,axis=0)
        thloslon = np.append(thloslon,drv5_lon,axis=0)
        thlosspd = np.append(thlosspd,drv5_spd,axis=0)
        if Vert_str5=="Pressure":
          thlosp = np.append(thlosp,drv5_prs,axis=0)
        else:
          thlosz = np.append(thlosz,drv5_hgt,axis=0)
      else:
        tlat = np.append(tlat,drv5_lat,axis=0)
        tlon = np.append(tlon,drv5_lon,axis=0)
        tspd = np.append(tspd,drv5_spd,axis=0)
        if Vert_str5=="Pressure":
          tp = np.append(tp,drv5_prs,axis=0)
        else:
          tz = np.append(tz,drv5_hgt,axis=0)

if np.size(thlosspd) > 0:
  print("3D HLOS wind")
  sslabel = "HLOS Wind Velocity (m/s)"
  sslabelfile = "WSPD_HLOS"
  slat  = thloslat
  slon  = thloslon
  twind = thlosspd
  if Vert_str=="Pressure":
    tvert = thlosp
    del thlosp
  else:
    tvert = thlosz
    del thlosz
  del thloslat,thloslon,thlosspd
else:
  print("3D wind speed")
  sslabel = "Wind Speed (m/s)"
  sslabelfile = "WSPD"
  slat  = tlat
  slon  = tlon
  twind = tspd
  if Vert_str=="Pressure":
    tvert = tp
    del tp
  else:
    tvert = tz
    del tz
  del tlat,tlon,tspd

print("D vert_str = "+Vert_str)  
#zzlabel = Vert_str
#tmp_name = x_name
#map_3d_prof(twind,slon,slat,tvert,x_name,tmp_name,dateIN,sslabel,zzlabel)
#outname = output_path+"MAP_3D."+str(sslabelfile)+"."+str(dateIN)+"."+str(x_name)+avgthin_str+".png"
#plt.savefig(outname)
#del outname,Vert_str,slon,slat,sslabel,sslabelfile,twind,tvert

	# DEPENDENT DATASETS - WIND SPEED
Vert_str = []

for i in range(ndset):
  if idx_file_vertstr[i] != "NO_MATCHES":			# check if DEPENDENT dataset has any matches to DRIVER. If yes, proceed.
    print("dens: i = "+str(i))
    if i == 0:
      tname = dset1_name
      txspd = drv1_spd
      tyspd = dset1_spd
      txlat = drv1_lat
      txlon = drv1_lon
      Vert_str = Vert_str1
      if Vert_str1=="Pressure":
        txp = drv1_prs
        typ = dset1_prs
      else:
        if np.isnan(drv1_hgt).all():
          txz = prs_to_hgt(drv1_prs)
        else:
          txz = drv1_hgt
        tyz = dset1_hgt
    if i == 1:
      tname = dset2_name
      txspd = drv2_spd
      tyspd = dset2_spd
      txlat = drv2_lat
      txlon = drv2_lon
      Vert_str = Vert_str2
      if Vert_str2=="Pressure":
        txp = drv2_prs
        typ = dset2_prs
      else:
        if np.isnan(drv2_hgt).all():
          txz = prs_to_hgt(drv2_prs)
        else:
          txz = drv2_hgt
        tyz = dset2_hgt
    if i == 2:
      tname = dset3_name
      txspd = drv3_spd
      txlat = drv3_lat
      txlon = drv3_lon
      tyspd = dset3_spd
      Vert_str = Vert_str3
      if Vert_str3=="Pressure":
        txp = drv3_prs
        typ = dset3_prs
      else:
        if np.isnan(drv3_hgt).all():
          txz = prs_to_hgt(drv3_prs)
        else:
          txz = drv3_hgt
        tyz = dset3_hgt
    if i == 3:
      tname = dset4_name
      txspd = drv4_spd
      txlat = drv4_lat
      txlon = drv4_lon
      tyspd = dset4_spd
      Vert_str = Vert_str4
      if Vert_str4=="Pressure":
        txp = drv4_prs
        typ = dset4_prs
      else:
        if np.isnan(drv4_hgt).all():
          txz = prs_to_hgt(drv4_prs)
        else:
          txz = drv4_hgt
        tyz = dset4_hgt
    if i == 4:
      tname = dset5_name
      txspd = drv5_spd
      txlat = drv5_lat
      txlon = drv5_lon
      tyspd = dset5_spd
      Vert_str = Vert_str5
      if Vert_str5=="Pressure":
        txp = drv5_prs
        typ = dset5_prs
      else:
        if np.isnan(drv5_hgt).all():
          txz = prs_to_hgt(drv5_prs)
        else:
          txz = drv5_hgt
        tyz = dset5_hgt
    
    if tname == 'Aeolus':
      tname += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str
      sslabel = "HLOS Wind Velocity (m/s)"
      sslabelfile = "WSPD_HLOS"
    else:
      sslabel = "Wind Speed (m/s)"
      sslabelfile = "WSPD"

	# plot 3D maps
    if Vert_str=="Pressure":
      txvert = txp
    else:
      if np.isnan(txz).all():
        print("txz all NaN")
        txvert = prs_to_hgt(txp)
      else:
        txvert = txz
      print("txvert hgt = "+str(txvert))

      		#WIND SPEED
    sslabel  = "Wind Speed (m/s)"
    zzlabel  = Vert_str
    			# DEPENDENT DATASET
    tspd     = tyspd
    tmp_name = tname
    map_3d_prof(tspd,txlon,txlat,txvert,x_name,tmp_name,dateIN,sslabel,zzlabel)
    outname = output_path+"MAP_3D.y_"+str(sslabelfile)+".colloc_with_DRV_"+str(x_name)+"."+str(dateIN)+"."+str(tmp_name)+avgthin_str+".png"
    plt.savefig(outname)
    del outname

  	# plot density scatterplots
		#WIND SPEED
    units = "m/s"
    density_scatter(txlat,txlon,txspd, tyspd, x_name, tname, units)
	# plot filename
    outname = output_path+"DENSITY_SCATTER.Wind."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
	# save plot
    plt.savefig(outname)
    del outname
    		#PRESSURE
    if Vert_str=="Pressure":
      units = "hPa"
      density_scatter(txlat,txlon,txp, typ, x_name, tname, units)
	# plot filename
      outname = output_path+"DENSITY_SCATTER.Pressure."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
	# save plot
      plt.savefig(outname)
      del txp,typ,outname
    else:
    		#HEIGHT
      units = "km"
      density_scatter(txlat,txlon,txz, tyz, x_name, tname, units)
	# plot filename
      outname = output_path+"DENSITY_SCATTER.Height."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
	# save plot
      plt.savefig(outname)
      del txz,tyz,outname
	
    del txspd,tyspd,tname,txlat,txlon
 
    print("end density_scatter") 

     
      
#==================================================================
    
print("***** END MAIN PROGRAM *****")
#########################################################################################################
# END SCRIPT
#########################################################################################################
