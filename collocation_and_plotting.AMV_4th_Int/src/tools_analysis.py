###########################################################################
#
# PYTHON 3 FUNCTIONS FOR quality_controls
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

import sys
import math
import numpy as np #....................................................... Array module
#import dask.array as da #.................................................. Dask module
#from dask.distributed import Client, LocalCluster #........................ Dask client modules (LocalCluster == runs on a local machine)
import datetime as dt #.................................................... Datetime module
#from geopy.distance import geodesic #...................................... Geodesic distance module
import time #.............................................................. Time module

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
#from mpl_toolkits.basemap import Basemap
#from scipy.stats import gaussian_kde

from statistics import mean

#
###########################################################################
#
# STATISTICAL ANALYSIS functions for collocation program
#

fill = -999.0

# -------------------------------------------------------------------------
# Compute Horizontal Line-of-Sight (HLOS) Wind for using Aeolus azimuth angle
#
#	INPUTS:
#		x_azm ............................... Aeolus azimuth angle
#		y_dir ............................... 
#		y_spd ............................... driver dataset name
#
#	OUTPUTS:
#		idx ................................. indices where dataset PASSES QC
#	 	qc_list ............................. string listing all QC criteria used
#
def compute_hlos(x_azm,y_dir,y_spd):

  sindir = math.sin(y_dir*(math.pi/180.0))
  cosdir = math.cos(y_dir*(math.pi/180.0))
  
  sinazm = -1.0*math.sin(x_azm*(math.pi/180.0))
  cosazm = -1.0*math.cos(x_azm*(math.pi/180.0))
  
  u = -1.0*y_spd*sindir
  v = -1.0*y_spd*cosdir
  
  hlos = u*sinazm + v*cosazm
  
  return hlos
    
# -------------------------------------------------------------------------
# Get Unique Values in List
#
def get_unique_values(list1):

    # initialize a null list
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    # print list
#    for x in unique_list:
#        print(x)

    return unique_list

# -------------------------------------------------------------------------
# Count Unique Values in List
#
def count_unique_values(list1):
    
    unique_list = get_unique_values(list1)

    count = np.size(unique_list)

    return count
    
# -------------------------------------------------------------------------
# Super-ob (average) matches
#
def superob_matches(Dlat,Dlon,Dprs,Dhgt,Dspd,Ddir,tlat,tlon,tprs,thgt,tspd,tdir):

    lat_uniq = get_unique_values(Dlat)
    lon_uniq = get_unique_values(Dlon)

    nlat_uniq = count_unique_values(Dlat)	# nlat_uniq should = nlon_uniq
    nlon_uniq = count_unique_values(Dlon)
    
    dDlatout   = []
    dDlonout   = []
    dDprsout   = []
    dDhgtout   = []
    dDspdout   = []
    dDdirout   = []
    tdset1_lat = []
    tdset1_lon = []
    tdset1_prs = []
    tdset1_hgt = []
    tdset1_spd = []
    tdset1_dir = []
    for jj in range(nlat_uniq):
      for kk in range(nlon_uniq):
        iloc = np.where((Dlat==lat_uniq[jj])*(Dlon==lon_uniq[kk]))
        print("num matches for lat "+str(lat_uniq[jj])+","+str(lon_uniq[kk])+" = "+str(np.size(iloc)))
    
        if np.size(iloc)!=0:  
          dDlatout.append(Dlat[jj])
          dDlonout.append(Dlon[jj])
          dDprsout.append(Dprs[jj])
          dDhgtout.append(Dhgt[jj])
          dDspdout.append(Dspd[jj])
          dDdirout.append(Ddir[jj])

          tdset1_lat.append(mean(tlat[iloc]))
          tdset1_lon.append(mean(tlon[iloc]))
          tdset1_prs.append(mean(tprs[iloc]))
          tdset1_hgt.append(mean(thgt[iloc]))
          tdset1_spd.append(mean(tspd[iloc]))
          tdset1_dir.append(mean(tdir[iloc]))
      
        del iloc

    Dlatout   = np.asarray(dDlatout)
    Dlonout   = np.asarray(dDlonout)
    Dprsout   = np.asarray(dDspdout)
    Dhgtout   = np.asarray(dDhgtout)
    Dspdout   = np.asarray(dDspdout)
    Ddirout   = np.asarray(dDdirout)

    dset1_lat = np.asarray(tdset1_lat)
    dset1_lon = np.asarray(tdset1_lon)
    dset1_prs = np.asarray(tdset1_prs)
    dset1_hgt = np.asarray(tdset1_hgt)
    dset1_spd = np.asarray(tdset1_spd)
    dset1_dir = np.asarray(tdset1_dir)

    return Dlatout,Dlonout,Dprsout,Dhgtout,Dspdout,Ddirout,dset1_lat,dset1_lon,dset1_prs,dset1_hgt,dset1_spd,dset1_dir

# -------------------------------------------------------------------------
# Get matched variables from matched indices
#
def matched_vars(avgthin_choice,driver_name,t_name,idx_drv_dset,idx_dset,tdrv_lat,tdrv_lon,tdrv_prs,tdrv_hgt,tdrv_spd,tdrv_dir,t_lat,t_lon,t_prs,t_hgt,t_spd,t_dir):

    drv_lat = tdrv_lat[idx_drv_dset]
    drv_lon = tdrv_lon[idx_drv_dset]
    drv_prs = tdrv_prs[idx_drv_dset]
    drv_hgt = tdrv_hgt[idx_drv_dset]
    drv_spd = tdrv_spd[idx_drv_dset]
    drv_dir = tdrv_dir[idx_drv_dset]
    
    dset_lat = t_lat[idx_dset]
    dset_lon = t_lon[idx_dset]
    dset_prs = t_prs[idx_dset]
    dset_hgt = t_hgt[idx_dset]
    dset_spd = t_spd[idx_dset]
    dset_dir = t_dir[idx_dset]
    
    	# omit missing values
    ispd = np.where((drv_spd != np.nan)*(drv_spd > fill)*(dset_spd != np.nan)*(dset_spd > fill))

    slat = drv_lat
    slon = drv_lon
    sprs = drv_prs
    shgt = drv_hgt
    sspd = drv_spd
    sdir = drv_dir
    del drv_lat,drv_lon,drv_prs,drv_hgt,drv_spd,drv_dir
    drv_lat = slat[ispd]
    drv_lon = slon[ispd]
    drv_prs = sprs[ispd]
    drv_hgt = shgt[ispd]
    drv_spd = sspd[ispd]
    drv_dir = sdir[ispd]
    del slat,slon,sprs,shgt,sspd,sdir
    slat = dset_lat
    slon = dset_lon
    sprs = dset_prs
    shgt = dset_hgt
    sspd = dset_spd
    sdir = dset_dir
    del dset_lat,dset_lon,dset_prs,dset_hgt,dset_spd,dset_dir
    dset_lat = slat[ispd]
    dset_lon = slon[ispd]
    dset_prs = sprs[ispd]
    dset_hgt = shgt[ispd]
    dset_spd = sspd[ispd]
    dset_dir = sdir[ispd]
    del slat,slon,sprs,shgt,sspd,sdir
    del ispd
    
    if t_name == 'Aeolus':
      print("dep Aeolus")
      print("size: "+str(np.size(dset_lat)))
      tmp_spd = drv_spd
      tmp_dir = drv_dir
      del drv_spd,drv_dir
      drv_spd  = np.nan * np.ones_like(drv_prs)
      drv_dir  = np.nan * np.ones_like(drv_prs)
      for jj in range(np.size(drv_prs)):
        drv_spd[jj] = compute_hlos(dset_dir[jj],tmp_dir[jj],tmp_spd[jj])
        drv_dir[jj] = dset_dir[jj]
      del tmp_spd,tmp_dir    

    if driver_name == 'Aeolus':
      print("drv Aeolus")
      tmp_spd = dset_spd
      tmp_dir = dset_dir
      del dset_spd,dset_dir
      dset_spd = np.nan * np.ones_like(dset_prs)
      dset_dir = np.nan * np.ones_like(dset_prs)
      for jj in range(np.size(dset_prs)):
        dset_spd[jj] = compute_hlos(drv_dir[jj],tmp_dir[jj],tmp_spd[jj])
        dset_dir[jj] = drv_dir[jj]
      del tmp_spd,tmp_dir
	  
	  # check if super-obbing or thinning is selected
    if avgthin_choice==-1:
      print("PLOT ALL MATCHES")
    elif avgthin_choice==0:
      print("SUPER-OB (AVERAGE) MATCHES")

      Dlat = drv_lat
      Dlon = drv_lon
      Dprs = drv_prs
      Dhgt = drv_hgt
      Dspd = drv_spd
      Ddir = drv_dir
      del drv_lat,drv_lon,drv_prs,drv_hgt,drv_spd,drv_dir
      tlat = dset_lat
      tlon = dset_lon
      tprs = dset_prs
      thgt = dset_hgt
      tspd = dset_spd
      tdir = dset_dir
      del dset_lat,dset_lon,dset_prs,dset_hgt,dset_spd,dset_dir 

      drv_lat,drv_lon,drv_prs,drv_hgt,drv_spd,drv_dir,dset_lat,dset_lon,dset_prs,dset_hgt,dset_spd,dset_dir = superob_matches(Dlat,Dlon,Dprs,Dhgt,Dspd,Ddir,tlat,tlon,tprs,thgt,tspd,tdir)
      del Dlat,Dlon,Dprs,Dhgt,Dspd,Ddir,tlat,tlon,tprs,thgt,tspd,tdir	      
    elif avgthin_choice==1:
      print("THIN MATCHES: TBA")
	  
    return drv_lat,drv_lon,drv_prs,drv_hgt,drv_spd,drv_dir,dset_lat,dset_lon,dset_prs,dset_hgt,dset_spd,dset_dir
	  
# -------------------------------------------------------------------------
# Convert pressure to height
#	Convert pressure to pressure altitude (height) following NWS formulation (https://www.weather.gov/media/epz/wxcalc/pressureAltitude.pdf)
#
# Inputs	prs = pressure in hPa
# Returns 	hgt = height in km
#
def prs_to_hgt(prs):	  

    	# check pressure units and convert to hPa
    if max(prs) > 10000.:
      prs = prs/100.

	# convert pressure to height
    hgt = np.nan * np.ones_like(prs)
    for i in range(np.size(prs)):
      hgt[i] = 145366.45 * (1.0 - (prs[i]/1013.25)**0.190284)			# convert hPa (mb) to feet
      hgt[i] = hgt[i] * 0.3048							# convert to meters
      hgt[i] = hgt[i] / 1000.0							# convert to km

      #print("prs_to_hgt: prs = "+str(prs[i])+" | hgt = "+str(hgt[i]))

    return hgt






#############################################################################
