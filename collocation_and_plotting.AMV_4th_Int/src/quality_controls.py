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
import numpy as np #....................................................... Array module
#import dask.array as da #.................................................. Dask module
#from dask.distributed import Client, LocalCluster #........................ Dask client modules (LocalCluster == runs on a local machine)
import datetime as dt #.................................................... Datetime module
#from geopy.distance import geodesic #...................................... Geodesic distance module
import time #.............................................................. Time module
from netCDF4 import Dataset #.............................................. netCDF module
#
###########################################################################
#
# QUALITY CONTROL functions for input datasets
#
# -------------------------------------------------------------------------
# QC Aeolus Winds
#	Based on QC recommended by ESA/ECMWF
#
#	INPUTS:
#		driver_type ......................... Aeolus wind type (i.e., RayClear, MieCloud)
#		prs ................................. pressure
#		err ................................. wind error
#		length .............................. accumulation length
#		htop ................................ height of top of range bin
#		hbot ................................ height of bottom of range bin
#
#	OUTPUTS:
#		idx ................................. indices where dataset PASSES QC
#	 	qc_list ............................. string listing all QC criteria used
#
def qc_aeolus(driver_type,prs,err,length,htop,hbot):
    if driver_type == 'RayClear':
        # QC limits
        plimit 		= 200.		#pressure level that separates different error limits. Units = hPa
        pmax 		= 800.		#max pressure allowed. Units = hPa
        emax_phigh 	= 8.5		#max error at lower levels. Units = m/s
        emax_plow 	= 12.		#max error at upper levels. Units = m/s
        lenmin 		= 60.		#min accumulation length. Units = km
        hmin 		= 0.3		#min height of range bin. Units = km

	# Indices where dataset PASSES QC
        #	Using np.where: To get new array size where input conditions are True, 
	#	replace "and" with * (multiplication), "or" with + (addition).
        idx = np.where((prs <= pmax)*(((prs <= plimit)*(err <= emax_plow)) + ((prs > plimit)*(err <= emax_phigh)))*(length >= lenmin)*(abs(htop - hbot) >= hmin))

        qc_list = "max pressure = "+str(pmax)+" hPa, max error at lower levels (p > "+str(plimit)+" hPa) = "+str(emax_phigh)+" m/s, max error at upper levels (p < "+str(plimit)+" hPa) = "+str(emax_plow)+" m/s, min accumulation length = "+str(lenmin)+" km, min height of range bin = "+str(hmin)+" km"

    if driver_type == 'MieCloud':
	# QC limits
        pmax            = 800.          #max pressure allowed. Units = hPa
        emax_p          = 5.            #max error at lower levels. Units = m/s

        # Indices where dataset PASSES QC
        #       Using np.where: To get new array size where input conditions are True,
        #       replace "and" with * (multiplication), "or" with + (addition).
        idx = np.where((prs <= pmax)*(err <= emax_p))
        
        qc_list = "max pressure allowed = "+str(pmax)+" hPa, max error allowed = "+str(emax_p)+" m/s"

    return idx,qc_list
    
# -------------------------------------------------------------------------
# QC Aircraft Winds
#	Based on QC used in NOAA/NCEP operations
#
#	INPUTS:
#		pres ................................ minimum pressure allowed
#
#	OUTPUTS:
#		idx ................................. indices where dataset PASSES QC
#		qc_list ............................. string listing all QC criteria used
#
def qc_aircraft(pres):
    print("QC AIRCRAFT: TBD")

    # QC limits
#    pmin 	= 126.		#min pressure allowed. Units = hPa

    # Indices where dataset PASSES QC
    #		pccf ... Units = %
    #		pct  ... Units = %
#    idx = np.where(pres >= pmin)
#
#    qc_list = "min pressure allowed = "+str(pmin)
#
#    return idx,qc_list

# -------------------------------------------------------------------------
# QC Atmospheric Motion Vectors (AMVs)
#	Based on QC recommended by AMV community
#
#	INPUTS:
#		pccf ................................ percent confidence (quality indicator, QI) in data point
#		pct ................................. min percent allowed
#
#	OUTPUTS:
#		idx ................................. indices where dataset PASSES QC
#		qc_list ............................. string listing all QC criteria used
#
def qc_amv(qi,pct):
    # Indices where dataset PASSES QC
    #		pccf ... Units = %
    #		pct  ... Units = %
    idx = np.where((qi >= pct)*(qi <= 100.))		#pccf cannot exceed 100%

    qc_list = str(pct)+" %"

    return idx,qc_list


###########################################################################
