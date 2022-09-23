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

import os
import sys
import csv
import numpy as np #....................................................... Array module
#import dask.array as da #.................................................. Dask module
#from dask.distributed import Client, LocalCluster #........................ Dask client modules (LocalCluster == runs on a local machine)
import datetime as dt #.................................................... Datetime module
#from geopy.distance import geodesic #...................................... Geodesic distance module
import time #.............................................................. Time module

import statistics as stats

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize 
import matplotlib.ticker as mticker

from scipy.interpolate import interpn
from scipy.stats import gaussian_kde
from scipy.stats import pearsonr

#os.system("module load miniconda/3.8-s4")
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D

##from gifly import gif_maker

import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker

from tools_analysis import get_unique_values
from tools_analysis import count_unique_values

#
###########################################################################

	# missing value
fill = -999

	# font sizes for plots
fonttitle  = 16
fontaxis   = 14
fontlegend = 12

###########################################################################
#
# PLOTTING functions for collocation program
#

# -------------------------------------------------------------------------
# Density Scatter Plot
#
#	INPUTS:
#		x ................................... driver dataset
#		y ................................... dependent dataset
#		x_name .............................. driver dataset name
#		y_name .............................. dependent dataset name
#
#	OUTPUTS:
#		ax .................................. plot
#
def density_scatter(Dlat,Dlon,tx,ty,x_name,y_name,units,ax=None,sort=True,bins=20,**kwargs):
   
    	# create plot
    fig = plt.figure(figsize=(11,8.5))		# (x,y) = (horiz,vert)
    
    if ax is None:
      #fig, ax = plt.subplots()
      ax = plt.subplot()
      
    x = tx
    y = ty
    
    xy = np.vstack([x,y])
    color_array = gaussian_kde(xy)(xy)
    
    idx = color_array.argsort()
    x, y, color_array = x[idx], y[idx], color_array[idx]
    
    norm = Normalize(vmin = np.min(color_array), vmax = np.max(color_array))
    cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    cbar.ax.set_ylabel("Normalized Density of Points")

      # plot scatter and fill
    ax.scatter(x, y, c=color_array, s=50)

    if units=="m/s":
      if x_name.find('Aeolus') != -1 or y_name.find('Aeolus') != -1: 
        label = "HLOS Wind Velocity"
        axismin = -100.0
        axismax = 100.0
        txpos	= -95
        typos	= 95
        diffpos = 10
      else:
        label = "Wind Speed"
        axismin = 0.0
        axismax = 100.0
        txpos	= 5
        typos	= 95
        diffpos = 5
    elif units=="hPa":
      label = "Pressure"
      axismin = 0.0
      axismax = 1000.0
      txpos   = 950
      typos   = 50
      diffpos = -50
    elif units=="km":
      label = "Height"
      axismin = 0.0
      axismax = 20.0
      txpos   = 2
      typos   = 19
      diffpos = 1
    ax.set_xlim([axismin,axismax])
    ax.set_ylim([axismin,axismax])
    
    ax.axhline(y=0,color="black")	      # horizontal line
    ax.axvline(x=0,color="black")	      # vertical line
    ax.axline((0,0),slope=1.0,color="black")  # one-to-one line

      # add text inside the plot
    xpos = txpos				      # value on x-axis where text will begin
    ypos = typos				      # value on y-axis where text will begin
    ftsz = 8					      # font size of text_* (see below)
    
    	      # get UNIQUE counts for DRIVER dataset (for legend labels)
    		# make lat/lon strings
    Dlatlonlist = []
    for j in range(np.size(x)):
      Dlatlonlist.append(str(x[j])+","+str(y[j]))
    		    # count unique values of strings
    nDuniq_list = count_unique_values(Dlatlonlist) 
    		    # add obs counts to dataset name for legend
    Dlegendlabel = x_name+" count = "+str(nDuniq_list)
    ftsz = 10						  # font size of text
    ypos = ypos - diffpos
    plt.text(xpos, ypos, Dlegendlabel, fontsize = ftsz)

    legendlabel = y_name+" count = "+str(np.size(y))
    ypos = ypos - diffpos
    plt.text(xpos, ypos, legendlabel, fontsize = ftsz)
        
	# compute and print stats of differences
    diff = y - x
    tcorr = np.corrcoef(x,y)	 # correlation
    corr = tcorr[0,1]
    avg  = np.mean(diff)		  # mean
    sd   = np.std(diff) 		  # standard deviation
    text = "Difference Stats ("+str(y_name)+"-"+str(x_name)+"):"
    ypos = ypos - diffpos
    plt.text(xpos, ypos, text, fontsize = ftsz)
    textr = "r = "+str(np.round_(corr,decimals=2))
    ypos = ypos - diffpos/2
    plt.text(xpos, ypos, textr, fontsize = ftsz)
    textm = "mean = "+str(np.round_(avg,decimals=2))
    ypos = ypos - diffpos/2
    plt.text(xpos, ypos, textm, fontsize = ftsz)
    texts = "stddev = "+str(np.round_(sd,decimals=2))
    ypos = ypos - diffpos/2
    plt.text(xpos, ypos, texts, fontsize = ftsz)
   
    # write output to csv text file 
    if 'outname' and 'dfile_out' in kwargs:
        rmse = np.sqrt( avg**2 + sd **2 )
        outname = kwargs.get('outname')
        dfile_out = kwargs.get('dfile_out')
        row1 = f'{outname} {x_name} {nDuniq_list} {fill} {fill} {fill}'
        row2 = f'{outname} {y_name} {np.size(y)} {avg:.2f} {sd:.2f} {rmse:.2f}'
        rows = [row1.split(), row2.split()]
        write_csv(rows=rows, dfile_out=dfile_out)
    
    del diff,tcorr,corr,avg,sd,text,textr,textm,texts
 
      # plot title and axis labels    
    ax.set_title(label+" "+str(y_name)+" vs "+str(x_name), fontsize=fonttitle)
    ax.set_xlabel(x_name+" ("+str(units)+")", fontsize=fontaxis)
    ax.set_ylabel(y_name+" ("+str(units)+")", fontsize=fontaxis)
   
    if units=="hPa":
      plt.gca().invert_xaxis()
      plt.gca().invert_yaxis()
 
    return ax
    
# -------------------------------------------------------------------------
# Histogram of Collocation Differences
#
#	INPUTS:
#		tmatch ............................... array of collocation dimensions (e.g., time, height, distance)
#		tname ................................ names of DEPENDENT datasets
#		alphaval ............................. transparency factor
#		ndset ................................ number of DEPENDENT datasets
#		match_str ............................ string of variable to plot
#
#	OUTPUTS:
#		fig .................................. plot
#
def hist_diffs(Dlat,Dlon,tmatch,x_name,tname,alphaval,match_str,dcolors,units,**kwargs):

	# create plot
    fig = plt.figure(figsize=(8.5,8.5))
    ax = plt.subplot()
	
	# histogram specs
    if match_str=="Time":
      label    = "Time Differences"
      hist_dir = "vertical"
      histmin  = -100
      histmax  = 100
      binsize  = 10
      bins   = int((2*int(histmax))/binsize)
      xpos = histmin + binsize
      ypos = histmax*10
      posdiff = binsize
    elif match_str=="Pressure":
      label    = "Pressure Differences"
      hist_dir = "horizontal"
      histmin  = -50
      histmax  = 50			# units = hPa
      binsize  = 5
      bins   = int((2*int(histmax))/binsize)
      xpos = 200
      ypos = histmax - 20
      posdiff = -1*binsize
    elif match_str=="Height":
      label    = "Height Differences"
      hist_dir = "horizontal"
      histmin  = -2000
      histmax  = 2000			# units = m
      binsize  = 100
      bins   = int((2*int(histmax))/binsize)
      xpos = 200
      ypos = histmax - 200
      posdiff = binsize
    elif match_str=="Distance":
      label    = "Collocation Distances"
      hist_dir = "vertical"
      histmin  = 0
      histmax  = 100			# units = km
      binsize  = 10
      bins     = int((2*int(histmax))/binsize)
      xpos = 75
      ypos = histmax*4
      posdiff = binsize

    if hist_dir == "vertical":
      xlabel = label+" ("+units+")"
      ylabel = "Number of Matched Pairs"
    elif hist_dir == "horizontal":
      ylabel = label+" ("+units+")"
      xlabel = "Number of Matched Pairs"

    trange = (histmin,histmax) 
    w      = 0.6			# width of bars

    xpos = xpos				# value on x-axis where text will begin
    ypos = ypos - posdiff 		# value on y-axis where text will begin

	# count unique DRIVER points
    Dlatlonlist = []
    xstat = []
    ystat = []
    for i in range(np.size(dcolors)):
      txx = Dlon[i]
      tyy = Dlat[i]
      print("size txx = "+str(np.size(txx)))
	    
      x   = txx
      y   = tyy
      xstat.append(x)
      ystat.append(y)
      del txx,tyy
              # count unique values of strings
      for j in range(np.size(x)):
        Dlatlonlist.append(str(x[j])+","+str(y[j]))

	# get unique number of DRIVER obs
    nDuniq_list = count_unique_values(Dlatlonlist)
    		# add obs counts to dataset name for legend
    Dlegendlabel = x_name+" count = "+str(nDuniq_list)
    ftsz = 10                                                   # font size of text
    t = plt.text(xpos, ypos, Dlegendlabel, fontsize = ftsz)
    t.set_bbox(dict(facecolor='white', alpha=0.75)) #, edgecolor='red'))

 	# get counts per dependent dataset for legend labels
    legendlabel = np.nan * np.ones(np.size(tname), dtype=list)
    rows = [] # rows to append to csv stats file
    for j in range(np.size(tname)):
      avg  = np.mean(tmatch[j])             # mean
      sd   = np.std(tmatch[j])              # standard deviation
      text = "mean = "+str(np.round_(avg,decimals=2))+" ... stddev = "+str(np.round_(sd,decimals=2))
      legendlabel[j] = tname[j]+" | count = "+str(np.size(tmatch[j]))+" ... "+text
      
      if 'outname' and 'dfile_out' in kwargs:
        outname = kwargs.get('outname')
        rmse = np.sqrt( avg**2 + sd **2 )
        row = f'{outname} {tname[j]} {np.size(tmatch[j])} {avg:.2f} {sd:.2f} {rmse:.2f}'
        rows.append(row.split())
	
    
    # write output to csv text file 
    if 'outname' and 'dfile_out' in kwargs:
        outname = kwargs.get('outname')
        dfile_out = kwargs.get('dfile_out')
        row = f'{outname} {x_name} {nDuniq_list} {fill} {fill} {fill}'
        rows.insert(0, row.split())
        write_csv(rows=rows, dfile_out=dfile_out)
	
    # plot
    ax.hist(tmatch, bins, trange, color=dcolors[0:np.size(tname)], histtype='bar', rwidth=w, alpha=alphaval, label=legendlabel, orientation=hist_dir)

    	# plot labels
    plt.title("Histogram of "+label+" relative to "+x_name+" Obs", fontsize=fonttitle)
    plt.xlabel(xlabel, fontsize=fontaxis)
    plt.ylabel(ylabel, fontsize=fontaxis)

    	# add legend to plot
	#	'bbox_to_anchor' = location of legend box. (x,y) = (left=0 and right=1, bottom=0 and top=1)
	# 	'center'         = center of bounding box at coords 'bbox_to_anchor'
	# 	'centerleft '    = center of left edge of bounding box at coords 'bbox_to_anchor'
	# 	'center right'   = center of right edge of bounding box at coords 'bbox_to_anchor'
    #leg = plt.legend(loc='center',bbox_to_anchor=(0.8,0))#,scatterpoints=1)
    leg = plt.legend(loc='upper right', prop={'size': fontlegend})

    for lh in leg.legendHandles: 
      lh.set_alpha(1)

    legendmarkersize = 30
    for lh in range(np.size(leg.legendHandles)):
      leg.legendHandles[lh]._sizes = [legendmarkersize]

    return fig

# -------------------------------------------------------------------------
# Map of Locations of Matched (Collocated) Observations
#
#	Cylindrical Equidistant Projection
#
#	INPUTS:
#		Ds ................................... driver wind speeds
#		s .................................... dependent wind speeds
#		Dx ................................... driver longitudes
#		x .................................... dependent longitudes
#		Dy ................................... driver latitudes
#		y .................................... dependent latitdues
#		x_name ............................... driver dataset name
#		y_name ............................... dependent dataset name
#		region_flag .......................... region to plot (0=global)
#		marksize ............................. marker size
#		alphaval ............................. transparency factor
#		dcolors .............................. array of marker colors
#
#	OUTPUTS:
#		fig .................................. plot
#
def map_locations_ce(Ds,ss,Dx,tx,Dy,ty,x_name,y_name,region_flag,marksize,imark,alphaval,dcolors,datein,**kwargs):

	# create plot
    fig = plt.figure(figsize=(11,8.5))

	# create map
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1,1,1,projection=proj)

    ax.coastlines()
    
    gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlines = True		# plot longitudes
    gl.ylines = True		# plot latitudes
    gl.xlabels_top=None
    gl.ylabels_right=None

    if region_flag==0:
	# set GLOBAL domain limits
      latmin = -90.0
      latmax = 90.0
      lonmin = -180.0
      lonmax = 180.0
      clat   = (latmax + latmin) / 2
      clon   = (lonmax + lonmin) / 2
      ax.set_extent([lonmin,lonmax,latmin,latmax], crs=proj)
      	# define x-axis tickmarks
      #ax.set_xticks(np.arange(lonmin,lonmax+1,60), crs=proj)
      gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
	# define x-axis tickmarks
      #ax.set_yticks(np.arange(latmin,latmax+1,30), crs=proj)
      gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])
		
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.ylabel_style = {'size': 15, 'color': 'gray'}

    	# plot DRIVER points
    Dlatlonlist = []
    Dlons = []
    Dlats = []
    for i in range(np.size(y_name)):
      tDs = Ds[i]
      tss = ss[i]
      
      tmp_x = Dx[i]
      tmp_y = Dy[i]
              # count unique values of strings
      for j in range(np.size(tmp_x)):
        Dlatlonlist.append(str(tmp_x[j])+","+str(tmp_y[j]))
      Dlons = np.append(Dlons,tmp_x,axis=0)
      Dlats = np.append(Dlats,tmp_y,axis=0)
      del tmp_x,tmp_y
      del tDs,tss  
	  
	# count unique values of strings
    nDuniq_list = count_unique_values(Dlatlonlist) 
    	# add obs counts to dataset name for legend
    Dlegendlabel = x_name+", count = "+str(nDuniq_list)
    	# plot DRIVER on map
    ax.scatter(Dlons, Dlats, c="black", s=marksize+5, marker="o", alpha=0.75, label=Dlegendlabel)
#      del x,y,tDx,tDy
#    del idx,Dxall,Dlatlonlist,nDuniq_list,Dlegendlabel

    rows = [] # rows of stats to be outputted in stats csv file
	# plot DEPENDENT points
    for i in range(np.size(y_name)):
      tname = y_name[i]

      tDs = Ds[i]
      tss = ss[i]
      
      x   = tx[i]
      y   = ty[i]
              # add obs counts to dataset name for legend
      legendlabel  = tname+", count = "+str(np.size(x)) #str(nuniq_list)

      if 'outname' and 'dfile_out' in kwargs:
        outname = kwargs.get('outname')
        row = f'{outname} {tname} {np.size(x)} {fill} {fill} {fill}'
        rows.append(row.split())

              # plot DEPENDENTS on map
      ax.scatter(x, y, c=dcolors[i], s=marksize, marker=imark[i], alpha=alphaval, label=legendlabel)
      del x,y,tDs,tss,legendlabel
    
    # write output to csv text file 
    if 'outname' and 'dfile_out' in kwargs:
      outname = kwargs.get('outname')
      dfile_out = kwargs.get('dfile_out')
      row = f'{outname} {x_name} {nDuniq_list} {fill} {fill} {fill}'
      rows.insert(0, row.split())
      write_csv(rows=rows, dfile_out=dfile_out)
 
    plt.title("Locations of Collocated Obs for "+str(datein), fontsize=fonttitle)

	# add legend to plot
	#	'bbox_to_anchor' = location of legend box. (x,y) = (left=0 and right=1, bottom=0 and top=1)
	# 	'center'         = center of bounding box at coords 'bbox_to_anchor'
	# 	'centerleft '    = center of left edge of bounding box at coords 'bbox_to_anchor'
	# 	'center right'   = center of right edge of bounding box at coords 'bbox_to_anchor'
    #leg = plt.legend(loc='center',bbox_to_anchor=(0.8,0))#,scatterpoints=1)
    leg = plt.legend(loc='lower right', prop={'size': fontlegend})

    for lh in leg.legendHandles: 
      lh.set_alpha(1)

    legendmarkersize = 30
    for lh in range(np.size(leg.legendHandles)):
      leg.legendHandles[lh]._sizes = [legendmarkersize]
  
    return fig

# -------------------------------------------------------------------------
# Map of Locations of Matched (Collocated) Observations
#
#	Orthographic Projection
#
#	INPUTS:
#		Ds ................................... driver wind speeds
#		s .................................... dependent wind speeds
#		Dx ................................... driver longitudes
#		x .................................... dependent longitudes
#		Dy ................................... driver latitudes
#		y .................................... dependent latitdues
#		x_name ............................... driver dataset name
#		y_name ............................... dependent dataset name
#		region_flag .......................... region to plot (0=global)
#		marksize ............................. marker size
#		alphaval ............................. transparency factor
#		dcolors .............................. array of marker colors
#
#	OUTPUTS:
#		fig .................................. plot
#
def map_locations_ortho(Ds,ss,Dx,tx,Dy,ty,x_name,y_name,marksize,imark,alphaval,dcolors,datein,center_lon,center_lat,**kwargs):

	# create plot
    fig = plt.figure(figsize=(8.5,8.5))

	# create map
    proj = ccrs.Orthographic(central_longitude=center_lon,central_latitude=center_lat,globe=None)
    geo  = ccrs.Geodetic()
    ax = plt.axes(projection=proj)
    
    	# add coastlines
    ax.coastlines(resolution='50m')
    ax.set_global()
    ax.gridlines()

	# plot DRIVER points
    Dlatlonlist = []
    Dlons = []
    Dlats = []
    for i in range(np.size(y_name)):
      tDs = Ds[i]
      tss = ss[i]
      
      tmp_x = Dx[i]
      tmp_y = Dy[i]

		# transform lat/lon points to orthographic points
      x = np.asarray(tmp_x)
      y = np.asarray(tmp_y)
      points  = proj.transform_points(geo, x, y)
      del x,y
      x = points[:,0]
      y = points[:,1]
      del points
                
		# count unique values of strings
      for j in range(np.size(tmp_x)):
        Dlatlonlist.append(str(tmp_x[j])+","+str(tmp_y[j]))
      del tmp_x,tmp_y

      Dlons = np.append(Dlons,x,axis=0)
      Dlats = np.append(Dlats,y,axis=0)
      del x,y
      del tDs,tss

        # count unique values of strings
    nDuniq_list = count_unique_values(Dlatlonlist)
        # add obs counts to dataset name for legend
    Dlegendlabel = x_name+", count = "+str(nDuniq_list)
        # plot DRIVER on map
    ax.scatter(Dlons, Dlats, c="black", s=marksize+5, marker="o", alpha=0.75, label=Dlegendlabel)

    rows = [] # rows of stats to be outputted to stats csv file
    # plot DEPENDENT points
    for i in range(np.size(y_name)):
      tname = y_name[i]

      tDs = Ds[i]
      tss = ss[i]
      
      idx_x = tx[i]
      idx_y = ty[i]

              # transform lat/lon points to orthographic points
      x = np.asarray(idx_x)
      y = np.asarray(idx_y)
      points  = proj.transform_points(geo, x, y)
      del idx_x,idx_y
			# add obs counts to dataset name for legend
      legendlabel  = tname+", count = "+str(np.size(x)) #str(nuniq_list)
      
      if 'outname' and 'dfile_out' in kwargs:
        outname = kwargs.get('outname')
        row = f'{outname} {tname} {np.size(x)} {fill} {fill} {fill}'
        rows.append(row.split())

      del x,y

      ttx = points[:,0]
      tty = points[:,1]
      x   = ttx
      y   = tty
              # plot DEPENDENTS on map
      ax.scatter(x, y, c=dcolors[i], s=marksize, marker=imark[i], alpha=alphaval, label=legendlabel)
      del x,y,tDs,tss,ttx,tty,legendlabel #,nuniq_list,latlonlist
  
    # write output to csv text file 
    if 'outname' and 'dfile_out' in kwargs:
      outname = kwargs.get('outname')
      dfile_out = kwargs.get('dfile_out')
      row = f'{outname} {x_name} {nDuniq_list} {fill} {fill} {fill}'
      rows.insert(0, row.split())
      write_csv(rows=rows, dfile_out=dfile_out)
 
    plt.title("Locations of Collocated Obs for "+str(datein), fontsize=fonttitle)

	# add legend to plot
	#	'bbox_to_anchor' = location of legend box. (x,y) = (left=0 and right=1, bottom=0 and top=1)
	# 	'center'         = center of bounding box at coords 'bbox_to_anchor'
	# 	'centerleft '    = center of left edge of bounding box at coords 'bbox_to_anchor'
	# 	'center right'   = center of right edge of bounding box at coords 'bbox_to_anchor'
    leg = plt.legend(loc='center',bbox_to_anchor=(0.8,0))#,scatterpoints=1)

    for lh in leg.legendHandles: 
      lh.set_alpha(1)

    legendmarkersize = 30
    for lh in range(np.size(leg.legendHandles)):
      leg.legendHandles[lh]._sizes = [legendmarkersize]

    return fig
    
# -------------------------------------------------------------------------
# Rotating Map of Locations of Matched (Collocated) Observations
#
#	Orthographic Projection
#
#	INPUTS:
#		Ds ................................... driver wind speeds
#		s .................................... dependent wind speeds
#		Dx ................................... driver longitudes
#		x .................................... dependent longitudes
#		Dy ................................... driver latitudes
#		y .................................... dependent latitdues
#		x_name ............................... driver dataset name
#		y_name ............................... dependent dataset name
#		region_flag .......................... region to plot (0=global)
#		marksize ............................. marker size
#		alphaval ............................. transparency factor
#		dcolors .............................. array of marker colors
#
#	OUTPUTS:
#		fig .................................. plot
#
def map_locations_ortho_rotate(outpath,outname,Ds,ss,Dx,tx,Dy,ty,x_name,y_name,marksize,imark,alphaval,dcolors,datein,center_lat,**kwargs):

	# central points
    center_lons = [*range(0,360,1)]		# should unpack all values from 0 to (360-1)=359 with an increment of 1
    size_lons = np.size(center_lons)
    #print("size lons = "+str(np.size(center_lons)))
    #print(center_lons)

	# create directory to store all images for gif creation
    gifdir 	= "gif_images"+str(datein)+"/"
    outpath_gif = outpath+gifdir

    os.system('if [ -d '+outpath_gif+' ]; then rm -Rf '+outpath_gif+'; fi')		# remove old gif directory before continuing
    os.system("mkdir -m 775 -p "+outpath_gif)						# make new empty gif directory
    
    tgifname 	= "image_"

	# create ortho plot for each central point
    print("create ROTATE images")
    for i in range(size_lons):
      #print("center point: lat = "+str(center_lat)+" | lon = "+str(center_lons[i]))
		# make each plot
      map_locations_ortho(Ds,ss,Dx,tx,Dy,ty,x_name,y_name,marksize,imark,alphaval,dcolors,datein,center_lons[i],center_lat)

		# assign number to each plot
      if (i+1) <= 100:
        gifnum = "0"+str(i)
        if (i+1) <= 10:
          gifnum = "00"+str(i)
      else:
        gifnum = str(i)

    		# save plot
      plt.savefig(outpath_gif+tgifname+gifnum+".png")

    del size_lons

	# make gif
    delay_sec = 15 #25 #50 #100				# time to display each image
    infiles   = outpath_gif+"image_*.png"	# all images to use to create gif
    outfile   = outname+".gif"			# filename for output gif
    print("infiles = "+infiles)
    print("outfile = "+outfile)    

    		# linux command to create gif
    cmd = "convert -delay "+str(delay_sec)+" -loop 0 "+str(infiles)+" "+str(outfile)
    print("cmd = "+cmd)
    		# run 'cmd'
    print("create ROTATE gif")
    os.system(cmd)
    
    #os.system("rm -r "+str(outpath_gif))

# -------------------------------------------------------------------------
# 3D Map of Matched (Collocated) Observations
#	
#	Cylindrical Equidistant Projection
#
#	INPUTS:
#		Ds ................................... driver wind speeds
#		s .................................... dependent wind speeds
#		Dx ................................... driver longitudes
#		x .................................... dependent longitudes
#		Dy ................................... driver latitudes
#		y .................................... dependent latitdues
#		x_name ............................... driver dataset name
#		y_name ............................... dependent dataset name
#		region_flag .......................... region to plot (0=global)
#		marksize ............................. marker size
#		alphaval ............................. transparency factor
#		dcolors .............................. array of marker colors
#
#	OUTPUTS:
#		fig .................................. plot
#
def map_3d_prof(ss,xx,yy,zz,x_name,tname,datein,sslabel,zzlabel,**kwargs):

	# longitude check: set range to 0=360
    for iP in range(np.size(xx)):
      if xx[iP] > 180.0:
        xx[iP] = xx[iP] - 360.

	# create plot
    fig = plt.figure(figsize=(11.5,8.5))

	# create map
    proj = "3d"
    ax = fig.gca(projection=proj)
    
    	# Define lower left, uperright lontitude and lattitude respectively
    extent = [-180, 180, -90, 90] 
    	# Create a basemap instance that draws the Earth layer
    bm = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[2],
             urcrnrlon=extent[1], urcrnrlat=extent[3],
             projection='cyl', resolution='l', fix_aspect=False, ax=ax)

    	# set z-axis range
    if zzlabel=="Pressure":
      zs_val = 1000
      ax.set_zlim(0., 1000.)
      plt.gca().invert_zaxis()
    elif zzlabel=="Height":
      zs_val = 0
      ax.set_zlim(0.,20.)
    
	# Add Basemap to the figure
    ax.add_collection3d(bm.drawcoastlines(linewidth=0.25),zs=zs_val)
    ax.add_collection3d(bm.drawcountries(linewidth=0.35),zs=zs_val)

	# Add meridian and parallel gridlines
    ax.set_xlabel('Longitude', labelpad=20)
    ax.set_ylabel('Latitude', labelpad=20)
    if zzlabel=="Pressure":
      zzunits = "(hPa)"
    elif zzlabel=="Height":
      zzunits = "(km)"
    ax.set_zlabel(zzlabel+" "+zzunits, labelpad=20)

    lon_step = 30
    lat_step = 30
    meridians = np.arange(extent[0], extent[1] + lon_step, lon_step)
    parallels = np.arange(extent[2], extent[3] + lat_step, lat_step)
    ax.set_yticks(parallels)
    ax.set_yticklabels(parallels)
    ax.set_xticks(meridians)
    ax.set_xticklabels(meridians)

    	# 3D viewpoint
    azim_pt = 240 #230=lower left corner is center; less is to the left, more is to the right
    elev_pt = 25  #50
    ax.view_init(azim=azim_pt, elev=elev_pt)

    print("size xx = "+str(np.size(xx))+" yy = "+str(np.size(yy))+" zz = "+str(np.size(zz)))

	# add counts
		# DRIVER: get UNIQUE counts for DRIVER dataset (for legend labels)
    Dlatlonlist = []
    for j in range(np.size(xx)):
      
      Dlatlonlist.append(str(xx[j])+","+str(yy[j])+","+str(zz[j]))
                      # count unique values of strings
    nDuniq_list = count_unique_values(Dlatlonlist)
                      # add obs counts to dataset name for legend
    Dlegendlabel = x_name+" count = "+str(nDuniq_list)
    ftsz = 10                                             # font size of text
    ax.text2D(0.05, 0.95, Dlegendlabel, transform=ax.transAxes)

    if tname!=x_name:  
		# DEPENDENT
      legendlabel = tname+" count = "+str(np.size(yy))
      ax.text2D(0.05, 0.90, legendlabel, transform=ax.transAxes)
    
    	# add scatterplot map based on lons 'xx', lats 'yy', heights 'zz', and colors based on 'ss'
    p = ax.scatter(xx, yy, zz, c=ss, cmap="jet")
    	# add colorbar to reference intensity of 'ss'
    fig.colorbar(p, label=sslabel)
    
    plt.title(tname+" "+sslabel+" for "+str(datein), fontsize=fonttitle)
    
        # write output to csv text file 
    if 'outname' and 'dfile_out' in kwargs:
        outname = kwargs.get('outname')
        dfile_out = kwargs.get('dfile_out')
        row1 = f'{outname} {x_name} {nDuniq_list} {fill} {fill} {fill}'
        row2 = f'{outname} {tname} {np.size(yy)} {fill} {fill} {fill}'
        rows = [row1.split(), row2.split()]
        write_csv(rows=rows, dfile_out=dfile_out)

    return fig
 
# -------------------------------------------------------------------------
# Scatter Plot (not density) of Matched Observations
#
#	INPUTS:
#		x ................................... driver dataset
#		y ................................... dependent dataset
#		x_name .............................. driver dataset name
#		y_name .............................. dependent dataset name
#
#	OUTPUTS:
#		ax .................................. plot
#
def scatter_matches(Dlat,Dlon,tx,ty,x_name,tname,marksize,alphaval,match_str,acolors,units,imark,**kwargs):

	# create plot
    fig = plt.figure(figsize=(9,8.5))
    ax = plt.subplot()

	# scatter specs
    if match_str=="HLOS":
      label = "HLOS Wind Velocity"
      axismin = -100.0
      axismax = 100.0
      diffpos = 5
    elif match_str=="Wind Speed":
      label = match_str
      axismin = 0.0
      axismax = 100.0
      diffpos = 5
    elif match_str=="Pressure":
      label = match_str
      axismin = 0.0
      axismax = 1000.0
      diffpos = -50
    elif match_str=="Height":
      label = match_str
      axismin = 0.0
      axismax = 20.0
      diffpos = 1

	# axis limits
    ax.set_xlim([axismin,axismax])
    ax.set_ylim([axismin,axismax])

    if match_str=="Pressure":
      xpos = axismax-50
      ypos = axismin+50
    elif match_str=="Height":
      xpos = axismin+1
      ypos = axismax-1
    else:
      xpos = axismin+10
      ypos = axismax-10

    	# get UNIQUE counts for DRIVER dataset (for legend labels)
    Dlatlonlist = []
    xstat = []
    ystat = []
    rows = [] # rows of data to append to csv stats file
    for i in range(np.size(acolors)):
      txx = tx[i]
      tyy = ty[i]
      print("size txx = "+str(np.size(txx)))

      x   = txx
      y   = tyy
      del txx,tyy

              # count unique values of strings
      for j in range(np.size(x)):
        Dlatlonlist.append(str(x[j])+","+str(y[j]))

          # adding text inside the plot
      if match_str=="Pressure":
        xpos = xpos
        ypos = ypos + 50
      elif match_str=="Height":
        xpos = xpos
        ypos = ypos - 1
      else:
        xpos = xpos						# value on x-axis where text will begin
        ypos = ypos						# value on y-axis where text will begin

              # compute and print stats of differences
      ftsz = 8
      diff = y - x
      tcorr = np.corrcoef(x,y)	 # correlation
      corr = tcorr[0,1]
      avg  = np.mean(diff)		  # mean
      sd   = np.std(diff) 		  # standard deviation
      text = "Difference Stats ("+str(tname[i])+"-"+str(x_name)+"):"
      ypos = ypos - diffpos
      plt.text(xpos, ypos, text, fontsize = ftsz)
      textr = "r = "+str(np.round_(corr,decimals=2))
      ypos = ypos - diffpos/2
      plt.text(xpos, ypos, textr, fontsize = ftsz)
      textm = "mean = "+str(np.round_(avg,decimals=2))
      ypos = ypos - diffpos/2
      plt.text(xpos, ypos, textm, fontsize = ftsz)
      texts = "stddev = "+str(np.round_(sd,decimals=2))
      ypos = ypos - diffpos/2
      plt.text(xpos, ypos, texts, fontsize = ftsz)
      
      if 'outname' and 'dfile_out' in kwargs:
        outname = kwargs.get('outname')
        rmse = np.sqrt( avg**2 + sd **2 )
        row = f'{outname} {tname[i]} {np.size(x)} {avg:.2f} {sd:.2f} {rmse:.2f}'
        rows.append(row.split())

      del diff,tcorr,corr,avg,sd,text,textr,textm,texts

	    # legend label
      legendlabel = tname[i]+" | count = "+str(np.size(x))

          # plot scatter
      plt.scatter(x, y, c=acolors[i], marker=imark[i], edgecolor="black", s=marksize, alpha=alphaval, label=legendlabel)

      del x,y

	# get unique number of DRIVER obs
    nDuniq_list = count_unique_values(Dlatlonlist)
    		# add obs counts to dataset name for legend
    Dlegendlabel = x_name+" count = "+str(nDuniq_list)
    ftsz = 10                                                   # font size of text
    plt.text(xpos, ypos-diffpos, Dlegendlabel, fontsize = ftsz)

    # write output to csv text file 
    if 'outname' and 'dfile_out' in kwargs:
      outname = kwargs.get('outname')
      dfile_out = kwargs.get('dfile_out')
      row = f'{outname} {x_name} {nDuniq_list} {fill} {fill} {fill}'
      rows.insert(0, row.split())
      write_csv(rows=rows, dfile_out=dfile_out)

	# plot labels
    plt.axhline(y=0,color="black")  	      # horizontal line
    plt.axvline(x=0,color="black")  	      # vertical line
    plt.axline((0,0),slope=1.0,color="black")     # one-to-one line

    plt.title(label+" Scatterplot: DEPENDENTS vs DRIVER", fontsize=fonttitle)
    plt.xlabel(x_name+" ("+units+")", fontsize=fontaxis)
    plt.ylabel("DEPENDENT dataset(s)", fontsize=fontaxis)

    if match_str=="Pressure":
      plt.gca().invert_xaxis()
      plt.gca().invert_yaxis()

      	# add legend to plot
	#	'bbox_to_anchor' = location of legend box. (x,y) = (left=0 and right=1, bottom=0 and top=1)
	# 	'center'         = center of bounding box at coords 'bbox_to_anchor'
	# 	'centerleft '    = center of left edge of bounding box at coords 'bbox_to_anchor'
	# 	'center right'   = center of right edge of bounding box at coords 'bbox_to_anchor'
    #leg = plt.legend(loc='center',bbox_to_anchor=(0.8,0))#,scatterpoints=1)
    leg = plt.legend(loc='lower right', prop={'size': fontlegend})

    for lh in leg.legendHandles: 
      lh.set_alpha(1)

    legendmarkersize = 30
    for lh in range(np.size(leg.legendHandles)):
      leg.legendHandles[lh]._sizes = [legendmarkersize]

    return fig


def write_csv(rows, dfile_out):
    file_exists = os.path.isfile(dfile_out)
    with open(dfile_out, 'a') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(['plot_filename', 'dataset', 'N', 'mean_diff', 'std_diff', 'rmse'])
        writer.writerows(rows)
# -------------------------------------------------------------------------
