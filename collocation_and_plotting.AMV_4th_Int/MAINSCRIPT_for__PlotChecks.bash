#!/bin/bash
##############################################################
# Main script to collocate wind datasets with Aeolus winds, and save to user-specified output directory
#
# 	In this script, the user chooses and assigns the following criteria for collocation: paths, dates, datasets, collocation limits, QC, HPC partition
#
# 	To run this script on the command line: ./NameOfThisScript.bash
#
#  echo 'dd = '$dd' | idy = '$idy' | ndd = '$ndd
# NOTES: 
#	1. This script runs all jobs in the background. If using S4 supercomputer, it is recommended to run everything on the
#	   's4' partition (6-hour runtime limit). The partition is assigned in this script.
#
# History:
#	2021-11-18	K.E.Lukens	NOAA/NESDIS/STAR ... Created
#	2022-06-08	K.E.Lukens	NOAA/NESDIS/STAR ... Expanded ingest capability and added comments throughout
#
##############################################################

#set -x

#==============================================
# BEGIN USER INPUT
#==============================================

#----------------------------------------------
# Set date range over which to run collocation
#	Collocation index files will be generated automatically for the dates specified below.
#	For computational efficiency, it is recommended to loop through one month at a time. 
#	Choose one year/month combination and an array of days.
#
# yyyy 	= year
# mm	= month
# ddarr = day array

	#`````````````````````````````````````
	# Year
yyyy=2019

	#`````````````````````````````````````
	# Month
mm=(10)

	#`````````````````````````````````````
	# Day array - choose ONE, comment out the rest

	# Months with 31 days: Jan, Mar, May, Jul, Aug, Oct, Dec
#ddarr=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31)
	# Months with 30 days: Apr, Jun, Sep, Nov
#ddarr=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
	# Feb - LEAP YEAR
#ddarr=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29)
	# Feb - NON-LEAP YEAR
#ddarr=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28)

	# Custom day array
ddarr=(20)

#----------------------------------------------
# Set full paths
#       !!! NOTE: Always end path names with a slash "/"

	# Set home path 'dir_home': This is the current working directory (i.e., where this script is located)
dir_home=/home/sreiner/Wind_Collocation_Software/collocation_and_plotting.AMV_4th_Int/
echo 'WORKING DIRECTORY = '${dir_home}

	# Set input path 'dir_in': This is where the collocation index files are located
dir_in=/data/users/klukens/for_CIMSS/Wind_Collocation_Software/index_files/
echo 'INPUT DIRECTORY = '${dir_in}

	# Set index file suffix
	#	This is the portion of the filename after the DATE. Include all punctuation.
	#	This should be the same for ALL index files.

		# Cloudy AMVs
#file_in='.drv_BRZ_Cloud__QC.dset1_Aeolus_MieCloud_ReprocB11__QC.dset2_Aircraft__NoQC.dset3_Radiosonde__NoQC.nc4'
#file_in='.drv_EUM_Cloud__QC.dset1_Aeolus_MieCloud_ReprocB11__QC.dset2_Aircraft__NoQC.dset3_Radiosonde__NoQC.nc4'
#file_in='.drv_JMA_Cloud__QC.dset1_Aeolus_MieCloud_ReprocB11__QC.dset2_Aircraft__NoQC.dset3_Radiosonde__NoQC.nc4'
#file_in='.drv_KMA_Cloud__QC.dset1_Aeolus_MieCloud_ReprocB11__QC.dset2_Aircraft__NoQC.dset3_Radiosonde__NoQC.nc4'
#file_in='.drv_NOA_Cloud__QC.dset1_Aeolus_MieCloud_ReprocB11__QC.dset2_Aircraft__NoQC.dset3_Radiosonde__NoQC.nc4'
#file_in='.drv_NWC_Cloud__QC.dset1_Aeolus_MieCloud_ReprocB11__QC.dset2_Aircraft__NoQC.dset3_Radiosonde__NoQC.nc4'
		# Clear AMVs
#file_in='.drv_BRZ_Clear__QC.dset1_Aeolus_RayClear_ReprocB11__QC.dset2_Aircraft__NoQC.dset3_Radiosonde__NoQC.nc4'
#file_in='.drv_EUM_Clear__QC.dset1_Aeolus_RayClear_ReprocB11__QC.dset2_Aircraft__NoQC.dset3_Radiosonde__NoQC.nc4'
#file_in='.drv_KMA_Clear__QC.dset1_Aeolus_RayClear_ReprocB11__QC.dset2_Aircraft__NoQC.dset3_Radiosonde__NoQC.nc4'
#file_in='.drv_NOA_Clear__QC.dset1_Aeolus_RayClear_ReprocB11__QC.dset2_Aircraft__NoQC.dset3_Radiosonde__NoQC.nc4'
#file_in='.drv_NWC_Clear__QC.dset1_Aeolus_RayClear_ReprocB11__QC.dset2_Aircraft__NoQC.dset3_Radiosonde__NoQC.nc4'

file_in=$1
	# Set output path 'dir_out': This is where the plots are saved
dir_out=/home/sreiner/Wind_Collocation_Software/output_plots/
echo 'OUTPUT DIRECTORY = '${dir_out}

        # Set path 'archive_parent': This is where the archive home directory is located
        #       Example: If the full path to /atmos-nc-dataset is /Full/Path/To/atmos-nc-dataset,
        #                then
        #                archive_prefix = /Full/Path/To/
        #
        #       Note: /atmos-nc-dataset and /aeolus... directories should be in the same 'archive_parent' directory
        #
archive_parent=/scratch/

#----------------------------------------------
# Select choice to super-ob (average) or thin multiple collocations per DRIVER observation, or plot all matches
#      -1 = plot all matches (i.e., do nothing)
#	0 = super-ob (average all matches)
#	1 = thin (for future use; capability TBA)

avgthin_choice=-1
#avgthin_choice=0

#----------------------------------------------
# Set HPC account, partition, and runtime limit for jobs

account="ssec"

partition="s4"
timelimit="06:00:00"
#partition="serial"
#timelimit="09:00:00"

#----------------------------------------------
# Load modules
#	This list is taylored for use on the S4 supercomputer.
#	If using another HPC, customize following its own protocol.

module purge

cd

module load license_intel/S4
module load intel/18.0.4  #emc-hpc-stack/2020-q3

module load hdf/
module load hdf5/
module load netcdf4/
module load udunits2
module load ncview
module load nco
module load ncl

#==============================================
# END USER INPUT
#==============================================

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!! USERS SHOULD NOT CHANGE ANYTHING BELOW THIS LINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#----------------------------------------------
# Find sizes of arrays

ndd=${#ddarr[@]}
#echo 'ndd = '$ndd

#----------------------------------------------
# Set working paths
#	Always end path names with a slash "/"

dir_plot=${dir_home}src/           #location of collocation source code
echo 'PLOTTING CODE DIRECTORY = '${dir_plot}

#----------------------------------------------
# Set name of main collocation code script (to pass as argument to run job script)

plotting_script=PLOT_MAIN.collocation_checks.py

#----------------------------------------------
# Set name of run job script (to pass as argument to run job script)

run_job_script=run_plotting_code.job


#==============================================
#==============================================
# LOOP thru days of selected month
#==============================================

#imm=0
#while [[ imm -lt $nmm ]]
#do
#  mm=${mmarr[$imm]}
#  echo 'mm = '$mm' | imm = '$imm' | nmm = '$nmm

idy=0
while [[ idy -lt $ndd ]]
do
  dd=${ddarr[$idy]}		#day

  echo "----- Date: $yyyy $mm $dd"
  yyyymmdd=$yyyy$mm$dd

  cd ${dir_home}

#----------------------------------------------
# Run Collocation (matching) algorithm using SLURM (sbatch command)
#	The following initiates collocations with Aeolus Rayleigh-clear and Mie-cloudy winds 
#	for all four 6-hour periods in a single day (+/- 3 hours around each center analysis 
#	time (00, 06, 12, and 18 UTC)). Each combination is run as a separate job, and is 
#	initiated at the same time; thus, 8 jobs will run at once (4 for each Aeolus wind type).
#
#	Run command format for each job:
#		sbatch $output_logfile_name $job_name $input_directory $colloc_script $arg1 /
#		$arg2 $arg3 $arg4 $arg5 $arg6 $arg7 $arg8 $partition $timelimit $run_job_script
#
#       Arguments (arg) to customize collocation:
#               1. Aeolus (driver) wind type: RayClear, MieCloud
#               2. Date in YYYYMMDDHH
#               3. Input directory for index files
#		4. Input index file suffix
#               5. Output directory
#		6. Archive parent path: path where home archive directory is located
#		7. Choice to super-ob, thin, or plot all matches

  #arg1: automatically set in loop below
  #arg2: automatically set in loop below
  arg3=${dir_in}
  arg4=${file_in}
  arg5=${dir_out}
  arg6=${archive_parent}
  arg7=${avgthin_choice}


if [[ "$file_in" == *"Aeolus"* ]]; then		# check dataset_names (all of them) if 'Aeolus' is present
	# 'Aeolus' is a listed dataset. Submit one job per 6-h analysis cycle (00, 06, 12, 18 UTC) per Aeolus wind type (Rayleigh-clear or Mie-cloudy) = 8 total jobs.

	#`````````````````````````````````````````
        # Collocate with Aeolus Rayleigh-clear

  arg1=RayClear

		#job logfile names
  Rout00=OUTplotRAY00_${dd}
  Rout06=OUTplotRAY06_${dd}
  Rout12=OUTplotRAY12_${dd}
  Rout18=OUTplotRAY18_${dd}

  rm $Rout00 $Rout06 $Rout12 $Rout18

		# hour = 00 UTC
  arg2=${yyyymmdd}00
  sbatch --output=${Rout00} --job-name=PlotRay00_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}
		# hour = 06 UTC
  arg2=${yyyymmdd}06
  sbatch --output=${Rout06} --job-name=PlotRay06_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}
		# hour = 12 UTC
  arg2=${yyyymmdd}12
  sbatch --output=${Rout12} --job-name=PlotRay12_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}
		# hour = 18 UTC
  arg2=${yyyymmdd}18
  sbatch --output=${Rout18} --job-name=PlotRay18_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}

	#`````````````````````````````````````````
	# Collocate with Aeolus Mie-cloud

  arg1=MieCloud

		#job logfile names
  Mout00=OUTplotMIE00_${dd}
  Mout06=OUTplotMIE06_${dd}
  Mout12=OUTplotMIE12_${dd}
  Mout18=OUTplotMIE18_${dd}

  rm $Mout00 $Mout06 $Mout12 $Mout18

                # hour = 00 UTC
  arg2=${yyyymmdd}00
  sbatch --output=${Mout00} --job-name=PlotMie00_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}
                # hour = 06 UTC
  arg2=${yyyymmdd}06
  sbatch --output=${Mout06} --job-name=PlotMie06_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}
                # hour = 12 UTC
  arg2=${yyyymmdd}12
  sbatch --output=${Mout12} --job-name=PlotMie12_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}
                # hour = 18 UTC
  arg2=${yyyymmdd}18
  sbatch --output=${Mout18} --job-name=PlotMie18_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}

else
	# 'Aeolus' is NOT a listed dataset. Submit one job per 6-h analysis cycle (00, 06, 12, 18 UTC) = 4 total jobs.

  arg1=" "	# leave empty

  out00=OUTplot00_${dd}
  out06=OUTplot06_${dd}
  out12=OUTplot12_${dd}
  out18=OUTplot18_${dd}

  rm $out00 $out06 $out12 $out18

  		# hour = 00 UTC
  arg2=${yyyymmdd}00
  sbatch --output=${out00} --job-name=Plot00_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}
                # hour = 06 UTC
  arg2=${yyyymmdd}06
  sbatch --output=${out06} --job-name=Plot06_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}
                # hour = 12 UTC
  arg2=${yyyymmdd}12
  sbatch --output=${out12} --job-name=Plot12_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}
                # hour = 18 UTC
  arg2=${yyyymmdd}18
  sbatch --output=${out18} --job-name=Plot18_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7} ${run_job_script}

fi

#----------------------------------------------

  let idy=idy+1
done

#==============================================
# END Date LOOP
#==============================================
#==============================================

# END
##############################################################
