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

	# Set home path 'dir_home': This is the current working directory (where this script is located)
dir_home=/home/klukens/s4-cardinal/scripts/run/Wind_Collocation/collocation_and_plotting.AMV_4th_Int/
#dir_home=/data/users/klukens/collocation/collocation_and_plotting/
echo 'WORKING DIRECTORY = '${dir_home}

	# Set output path 'dir_out': This is where the collocation index files are saved
dir_out=/data/users/klukens/collocation/AMV_4th_intercomparison/
#dir_out=/data/users/klukens/collocation/longterm_anl/
echo 'OUTPUT DIRECTORY = '${dir_out}

	# Set path 'archive_parent': This is where the archive home directory is located
	#	Example: If the full path to /atmos-nc-dataset is /Full/Path/To/atmos-nc-dataset,
	#	         then 
	#	 	 archive_prefix = /Full/Path/To/
	#
	#	Note: /atmos-nc-dataset and /aeolus... directories should be in the same 'archive_parent' directory
	#
archive_parent=/scratch/

#----------------------------------------------
# Set datasets to collocate
#	Users must choose a driver dataset.
#	Users must choose at least one dependent dataset.
#	Users should only choose wind datasets from the following list, as written:
#		Aeolus		... Refers to NetCDF Aeolus winds (both Rayleigh-clear and Mie-cloudy) converted from ESA's Earth Explorer (EE) files.
#		Aircraft	... Refers to NetCDF Aircraft data converted from NCEP aircft prepBUFR files.
#		AMV_4th_Int	... Refers to NetCDF AMVs for 4th International Comparison
#		AMV_NCEP	... Refers to NetCDF AMVs converted from NCEP satwnd prepBUFR files.
#		Loon		... Refers to NetCDF Loon winds.
#		Radiosonde	... Refers to NetCDF Radiosonde winds converted from NCEP adpupa prepBUFR files.

	#set dataset names to use. += appends the names onto the end of the variable 'dataset_names'
	#	NOTES:
	#		1. The first name should be the driver dataset. All subsequent names should be dependent datasets.
	#		2. Make sure to add a comma at the end of each name <-- this is important for the collocation program.
	#		3. BEWARE! While the order of the DEPENDENT dataset names does not matter, the program assumes that 
	#		   the 'qc_flags' order corresponds to the order of 'dataset_names'.
	
#dataset_names=\"Aeolus,Aircraft,AMV_NCEP,Loon,Radiosonde\"
#dataset_names=\"Aeolus,Aircraft,Loon,Radiosonde\"
#dataset_names=\"Aircraft,Aeolus,Loon,Radiosonde\"
#dataset_names=\"Aeolus,Aircraft\"
dataset_names=\"AMV_4th_Int,Aeolus,Aircraft,Radiosonde\"

#AMV_center="BRZ"
#AMV_center="EUM"
#AMV_center="JMA"
#AMV_center="KMA"
AMV_center="NOA"
#AMV_center="NWC"

#----------------------------------------------
# Define collocation criteria
#       Each criterion should have the same number of elements, and that number should equal the number of DEPENDENT datasets

#dst_max=\"100.0,100.0,100.0,100.0,150.0\"          # collocation distance maximum in km ... horizontal distance
#prs_max=\"0.04,0.04,0.04,0.04,0.04\"           # colloation pressure difference maximum in log10(hPa) ... for Datasets with given pressures
#tim_max=\"60.0,60.0,60.0,60.0,90.0\"           # collocation time difference maximum in minutes
#hgt_max=\"1.0,1.0,1.0,1.0,1.0\"            # collocation height difference maximum in km ... vertical distance ... for Datasets given heights

dst_max=\"100.0,100.0,100.0,150.0\"          # collocation distance maximum in km ... horizontal distance
prs_max=\"0.04,0.04,0.04,0.04\"           # colloation pressure difference maximum in log10(hPa) ... for Datasets with given pressures
tim_max=\"60.0,60.0,60.0,90.0\"           # collocation time difference maximum in minutes
hgt_max=\"1.0,1.0,1.0,1.0\"

#dst_max=\"100.0\"          # collocation distance maximum in km ... horizontal distance
#prs_max=\"0.04\"           # colloation pressure difference maximum in log10(hPa) ... for Datasets with given pressures
#tim_max=\"60.0\"           # collocation time difference maximum in minutes
#hgt_max=\"1.0\"

#----------------------------------------------
# Set quality control (QC) flags for each dataset
#	QC flags should be set for each dataset named above and written as a single string, 
#	with each flag separated by a comma "," <-- this is important for the collocation program.
#
#	Set to 0 if no QC is to be applied
#	Set to 1 if QC is to be applied

#qc_flags=\"1,0,1,0,0\"
qc_flags=\"1,1,0,0\"
#qc_flags=\"1,0\"

#----------------------------------------------
# Set quality indicator (QI) in percent (%) for AMV dataset(s)
#	QI flags should be set for each AMV dataset named above and written as a single string, 
#	with each flag separated by a comma "," <-- this is important for the collocation program.
#
#	Set to value >= 0 and <= 100
#
#	NOTES: 
#		1. This variable must be set even if AMVs are not used.
#		2. If more than 1 AMV dataset is used, list each QI and separate with comma "," delimiter
#			For example: amv_qi_flags='80,60'
#		   Keep in mind that order matters! First (second) QI listed refers to first (second) AMV dataset listed in 'dataset_names', etc.

#amv_qi_flags=\"0\"
amv_qi_flags=\"60\"
#amv_qi_flags=\"80\"

#----------------------------------------------
# Set quality indicator (QI) choice(s)
#	'qi_choices' is a string
#
#	NO_FC (default value) --> QI without forecast
#	YES_FC                --> QI with forecast
#
#       NOTES:
#               1. This variable must be set even if AMVs are not used.
#               2. If more than 1 AMV dataset is used, list each QI and separate with comma "," delimiter
#                       For example: qi_choices='NO_FC,YES_FC'
#                  Keep in mind that order matters! First (second) QI listed refers to first (second) AMV dataset listed in 'dataset_names', etc.

qi_choices=\"NO_FC\"

#----------------------------------------------
# Number of collocations allowed per DRIVER observation

n_max=50

#----------------------------------------------
# Number of processors to use during parallelization

nproc=50
#nproc=100

#----------------------------------------------
# Set HPC account, partition, and runtime limit for jobs

account="star"

partition="s4"
timelimit="06:00:00"
#partition="serial"
#timelimit="09:00:00"

#----------------------------------------------
# For Aeolus dataset only: Choose if using reprocessed or non-reprocessed (original) dataset
#	'bline' and 'dtype' must be specificed.

	#orig = original = Aeolus dataset not reprocessed
#bline=orig
#dtype=original

	#reprocessed Aeolus dataset: use baseline number as bline (e.g., B10)
#bline=B10
bline=B11
dtype=reprocessed/2${bline}

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

dir_coll=${dir_home}src/           #location of collocation source code
echo 'COLLOCATION CODE DIRECTORY = '${dir_coll}

#----------------------------------------------
# Set name of main collocation code script (to pass as argument to run job script)

colloc_script=MAIN.match_driver_dependents.py

#----------------------------------------------
# Set name of run job script (to pass as argument to run job script)

run_job_script=run_collocation_code.job

#----------------------------------------------
# Load python 2.7 (for CODA)

cd

module load miniconda/2.7-base

__conda_setup="$('/opt/miniconda/2.7/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/miniconda/2.7/etc/profile.d/conda.sh" ]; then
        . "/opt/miniconda/2.7/etc/profile.d/conda.sh"
    else
        export PATH="/opt/miniconda/2.7/bin:$PATH"
    fi
fi
unset __conda_setup

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

#----------------------------------------------
# Activate .yml file to run python codes

#echo "... ACTIVATE YML -----"

  cd ${dir_coll}

  source activate bhoover-Universal_AMV_Matching &		# activates conda environment needed to run collocation code, and run in background (&)

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
#               2. QC applied to each dataset? 0=no, 1=yes. Separated by comma "," delimiter
#               3. If using AMVs, this is the lower limit of AMV QI in %, for QC
#               4. If using AMVs, this is the choice of QI type to use, for QC
#               5. Date in YYYYMMDDHH
#               6. Aeolus dataset type: for non-reprocessed data, use orig (original); for reprocessed files, use baseline number (example: B10)
#               7. Output directory
#		8. Collocation distance maximum
#	        9. Collocation pressure difference maximum
#	       10. Collocation time difference maximum
#	       11. Collocation height difference maximum
#	       12. Names of all datasets to use for collocation, separated by comma "," delimiter
#	       13. Archive parent path: path where home archive directory is located
#              14. Maximum number of matches allowed per data point
#              15. Number of processors to use during parallelization
#	       16. Abbreviation of wind-producing center

  #arg1: automatically set in loop below
  arg2=${qc_flags}
  arg3=${amv_qi_flags}
  arg4=${qi_choices}
  #arg5: automatically set in loop below
  arg6=${bline}
  arg7=${dir_out}
  arg8=${dst_max}
  arg9=${prs_max}
  arg10=${tim_max}
  arg11=${hgt_max}
  arg12=${dataset_names}
  arg13=${archive_parent}
  arg14=${n_max}
  arg15=${nproc}
  arg16=${AMV_center}


if [[ "$dataset_names" == *"Aeolus"* ]]; then		# check dataset_names (all of them) if 'Aeolus' is present
	# 'Aeolus' is a listed dataset. Submit one job per 6-h analysis cycle (00, 06, 12, 18 UTC) per Aeolus wind type (Rayleigh-clear or Mie-cloudy) = 8 total jobs.

	#`````````````````````````````````````````
        # Collocate with Aeolus Rayleigh-clear

  arg1=RayClear

		#job logfile names
  Rout00=OUTray00_${dd}
  Rout06=OUTray06_${dd}
  Rout12=OUTray12_${dd}
  Rout18=OUTray18_${dd}

  rm $Rout00 $Rout06 $Rout12 $Rout18

		# hour = 00 UTC
  arg5=${yyyymmdd}00
  sbatch --output=${Rout00} --job-name=CollRay00_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
		# hour = 06 UTC
  arg5=${yyyymmdd}06
  sbatch --output=${Rout06} --job-name=CollRay06_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
		# hour = 12 UTC
  arg5=${yyyymmdd}12
  sbatch --output=${Rout12} --job-name=CollRay12_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
		# hour = 18 UTC
  arg5=${yyyymmdd}18
  sbatch --output=${Rout18} --job-name=CollRay18_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}

	#`````````````````````````````````````````
	# Collocate with Aeolus Mie-cloud

  arg1=MieCloud

		#job logfile names
  Mout00=OUTmie00_${dd}
  Mout06=OUTmie06_${dd}
  Mout12=OUTmie12_${dd}
  Mout18=OUTmie18_${dd}

  rm $Mout00 $Mout06 $Mout12 $Mout18

                # hour = 00 UTC
  arg5=${yyyymmdd}00
  sbatch --output=${Mout00} --job-name=CollMie00_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 06 UTC
  arg5=${yyyymmdd}06
  sbatch --output=${Mout06} --job-name=CollMie06_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 12 UTC
  arg5=${yyyymmdd}12
  sbatch --output=${Mout12} --job-name=CollMie12_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 18 UTC
  arg5=${yyyymmdd}18
  sbatch --output=${Mout18} --job-name=CollMie18_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}

else
	# 'Aeolus' is NOT a listed dataset. Submit one job per 6-h analysis cycle (00, 06, 12, 18 UTC) = 4 total jobs.

  arg1=""	# leave empty

  out00=OUT00_${dd}
  out06=OUT06_${dd}
  out12=OUT12_${dd}
  out18=OUT18_${dd}

  rm $out00 $out06 $out12 $out18

  		# hour = 00 UTC
  arg5=${yyyymmdd}00
  sbatch --output=${out00} --job-name=Coll00_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 06 UTC
  arg5=${yyyymmdd}06
  sbatch --output=${out06} --job-name=Coll06_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 12 UTC
  arg5=${yyyymmdd}12
  sbatch --output=${out12} --job-name=Coll12_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 18 UTC
  arg5=${yyyymmdd}18
  sbatch --output=${out18} --job-name=Coll18_${dd} --account=${account} --partition=${partition} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}

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
