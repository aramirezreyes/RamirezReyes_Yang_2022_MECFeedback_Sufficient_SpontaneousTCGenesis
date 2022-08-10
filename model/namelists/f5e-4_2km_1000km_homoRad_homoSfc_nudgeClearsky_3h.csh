#! /bin/csh -f
#SBATCH --job-name nudgeClearsky5h
#SBATCH --qos=regular
#SBATCH --constraint=knl
#SBATCH --nodes=8
#SBATCH --tasks-per-node=64
#SBATCH --time=48:00:00
#SBATCH --license=project,SCRATCH
#SBATCH --mail-user=aramirezreyes@ucdavis.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module swap craype-haswell craype-mic-knl
source /etc/profile.d/modules.csh
module load cray-netcdf

#set loadedmodules=`module list`
echo Working directory is $cwd
set sst = 300
set run_name = f5e-4_2km_1000km_homoRad_homoSfc_nudgeClearsky_3h		# label for the run
set exp_home = $cwd:h
set exp_name = $exp_home:t
set sam_home = $cwd:h:h:h

#echo "${loadedmodules} "
echo "*** Running ${run_name} on ${HOSTNAME} ***"
date

# Tell me which nodes it is run on; for sending messages to help-hpc          
echo " "
echo This jobs runs on the following processors:
echo `cat ${SLURM_JOB_NODELIST}`
echo " "

#----------------------------------------------------------------
cd $exp_home

# define variables
 set data_dir    = $SCRATCH/sam3d/${exp_name}    # temporary directory for model workdir, output, etc
 set run_dir     = $data_dir/$run_name                  # tmp directory for current run
 set workdir     = $run_dir/workdir                   # where model is run and model output is produced
 set output_dir  = $run_dir/output                    # output directory will be created here
#
# # note the following init_cond's are overwritten later if reload_commands exists
 set init_cond   = ""
#
 set pathnames   = $exp_home/input/path_names         # path to file containing list of source paths
 set namelist    = $exp_home/input/namelists          # path to namelist file
 set diagtable   = $exp_home/input/diag_table        # path to field table (specifies tracers)
 set fieldtable  = $exp_home/input/field_table        # path to field table (specifies tracers)
 set execdir     = $run_dir/exe.sam                    # where code is compiled and executable is created
 set sourcedir   = $sam_home/src                      # path to directory containing model source code

 set casename    = RCE_Muller

# set time_stamp  = $fms_home/bin/time_stamp.csh       # generates string date for file name labels
#
# set ireload     = 1                                  # counter for resubmitting this run script
# set irun        = 1                                  # counter for multiple model submissions within this script
# #--------------------------------------------------------------------------------------------------------
#
cat > CaseName <<EOF
${casename}
EOF


# setup directory structure
 if ( ! -d $execdir ) mkdir -p $execdir

 if ( ! -e $workdir ) then
   mkdir -p $workdir $workdir/INPUT $workdir/RESTART
 else
   rm -rf $workdir
   mkdir -p $workdir $workdir/INPUT $workdir/RESTART
         echo "WARNING: Existing workdir $workdir removed."
 endif

 if ( ! -d $output_dir )  then
    mkdir -p $output_dir
    mkdir -p $output_dir/combine
    mkdir -p $output_dir/logfiles
    mkdir -p $output_dir/restart
 endif

 if ( ! -d $data_dir )  then
    mkdir -p $exp_dir
    mkdir -p $data_dir
    mkdir -p $data_dir/combine
 endif



# compile the model code and create executable
#
cp $sam_home/Makefile $execdir/
cp -r $sam_home/SRC   $execdir/
cp CaseName $execdir/               # overwrites the existing CaseName
cp -r $sam_home/RUNDATA $execdir/
cp $exp_home/srcmods/*       $execdir/SRC     # overwrite by using the modified source files. 
mkdir $execdir/$casename
cp $exp_home/input/*         $execdir/$casename
cp $exp_home/init_snd/snd${sst}         $execdir/$casename/snd
cd $execdir


#----------------------Specifying modules in SAM--------------
setenv SAM_SCR `pwd`
# ----------------------------------
# # specify scalar-advection directory in SRC
setenv ADV_DIR ADV_MPDATA
# specify SGS directory in SRC
setenv SGS_DIR SGS_TKE
# # specify radiation directory in SRC
setenv RAD_DIR RAD_CAM
# #setenv RAD_DIR RAD_RRTM
# specify microphysics directory in SRC
setenv MICRO_DIR MICRO_SAM1MOM
# #setenv MICRO_DIR MICRO_M2005
# #setenv MICRO_DIR MICRO_DRIZZLE
# ----------------------------------
# # specify (GNU) make utility
# #setenv GNUMAKE 'gnumake -j8'
 setenv GNUMAKE 'gmake -j8'
# #--------------------------------------------
# #--------------------------------------------

setenv SAM_DIR  `pwd`
setenv SAM_OBJ  $SAM_SCR/OBJ
setenv SAM_SRC  `pwd`/SRC

if !(-d $SAM_SCR) mkdir -p $SAM_SCR
if !(-d $SAM_SCR/OUT_2D) mkdir $SAM_SCR/OUT_2D
if !(-d $SAM_SCR/OUT_3D) mkdir $SAM_SCR/OUT_3D
if !(-d $SAM_SCR/OUT_MOMENTS) mkdir $SAM_SCR/OUT_MOMENTS
if !(-d $SAM_SCR/OUT_STAT) mkdir $SAM_SCR/OUT_STAT
if !(-d $SAM_SCR/OUT_MOVIES) mkdir $SAM_SCR/OUT_MOVIES
if !(-d $SAM_SCR/RESTART) mkdir $SAM_SCR/RESTART
if !(-d $SAM_OBJ) mkdir $SAM_OBJ

if !(-d OUT_2D) ln -s $SAM_SCR/OUT_2D  OUT_2D
if !(-d OUT_3D) ln -s $SAM_SCR/OUT_3D  OUT_3D
if !(-d OUT_MOMENTS) ln -s $SAM_SCR/OUT_MOMENTS OUT_MOMENTS
if !(-d OUT_STAT) ln -s $SAM_SCR/OUT_STAT  OUT_STAT
if !(-d OUT_MOVIES) ln -s $SAM_SCR/OUT_MOVIES  OUT_MOVIES
if !(-d RESTART) ln -s $SAM_SCR/RESTART RESTART
if !(-d OBJ) ln -s $SAM_OBJ  OBJ
#--------------------------------------------
#bloss: add "make clean" if MICRO or RAD options
#        have changed.
cat > MICRO_RAD_OPTIONS.new <<EOF
$ADV_DIR
$SGS_DIR
$MICRO_DIR
$RAD_DIR
EOF
if (-e $SAM_OBJ/MICRO_RAD_OPTIONS) then
  # use of cmp suggested by http://docs.hp.com/en/B2355-90046/ch14s03.html
  cmp -s $SAM_OBJ/MICRO_RAD_OPTIONS MICRO_RAD_OPTIONS.new
  if ($status != 0) then
    # the file has changed -- remove everything from SAM_OBJ
    #   so that we get a fresh compile of the model
    echo "MICRO or RAD option changed in Build.  Removing all object files from OBJ/"
    rm -f $SAM_OBJ/*
  endif
endif
# move the new options into $SAM_OBJ/MICRO_RAD_OPTIONS
mv -f MICRO_RAD_OPTIONS.new $SAM_OBJ/MICRO_RAD_OPTIONS
#--------------------------------------------


cd $SAM_OBJ

if ( !(-e Filepath) ) then
cat >! Filepath << EOF
$SAM_SRC
$SAM_SRC/$ADV_DIR
$SAM_SRC/$SGS_DIR
$SAM_SRC/$RAD_DIR
$SAM_SRC/$MICRO_DIR
$SAM_SRC/SIMULATORS
$SAM_SRC/TIMING
EOF
endif
$GNUMAKE -f $SAM_DIR/Makefile

### set parameter list
cd $execdir/$casename
cat > prm << EOF
&SGS_TKE
dosmagor = .true.
/

&MICRO_M2005
/

&PARAMETERS

caseid ='expt8.001',

nrestart = 0,
nrestart_skip = 119, !write restarts every 120 stat outputs (except at the end of the run when restarts are written regardless of nrestart_skip).

CEM = .true.,
OCEAN = .true.,
dosgs           = .true.,
dodamping       = .true.,
doupperbound    = .true.,
docloud         = .true.,
doprecip        = .true.,
dolongwave      = .true.,
doshortwave     = .true.,
dosurface       = .true.,
dosfchomo       = .true.,
doradhomo       = .true.,
dolargescale    = .false.,
doradforcing    = .false.,  !note that when this is true, you need two rads&radc files
dosfcforcing    = .true.,  !keep forcing the SST=300K
docoriolis      = .true.,
donudging_uv    = .true.,
donudging_tq    = .false.,
donudging_q_clearsky = .true.,
doperpetual     = .true.,
dofplane        = .true.,
fcor            = 4.97e-4, !0.253s-1 corresponds to 10 deg latitude 0.497e-4 20deg
tauls = 7200.,
tautqls = 84600.,
tauq_clearsky = 10800
water_loading_coeff = 0,
cold_pool_level = 1,
do_q_homo       = .false.,
do_T_homo       = .false.,
!z1 refers to lower boundary, z2 refers to upper boundary (meters)
homo_Tq_z1      = 10,
homo_Tq_z2      = 80,
homo_rad_z1     = 0,
homo_rad_z2     = 35000,

SFC_FLX_FXD    = .false.,
SFC_TAU_FXD    = .false.,

dx =    2000.,
dy =    2000.,
dt =    10.,


latitude0 = 20.,
longitude0 = -23.5,
nrad = 15, !frequency (in time step) of radiation computation (3-5 min should be ok)
day0=0.0,

nstop    = 432000,! (259200  30 days, 86400 = 10 days 432000 = 50 days,1296000=150) dt=10
nprint   = 360,
nstat    = 360,
nstatfrq = 360,

doSAMconditionals     = .false.,
dosatupdnconditionals = .false.,

output_sep      = .true., 
ncycle_max      = 8,
nsave2D         = 360,
nsave2Dstart    = 0,
nsave2Dend      = 99999999,
save2Dsep       = .true.,

nsave3D         = 720,
nsave3Dstart    = 0, 
nsave3Dend      = 99999999,
save3Dsep       = .true. /

EOF

cat > sfc << EOF
  day sst(K) H(W/m2) LE(W/m2) TAU(m2/s2)
     0.000    ${sst}.000     8.750     96.88      0.000
  1000.000    ${sst}.000     8.750     96.88      0.000
EOF

cd $execdir
##set NPROCS=`wc -l < $SLURM_JOB_NODELIST`

# mpirun ./ADV_MPDATA_SGS_TKE_RAD_CAM_MICRO_SAM1MOM > log.txt
srun  -n 512 -c 4 --cpu_bind=cores  ./ADV_MPDATA_SGS_TKE_RAD_CAM_MICRO_SAM1MOM > log.txt
