#!/bin/bash

# sbatch options
#SBATCH --job-name=ICsno
#SBATCH --output=/home/chepplew/logs/sbatch/ICsno-slurm-%N.%A.%a.out
#SBATCH --error=/home/chepplew/logs/sbatch/ICsno-slurm-%N.%A.%a.err
#SBATCH --partition=high_mem
#SBATCH --qos=medium+
#SBATCH --account=pi_strow
######## #SBATCH --constraint=lustre
#SBATCH --mem=10000
#SBATCH -s
######## #SBATCH --array=1-21
######## #SBATCH --exclude=cnode[202,203,204,212,213,236,240,250,260,282,284]

# matlab options
# #MATLAB='/usr/cluster/matlab/2016b/bin/matlab'
MATLAB='/usr/ebuild/software/MATLAB/2020a/bin/matlab'
MATOPTS=' -nodisplay -nojvm -nosplash'

# calling routine
echo "Calling drive_IASI_CRIS_SNO_frmL1.m"

$MATLAB $MATOPS -r "addpath('/home/chepplew/gitLib/asl_sno/sno_maker'); drive_IASI_CRIS_SNO_frmL1; exit"

echo "Completed drive_IASI_CRIS_SNO_frmL1.m"
