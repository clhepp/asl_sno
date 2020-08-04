#!/bin/bash

# sbatch options
#SBATCH --job-name=ACsno
#SBATCH --output=/home/chepplew/logs/sbatch/ACsno-slurm-%N.%A.%a.out
#SBATCH --error=/home/chepplew/logs/sbatch/ACsno-slurm-%N.%A.%a.err
#SBATCH --partition=high_mem
#SBATCH --qos=medium+
#SBATCH --account=pi_strow
#SBATCH --mem=16000
####SBATCH --constraint=lustre
#SBATCH -s
####SBATCH --exclude=cnode026,cnode239,cnode241,cnode242,,cnode260,cnode267,cnode204,cnode279,cnode282,cnode284
#########SBATCH --array=1-21

# matlab options
##MATLAB='/usr/cluster/matlab/2018b/bin/matlab'
MATLAB='/usr/ebuild/software/MATLAB/2020a/bin/matlab'
MATOPTS=' -nodisplay -nojvm -nosplash'

# calling routine
echo "Calling drive_AIRS_CRIS_SNO_frmL1.m"

$MATLAB $MATOPS -r "addpath('/home/chepplew/gitLib/asl_sno/sno_maker'); drive_AIRS_CRIS_SNO_frmL1; exit"

echo "Completed drive_AIRS_CRIS_SNO_frmL1.m"
