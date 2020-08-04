#!/bin/bash

# sbatch options
#SBATCH --job-name=AIsno
#SBATCH --output=/home/chepplew/logs/sbatch/AIsno-slurm-%N.%A.%a.out
#SBATCH --error=/home/chepplew/logs/sbatch/AIsno-slurm-%N.%A.%a.err
#SBATCH --partition=batch
#SBATCH --qos=medium+
#SBATCH --account=pi_strow
#SBATCH --mem=11000
###SBATCH --constraint=lustre
#SBATCH -s
###SBATCH --exclude=cnode026,cnode239,cnode241,cnode242,,cnode260,cnode267,cnode204
#########SBATCH --array=1-21

# matlab options
####MATLAB='/usr/cluster/matlab/2016b/bin/matlab'
MATLAB='/usr/ebuild/software/MATLAB/2018b/bin/matlab'
MATOPTS=' -nodisplay -nojvm -nosplash'

# calling routine
echo "Calling drive_AIRS_IASI_SNO_frmL1.m"

$MATLAB $MATOPS -r "addpath('/home/chepplew/projects/sno/makeSNO'); drive_AIRS_IASI_SNO_frmL1; exit"

echo "Completed drive_AIRS_IASI_SNO_frmL1.m"
