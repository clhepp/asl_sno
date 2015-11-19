function mstr_airs_to_cris_mat(srcD, fPatt)

% SYNOPSIS: mstr_airs_to_cris(srcD, fPatt)
%
% DESCRIPTION: Master control script for running batch AIRS to CRIS deconvolution
%      routine 'proc_airs_to_cris_mat.m'
%
% PARAMETERS: srcD: the source directory of radiance spectra
%             fPatt: file name pattern for source files, not including
%             the .filetype.
%
% DEPENDENCIES: Motteler's deconvolution routines
%
% NOTES:  The source files are .mat files and can be written to.
%
%

%%cd /home/chepplew/projects/airs/AIRStoCRIS/batchJobs

addpath /asl/packages/ccast/source                %
addpath /asl/matlib/h4tools                       %  
addpath /asl/packages/airs_decon/source           %

disp(srcD)
disp(fPatt)
% check input parameters:
res = exist(srcD,'dir');
if(res ~= 7) fprintf(1,'Error: source directory %s not found\n',srcD); exit; end
clear res;
res = dir(strcat(srcD,'/',fPatt,'.mat'));
if(numel(res) < 1) fprintf(1,'Error: no files with this pattern found\n'); exit; end
fLst = res; clear fnams;

njobs = 2; clear sJobs;
for i=1:njobs sJobs{i} = [srcD fLst(i).name]; end

% create the list of files to process and to be used by proc_airs_to_cris_mat.m
FH = fopen('JobList.txt','w+');
fprintf(FH,'%s\r\n',sJobs{:});
fclose(FH);

% create the slurm batch control shell script
batch = './batch_Airs2Cris.slurm';
FH = fopen(batch,'w+');
fprintf(FH,'#!/bin/bash\n\n');

fprintf(FH,'#SBATCH --job-name=Ar2Cr\n');
fprintf(FH,'#SBATCH --partition=batch\n');
fprintf(FH,'#SBATCH --qos=medium\n');
fprintf(FH,'#SBATCH --account=pi_strow\n');
fprintf(FH,'#SBATCH -N1\n');
fprintf(FH,'#SBATCH --mem-per-cpu=12000\n');
fprintf(FH,'#SBATCH --cpus-per-task 1\n');
fprintf(FH,'#SBATCH --array=1-%d\n\n',njobs);

fprintf(FH,'MATLAB=''/usr/cluster/matlab/2014b/bin/matlab''\n');
fprintf(FH,'MATOPTS='' -nodisplay -nojvm -nosplash''\n\n');

junk = sprintf('srun $MATLAB $MATOPTS -r "proc_airs_to_cris_mat; exit"');
fprintf(FH,'%s\n',junk);

fclose(FH);


% now submit the slurm batch controller

disp('submit the slurm controller')
command = 'sbatch batch_Airs2Cris.slurm'
[status, cmdout] = system(command);
disp(status)

end
