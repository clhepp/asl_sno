function drive_IASI_AIRS_SNO_frmL1()
%
% Driver file for master script: make_IASI_CRIS_SNO_frmL1.m
% Controlled as required by shell script: run_IASI_CRIS_sno_frmL1.sh
% Submit batch jobs as follows: 
%        sbatch --array=0-N ./scripts/run_IASI_CRIS_sno_frmL1.sh
%
% C.Hepplewhite Jan 2018
% Modified from batch_IASI_CRIS_SNO_frmL1.m

addpath /home/chepplew/projects/sno/makeSNO

cd /home/chepplew/gitLib/asl_sno/sno_maker/

% Get the driver file of job dates - should be one calender month
 fh = fopen('./driver_files/jobDates.drv','r');
 fdates = textscan(fh,'%s');                 % cell array
 fclose(fh);
 njobs = numel(fdates{1});
 disp(['number of jobs: ' num2str(njobs)]);

% Get the slurm array job and assign current job.
sindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));   % 0-19999
%sindx=0

thisJob = cell2mat(fdates{1}(sindex+1));
  disp(['this job: ' thisJob]);

% COnfigure job settings and options {SDR source, SNO output}
cris_res = {'high','mid'};
iasi     = 1;
cris     = 1;
vers     = 'v20d';
% Options to control iasi2cris conversion (must not be called opts, conflict with CCAST SDR)
opt1 = struct;
opt1.hapod    = 1;
opt1.nguard   = 2;
opt1.user_res = 'midres';

% Call the master script
disp(['Starting master script ' char(datetime('now'))]);

  make_IASI_CRIS_SNO_frmL1(thisJob, cris_res, iasi, cris, opt1, vers);
  
fprintf(1, '*** Task run end %s\n', char(datetime('now')));

