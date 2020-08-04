function drive_AIRS_IASI_SNO_frmL1()
%
% Driver file for master script: make_AIRS_IASI_SNO_frmL1.m
% Controlled as required by shell script: run_AIRS_IASI_sno_frmL1.sh
% Submit batch jobs as follows: 
%        sbatch --array=0-N ./scripts/run_AIRS_IASI_sno_frmL1.sh
%
% C.Hepplewhite Jan 2018
% Modified from batch_AIRS_IASI_SNO_frmL1.m

addpath /home/chepplew/projects/sno/makeSNO

cd /home/chepplew/projects/sno/makeSNO/batchJobs

% Get the driver file of job dates - should be one calender month
 fh = fopen('jobDates.drv','r');
 fdates = textscan(fh,'%s');                 % cell array
 fclose(fh);
 njobs = numel(fdates{1});
 disp(['number of jobs: ' num2str(njobs)]);

% Get the slurm array job and assign current job.
sindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));   % 0-19999
%sindx=0
disp(['sindex: ' num2str(sindex)])

thisJob = cell2mat(fdates{1}(sindex+1));
  disp(['this job: ' thisJob]);

% COnfigure job settings and options
airs     = 'L1C';
iasi     = 2;
% Options to control iasi2airs conversion
opts = struct;
opts.sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
%opts.hapod   = 1;
%opts.resmode = 'lowres';
%opts.nguard  = 2;

% Call the master script
disp(['Starting master script ' char(datetime('now'))]);

  make_AIRS_IASI_sno_frmL1(thisJob, airs, iasi, opts);
  
fprintf(1, '*** Task run end %s\n', char(datetime('now')));
