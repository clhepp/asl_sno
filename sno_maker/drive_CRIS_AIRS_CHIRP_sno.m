function drive_CRIS_AIRS_CHIRP_sno()
%
% Driver file for master script: make_AIRS_CRIS_SNO_frmL1.m
% Controlled as required by shell script: run_AIRS_CRIS_sno_frmL1.sh
% Submit batch jobs as follows: 
%        sbatch --array=0-N ./scripts/run_AIRS_CRIS_SNO_frmL1.sh
%
% C.Hepplewhite Jan 2018
% Modified from batch_AIRS_CRIS_SNO_frmL1.m

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

thisJob = cell2mat(fdates{1}(sindex+1));
  disp(['this job: ' thisJob]);

% COnfigure job settings and options
cris_res = {'high','mid'};   % 'low','mid','high'. {SDR source,  SNO product}
cris     = 1;                %  1: NPP, 2: J01

% Version reference to add to SNO save file, if none enter sver=''; Check with L1 source granules
vers = 'v01a';      % 'v01a': CHIRP v1.

% Options to control airs2cris conversion
opts = struct;
opts.hapod    = 1;
opts.nguard   = 2;
opts.dvb      = 0.1;              % deconvolution frequency step
opts.bfile    = '/tmp/bconv.mat'; %
opts.inst_res = 'hires3';         % 'lowres' (default), 'hires1-3' or 'hi3to2'
opts.user_res = 'midres';         % {'lowres','midres', 'hires'}
opts.scorr    = 0;                % statistical correction (only for MSR & NSR)

% Call the master script
disp(['Starting master script ' char(datetime('now'))]);
disp(cris_res)

  make_CRIS_AIRS_CHIRP_sno(thisJob, cris_res, cris, opts, vers);
  
fprintf(1, '*** Task run end %s\n', char(datetime('now')));
