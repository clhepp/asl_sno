function drive_AIRS_CRIS_SNO_frmL1()
%
% Driver file for master script: make_AIRS_CRIS_SNO_frmL1.m
% Controlled as required by shell script: run_AIRS_CRIS_sno_frmL1.sh
% Submit batch jobs as follows: 
%        sbatch --array=0-N ./scripts/run_AIRS_CRIS_SNO_frmL1.sh
%
% C.Hepplewhite Jan 2018
% Modified from batch_AIRS_CRIS_SNO_frmL1.m

addpath /home/chepplew/gitLib/asl_sno/sno_maker

% Get the driver file of job dates - should be one calender month
 fh = fopen('driver_files/jobDates.drv','r');
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
cris     = 2;                %  1: NPP, 2: J01, 3: J02.

% Version reference to add to SNO save file, if none enter vers=''; Check with L1 source granules
% 'a2v4_ref': CRIS-2 hires; 'v20{a,d}': CRIS-1 & -2 from 2018, or '': CRIS-2 lowres, 'noaa'
vers = 'v20d';

% Options to control airs2cris conversion
opt1 = struct;
opt1.hapod    = 1;
opt1.nguard   = 2;
opt1.dvb      = 0.1;              % deconvolution frequency step
opt1.bfile    = '/tmp/bconv.mat'; %
opt1.inst_res = 'hires3';         % 'lowres' (default), 'hires1-3' or 'hi3to2'
opt1.user_res = 'midres';         % {'lowres','midres', 'hires'}
opt1.scorr    = 1;                % statistical correction (only for MSR & NSR)

% Call the master script
disp(['Starting master script ' char(datetime('now'))]);
disp(cris_res)

  make_AIRS_CRIS_sno_frmL1c(thisJob, cris_res, cris, opt1, vers);
  
fprintf(1, '*** Task run end %s\n', char(datetime('now')));
