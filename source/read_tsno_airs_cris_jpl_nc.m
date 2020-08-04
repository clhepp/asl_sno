function sd = read_tsno_airs_cris_jpl_nc(cchns)

%
% Function read_tsno_airs_cris_jpl_nc() reads in the Airs CrIS SNO tropical
%   subset netcdf data created by JPL from the CCAST CrIS files. 
%   cchns: Selected channels are supplied as input. 
%   Relevant fields are returned in the structure sd.
%
% Synpopsis: cchns: eg [403 499 737 884 905 998 1021 1297]

% Notes: reads in data from all available files. Maximum number of select channels
%   is 10 due to memory limitations.
%  
% Assumptions and Dependencies: Assumes fixed source directory at UMBC for the
%  input data files, and they have a fixed filename syntax. Assumes the AIRS has
%  been converted to the CrIS spectral grid (see analyze_tsno.m) and the data
%  are available in separate .mat files.
%  Dependent proedure:  read_netcdf_lls.m.
%  Dependent data: The Airs, CrIS and deconvolved spectral grids.
%
% Author: C. L. Hepplewhite. ASL UMBC/JCET.
%
% Version: Initial 02-May-2015.
%
cd /home/chepplew/projects/sno/tsno
addpath /home/strow/Matlab                        % read_netcdf_lls

% General initialize
load('~strow/Matlab/Airs/airs_f'); fa=f;
load('~strow/Matlab/Cris/Data/freq_cris');             % 1305 chns (no guard chans)
load('~chepplew/cris_f1317.mat');  fc = cfchan;        % 1317 chns (12 guard chans)
fd = fcris([1:1037 1158:1305]);                        % prior knowledge from Howards decon routine

% get list of good AIRS channels (nig) to use, & bad (nib) to avoid
load('/home/chepplew/projects/airs/master_nig_01_2009.mat');   % nig [1 x 1535]
% nig = importdata('/home/strow/Work/Airs/good_chan_list');
junk = ismember([1:2378], nig);  nib = find(junk == 0);  clear junk;

snoD   = '/asl/s1/strow/tsno/';
snoLst = dir(strcat(snoD,'tsno.l1c.2014-04*ccastumbc.nc'));
fprintf(1,'Found %d tsno files\n',numel( snoLst));
savD   = '/asl/s1/chepplew/projects/sno/tsno/';

% Check channel numbers entered correctly
if(length(cchns) > 10 || length(cchns) < 1 ) fprintf(1,'Wrong number channels\n'); end
if(min(cchns) < 1 || max(cchns) > 1317 ) fprintf(1,'Wrong channel numbers used\n'); end

% Screen the channel selection for AIRS bad channels and update if necessary:
cWavs = fc(cchns);
for i=1:numel(cWavs)
  tmp     = find(fa  > cWavs(i)-0.1,1);
  aref    = find(nig > tmp,1);  achn(i) = nig(aref);
end
cWavs = fa(achn);
for i=1:numel(cWavs)
  cchn(i) = find(fc  > cWavs(i)-0.25, 1);
  dchn(i) = find(fd  > cWavs(i)-0.25, 1);
end
for i=1:numel(cWavs) sWavs{i}  = sprintf('%6.2f',cWavs(i)); end

sd.td    = [];  sd.arad = [;]; sd.crad = [;]; sd.drad = [;]; sd.ctim = [];  sd.atim = []; 
sd.arlat = []; sd.arlon = [];   sd.dsn = []; sd.crlat = []; sd.crlon = []; sd.asolz = []; 
sd.csolz = [];  sd.nSam = []; sd.alnfr = []; sd.ocean = [];  sd.avrd = [;]; sd.avra = [;]; 
sd.avrc = [;]; sd.sdra = [;]; sd.sdrc = [;]; sd.sdrd  = [;]; sd.cifv = [];  sd.l1cp = [];  
sd.l1cr = [];

for ifn = 1:numel(snoLst);
  d = read_netcdf_lls( strcat(snoD,snoLst(ifn).name) );
  junk = snoLst(ifn).name(10:19);
  fmat = ['airs2cris.' junk '.mat'];
  g = load(strcat(savD,fmat));
  
  sd.arad  = [sd.arad,  d.L1bAIRS.rad(achn,:)];                % [arad, [ra(achn,:); avaw]]; etc
  sd.crad  = [sd.crad,  d.L1bCrIS.rad(cchn,:)];                % 1317 chns (12 guard chans)
  sd.drad  = [sd.drad,  g.drad(dchn,:)];                       % 
  sd.atim  = [sd.atim;  d.L1bAIRS.time];
  sd.ctim  = [sd.ctim;  d.L1bCrIS.time];
  sd.arlat = [sd.arlat; d.L1bAIRS.lat];        sd.arlon = [sd.arlon; d.L1bAIRS.lon];
  sd.crlat = [sd.crlat; d.L1bCrIS.lat];        sd.crlon = [sd.crlon; d.L1bCrIS.lon];
  sd.cifv  = [sd.cifv; d.L1bCrIS.fov];
  sd.td    = [sd.td; d.L1bAIRS.time - d.L1bCrIS.time];           % use [x,xd'] for JPL SNO
  sd.l1cp  = [sd.l1cp, d.L1bAIRS.L1cProc];
  sd.l1cr  = [sd.l1cr, d.L1bAIRS.L1cSynthReason];
  %%dsn   = [dsn,dist];
  sd.asolz = [sd.asolz; d.L1bAIRS.solzen];   sd.csolz = [sd.csolz; d.L1bCrIS.solzen];
  sd.alnfr = [sd.alnfr; d.L1bAIRS.landfrac]; sd.ocean = [sd.ocean; d.L1bCrIS.is_ocean]; % 
  sd.nSam  = [sd.nSam,size(d.L1bAIRS.rad,2)];
  sd.avra  = [sd.avra,nanmean(d.L1bAIRS.rad,2)];  sd.sdra = [sd.sdra,nanstd(d.L1bAIRS.rad,1,2)]; 
  sd.avrc  = [sd.avrc,nanmean(d.L1bCrIS.rad,2)];  sd.sdrc = [sd.sdrc,nanstd(d.L1bCrIS.rad,1,2)];
  sd.avrd  = [sd.avrd,nanmean(g.drad,2)];         sd.sdrd = [sd.sdrd,nanstd(g.drad,1,2)];
  
  fprintf(1,'.'); 
end
sd.achn = achn; sd.cchn = cchn; sd.dchn = dchn; sd.cWavs = cWavs;

sd

end
