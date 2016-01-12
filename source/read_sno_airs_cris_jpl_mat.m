function s = read_sno_airs_cris_jpl_mat(sdate1, sdate2, cchns)
%
% function read_sno_airs_cris_jpl_mat() reads data from the JPL SNO AIRS, CrIS mat files
%   and the radiances for a selected number of channels, specified by CrIS channel, 
%   and for the specified year and months. Unlike the sister function 'load_sno_airs_cris_jpl_mat'
%   it does not calculate statistics during load.
%
% Synopsis: read_sno_airs_cris_jpl_mat('date1','date2',[chan1...chan10]);
%           date1: first month as string: 'YYYY/MM/DD'
%           date2: last month as string: 'YYYY/MM/DD'
%           [chan1, ...]: numeric list of CrIS channels to load (max 10).
%           eg [403 499 737 884 905 998 1021 1297]
%
% Output: s. A structure containing all accummulated fields.
%
% Notes: If the CrIS channel is associated with a bad AIRS channel or has been
%        modified by AIRS L1C (clean and fill) then the next neares good channel
%        is substituted.
%
% Dependencies: i) AIRS good channel list; ii) nominal AIRS and CrIS frequency grids.
%    iii) fixed path and file name syntax for SNO files.
%
% Notes: i) No QA is applied. ii) time separation of SNO pairs from file is positive-only
%    so is recomputed here.
%

addpath /home/chepplew/gitLib/asl_sno/source
addpath /home/chepplew/gitLib/asl_sno/data

s = struct;
 
% Process the date strings
posYrs = [2002:2015];
posMns = [1:12];
% check dates are entered correctly:
if (length(sdate1) ~= 10 || length(sdate2) ~= 10) fprintf(1,'Error in dates\n'); exit; end
syr1 = sdate1(1:4);  smn1 = sdate1(6:7);  sdy1 = sdate1(9:10);
syr2 = sdate2(1:4);  smn2 = sdate2(6:7);  sdy2 = sdate2(9:10);
junk = ismember(posYrs, str2num(syr1)); if(isempty(~find(junk))) fprintf('invalid year\n'); end
junk = ismember(posMns, str2num(smn1)); if(isempty(~find(junk))) fprintf('invalid month\n'); end
junk = ismember([1:31], str2num(sdy1)); if(isempty(~find(junk))) fprintf('invalid day\n'); end
junk = ismember(posYrs, str2num(syr2)); if(isempty(~find(junk))) fprintf('invalid year\n'); end
junk = ismember(posMns, str2num(smn2)); if(isempty(~find(junk))) fprintf('invalid month\n'); end
junk = ismember([1:31], str2num(sdy2)); if(isempty(~find(junk))) fprintf('invalid day\n'); end
nyr1 = str2num(syr1); nmn1 = str2num(smn1);  ndy1 = str2num(sdy1);
nyr2 = str2num(syr2); nmn2 = str2num(smn2);  ndy2 = str2num(sdy2);
junk = sprintf('%4d/%02d/%02d',nyr1-1,12,31);
jdy1 = datenum(sdate1)-datenum(junk);  clear junk;           % needed for data directory
junk = sprintf('%4d/%02d/%02d',nyr2-1,12,31);
jdy2 = datenum(sdate2)-datenum(junk);  clear junk;           % needed for data directory

% Check channel numbers entered correctly
if(length(cchns) > 10 || length(cchns) < 1 ) fprintf(1,'Wrong number channels\n'); end
if(min(cchns) < 1 || max(cchns) > 1317 ) fprintf(1,'Wrong channel numbers used\n'); end

% get list of good AIRS channels (nig) to use, & bad (nib) to avoid
load('/home/chepplew/projects/airs/master_nig_01_2009.mat');   % nig [1 x 1535]
% nig = importdata('/home/strow/Work/Airs/good_chan_list');
junk = ismember([1:2378], nig);  nib = find(junk == 0);  clear junk;

xx=load('/asl/data/airs/airs_freq.mat'); fa=xx.freq; clear xx;
xx=load('cris_freq_nogrd.mat'); f_cris=xx.vchan; clear xx;  % 1305 chns (no guard chans)
xx=load('cris_freq_2grd.mat');  fc = xx.vchan; clear xx;    % 1317 chns (12 guard chans)
fd = f_cris([1:1037 1158:1305]);                            % prior knowledge from Howards decon routine

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
s.Wavs  = cWavs; 
s.sWavs = sWavs;

% ************* load up SNO data ********************
dp     = '/asl/s1/chepplew/projects/sno/airs_cris/JPL/v10_0_0/'; % or /JPL/standard/';
snoLst = dir(strcat(dp,'sno_airs_cris_*.mat'));
fprintf(1,'Found %d SNO files\n',numel(snoLst));

% subset range by date as requested:
dstart = datenum([nyr1 nmn1 ndy1]);
dlast  = datenum([nyr2 nmn2 ndy2]);
ifn1   = 1;
for i=1:numel(snoLst)
  junk = snoLst(i).name(15:22);
  thisdat = datenum( [str2num(junk(1:4)) str2num(junk(5:6)) str2num(junk(7:8))] );
  if(thisdat <= dstart) ifn1 = i; end
  if(thisdat <= dlast) ifn2 = i; end
end
fprintf(1,'Prcessing SNO files from: %s to %s\n',snoLst(ifn1).name, snoLst(ifn2).name);

% Load requested SNO data
s.td    = [];  s.arad = [;]; s.crad = [;]; s.drad = [;]; s.ctim = [];  s.atim = []; 
s.arlat = []; s.arlon = [];  s.dsn  = []; s.crlat = []; s.crlon = []; s.csolz = [];  
s.nSam  = []; s.alnfr = []; s.clnfr = [];  s.avrd = [;]; s.avra = [;]; s.avrc = [;]; 
s.sdra  = [;]; s.sdrc = [;]; s.sdrd = [;];s.cifv  = [];

for ifnum = ifn1:ifn2
  vars = whos('-file',strcat(dp,snoLst(ifnum).name));
  if( ismember('raDecon', {vars.name}) )              % was raDecon
    if(snoLst(ifnum).bytes > 1.0E4)
      g = load(strcat(dp,snoLst(ifnum).name));
      if(size(g.alat,1) < size(g.alat,2)) fprintf(1,'Unexpected row/column swapped\n'); end
      s.arad  = [s.arad, g.ra(achn,:)];                % [arad, [ra(achn,:); avaw]]; etc
      s.crad  = [s.crad, g.rc(cchn,:)];                % 1317 chns (12 guard chans)
      s.drad  = [s.drad, g.raDecon(dchn,:)];           %
      s.atim  = [s.atim,g.atime'];
      s.ctim  = [s.ctim,g.ctime'];
      s.arlat = [s.arlat,g.alat'];        s.arlon = [s.arlon,g.alon'];
      s.crlat = [s.crlat,g.clat'];        s.crlon = [s.crlon,g.clon'];
      s.cifv  = [s.cifv,g.cifov'];
      s.td    = [s.td,g.tdiff'];                       %
      s.dsn   = [s.dsn,g.dist'];
      s.csolz = [s.csolz,g.csolzen'];
      %%%s.alnfr = [s.alnfr, g.alandfrac'];   s.clnfr = [s.clnfr,g.clandfrac']; 
      s.nSam  = [s.nSam,single(size(g.ra,2))];

      s.avra  = [s.avra,nanmean(g.ra,2)];       s.sdra = [s.sdra,nanstd(g.ra,1,2)];  
      s.avrc  = [s.avrc,nanmean(g.rc,2)];       s.sdrc = [s.sdrc,nanstd(g.rc,1,2)];
      s.avrd  = [s.avrd,nanmean(g.raDecon,2)];  s.sdrd = [s.sdrd,nanstd(g.raDecon,1,2)];

    else
      fprintf('%d insufficient number samples\n',ifnum)
    end
  else
    fprintf('skip %s ',snoLst(ifnum).name(15:22) );
  end
  fprintf('.');
end    
fprintf('\n');
s.td2 = s.atim - s.ctim;              % tdiff from JPL are exclusively real positive -WRONG



%{
% plot options
addpath /asl/matlab2012/aslutil                         % simplemap
addpath /asl/matlib/plotutils                           % aslprint
abt   = real(rad2bt(s.Wavs,s.arad));   abtm  = real(rad2bt(fa,s.avra));
cbt   = real(rad2bt(s.Wavs,s.crad));   cbtm  = real(rad2bt(fc,s.avrc));
arr1  = ones(1,size(s.td2,2)); %junk1 = s.td2*86400; junk2 = s.dsn;

figure(1);clf;simplemap(s.arlat, s.arlon, s.td2/60); title('2013 AIRS CRIS stnd SNO dist (km)');
  % aslprint('./figs/2013_AirsCris_jplStnd_mapDelay.png');
figure(2);clf;semilogy(fa,s.avra(:,1),'b',fc,s.avrc(:,1),'m',fd,s.avrd(:,1),'c');grid on;
figure(2);clf;plot(fa,abtm(:,19),fc,cbtm(:,19));grid on;axis([-Inf Inf -Inf Inf]);
figure(3);clf;hist(s.arlat,180),xlim([-90 90]);
   xlabel('latitude bin');ylabel('population');title('AIRS CRIS stnd SNO All data Distribution');
   % aslprint('AirsCris_AllSno_lat_hist.png')
figure(3);clf;hist(abt(1,:),120); xlim([180 330]);grid on;
  xlabel('B.T. (K)');ylabel('population');title('2013 Stnd AIRS 900 wn population');
   % aslprint('./figs/2013_AC_stnd_Airs_900wn_hist.png')

%}

end
