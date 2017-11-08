function [s] = load_sno_airs_cris_jpl_mat(sdate1, sdate2, xchns)

% Load and prep stats for all channels from a 'SMALL' SNO set (see notes below)
%
% INPUT: sdate1 calendar date to start collection
%        e.g. sdate1='2013/01/01';
%        xchns: the channels corresponding to CrIS lo-res bands w/ no guard channels.
%               valid range 1:1305. or: LW: [1:713], MW: [714:1146], SW: [1157:1305]
%     
% OUTPUTS:
%        s: structure with fields of data ** BEFORE SUBSETTING **
%        s.clat,  s.clon        CrIS Latitude and Longitude
%        s.alat,  s.alon        AIRS Latitude and Longitude
%        s.ctim,  s.atim        CrIS and AIRS time
%        s.czols, s.azols       CrIS and AIRS solar zenith angle
%        s.cifov                CrIS IFOV [1..9]
%        s.td                   Time delay between Obs (secs)
%        s.dist                 Separation between Obs (km)
%        s.crad,  s.drad        Cris and AIRS raw radiance spectra.
%        s.fa, s.fc, s.fd       AIRS [2378] CrIS [1317] Common [1185] channels
%        s.ichns                The channels loaded in this session.
%        s.prf                  The quantil profiler to use.
%        s.fnames               The SNO filenames loaded.
%        s.inh, s.ish           Indexes of samples in the northern(southern) hemisphere.
%        s.idn, s.idd           Indexes of samples obtained during the day(night).
%
%
% VERSION: 1: Gulps subset of channels of SNO data for selected period.
%          Replaced dir() with unix(find...) which generates temporary file.
%
% NOTES: allow >20 GB on compute node.
%        default end date is dlast = 2013/12/31;
%
% C Hepplewhite. November 2015.
% CLH. Mar 2017. Update path to SNO data
%                changed grouping spectral channels to match CrIS LW, MW and SW.

clearvars s a wmstats g -except sdate igrp;

cd /home/chepplew/gitLib/asl_sno/run
addpath /home/chepplew/gitLib/asl_sno/source
addpath /home/chepplew/gitLib/asl_sno/data        % cris_freq_*.mat
addpath /home/chepplew/myLib/matlib/aslutil       % rad2bt.m
addpath /asl/packages/airs_decon/source           % hamm_app.m seq_match.m
addpath /home/chepplew/myLib/matlib/math          % remove_6sigma
addpath /home/chepplew/projects/cris              % cris_freq*.mat

s = struct;

% ---------------------------------
% Check the requested channel list:
if(length(xchns) > 1305 | min(xchns) < 1 | max(xchns) > 1305)
 error('Invalid channel selection, valid range 1:1305'); end
 
% --------------------------------
% Check & Process the date strings
posYrs = [2002:2017];
posMns = [1:12];
whos sdate1; disp([sdate1 ' to ' sdate2]); fprintf('\n');
try 
   D1 = datenum(sdate1,'yyyy/mm/dd');
   D2 = datenum(sdate2,'yyyy/mm/dd');
catch
   error('Incorrect Date Format')
   return
end
[nyr1 nmn1 ndy1] = datevec(D1);
[nyr2 nmn2 ndy2] = datevec(D2);
if(nyr1 ~= nyr2) error('Use same year only'); return; end
cyr1   = sdate1(1:4);     cmn1 = sdate1(6:7);     cdy1 = sdate1(9:10);
cyr2   = sdate2(1:4);     cmn1 = sdate2(6:7);     cdy1 = sdate2(9:10);
  junk = sprintf('%4d/%02d/%02d',nyr1-1,12,31);
jdy1   = datenum(sdate1)-datenum(junk);  clear junk;           % needed for data directory
  junk = sprintf('%4d/%02d/%02d',nyr2-1,12,31);
jdy2   = datenum(sdate2)-datenum(junk);  clear junk;           % needed for data directory

% ------------------------------------
% Get SNO files by dates as requested:
dp     = ['/home/chepplew/data/sno/airs_cris/JPL/standard/' cyr1 '/'];
snoLst = dir(strcat(dp, 'sno_airs_cris_*.mat'));

ifn1 = 1;             % default start with first file unless later.
for i=1:numel(snoLst)
  junk = regexp(snoLst(i).name,'[0-9]','match');
  junk = cell2mat(junk);
  thisdat = datenum(junk,'yyyymmdd');
  if(thisdat <= D1) ifn1 = i; end
  if(thisdat <= D2) ifn2 = i; end
end
fprintf(1,'Loading %d SNO files from: %s to %s\n',(ifn2-ifn1+1),snoLst(ifn1).name, ...
        snoLst(ifn2).name);

% --------------------------------------------------
% Get number of SNO pairs used to initialize arrays:
clear fs na;
for i=ifn1:1:ifn2;
  vars = whos('-file',strcat(dp,snoLst(i).name));
  if( ismember('raDecon', {vars.name}) )                   % AIRS->CrIS is present
    alat = load([dp snoLst(i).name],'alat');
    fs(i) = length(alat.alat);
    fprintf(1,'.')
  end
end
na = sum(fs); clear alat;
disp(['total SNO pairs to load: ' num2str(na)]);

% ------------------------------------------------------------
% get channel lists (need to know in advance which are stored)
load('/home/chepplew/projects/airs/master_nig_01_2009.mat');   % nig [1 x 1535]
% nig = importdata('/home/strow/Work/Airs/good_chan_list');
junk = ismember([1:2378], nig);  nib = find(junk == 0);  clear junk;

xx = load('/home/chepplew/projects/airs/airs_f.mat'); 
  f2645 = xx.fairs;  f2378 = xx.f; clear xx;
xx = load('cris_freq_nogrd.mat');     fc_ng = xx.vchan;       % 1305 chns (no guard chans)
xx = load('cris_freq_2grd.mat');      fc_2g = xx.vchan;       % 1317 chns (12 guard chans)
load('/home/chepplew/projects/sno/airs_cris/fa2c_x.mat');     % fa2c: 1185 frm HM decon routine

% rc (CrIS) radiances are supplied on 1317 channel grid (2 guard channels)
% Get the channels to load rc.
[xf, cchns] = intersect(fc_2g, fc_ng(xchns));

% ra2c (airs2cris) is supplied on 1185 channel grid
% Get the channels to load ra2c
[~, dchns] = intersect(fa2c, xf);                        %  maps fc onto fa2c

% the A2C channels are a subset of the CrIS channels
[iwant,~] = seq_match(fc_2g(cchns), fa2c(dchns));

 % map AIRS good chans to 2645 grid (not used here)
[zi zj] = seq_match(sort(f2378(nig)), f2645);
  
% Get the AIRS channels to load and record CrIS and A2C frequencies.
achns = [find(f2645 >= fc_ng(xchns(1)),1):find(f2645 >= fc_ng(xchns(end)),1)-1];
s.fc   = fc_2g(cchns(iwant))';  
s.fa   = f2645(achns);
s.fa2c = fa2c(dchns);

% for plotting at the end:
bands = [640, 900; 900, 1320; 1300 2560];                  % fcris(ichns(1):ichns(end));
bands = [640, 1100; 1200, 1760; 2150 2550];                % fcris(ichns(1):ichns(end));
bands = [650, 1095; 1210, 1615; 2182 2550];


% ------------- load data into memory --------------------- %
clear g;
n1 = 1;
s.ra    = single(zeros(length(achns),na)); 
s.rc    = single(zeros(length(iwant),na)); 
s.ra2c  = single(zeros(length(dchns),na)); 
s.ctim  = [na]; s.atim = [na]; s.cifov = single(zeros(1,na));
s.alat  = single(zeros(1,na)); s.alon  = single(zeros(1,na));  
s.clat  = single(zeros(1,na)); s.clon  = single(zeros(1,na)); 
s.csolz = single(zeros(1,na)); s.asolz = single(zeros(1,na)); 
s.cland = single(zeros(1,na)); s.aland = single(zeros(1,na));
s.td    = single(zeros(1,na)); s.dsn   = single(zeros(1,na));
for ifnum = ifn1:1:ifn2
  n2 = n1 + fs(ifnum)-1;
  vars = whos('-file',strcat(dp,snoLst(ifnum).name));
  if( ismember('raDecon', {vars.name}) )                   % AIRS->CrIS is present
     g = load([snoLst(ifnum).folder '/' snoLst(ifnum).name]);
     if( size(g.ra,2) == size(g.raDecon,2) )               % can get incomplete SNO pairs
        s.ra(:,n1:n2)   = g.raScaf(achns,:);               % the L1C radiances
        rc = single(hamm_app(double(g.rc(cchns(iwant),:))));
	rd = single(hamm_app(double(g.raDecon(dchns,:))));
	s.rc(:,n1:n2)   = rc;
	s.ra2c(:,n1:n2) = rd;
        s.csolz(n1:n2)  = g.csolzen;
	s.asolz(n1:n2)  = g.asolzen;
        s.td(n1:n2)     = g.tdiff;
        s.dsn(n1:n2)    = g.dist;
	s.clat(n1:n2)   = g.clat;      s.clon(n1:n2)  = g.clon;
	s.alat(n1:n2)   = g.alat;      s.alon(n1:n2)  = g.alon;
	s.ctim(n1:n2)   = g.ctime;     s.atim(n1:n2)  = g.atime;
	s.cifov(n1:n2)  = g.cifov;
	s.cland(n1:n2)  = g.clandfrac; s.aland(n1:n2) = g.alandfrac;
      else
        disp(['skip: ra .ne. raDecon. file: ' num2str(ifnum)]);
      end
  end
  fprintf(1,'%03d ',ifnum);
  n1 = n1 + fs(ifnum);
end
fprintf(1,'\n');
clear g rc rd;
[nx ny] = size(s.rc);
fprintf(1,'number of SNO pairs: %d\n', ny);

% ------------- Quality Control ---------------------
disp(['Finding outliers: stored in s.ing']);
acbias = single(s.rc - s.ra2c);
clear gx;
for i=1:length(dchns)
   n  = single(remove_6sigma(acbias(i,:)));
   nn = single(remove_6sigma(acbias(i,n)));
   gx(i).n = n(nn);
end
% Now find unique set of bad SNO samples
ux = [];
[~, psz] = size(acbias);
for i=1:length(dchns)
   ux = [ux setdiff(1:psz,gx(i).n)];
end
ux  = unique(ux);
s.ing = single(setdiff(1:psz,ux));
disp(['  ' num2str(numel(ux)) ' outliers located']);
clear gx n nn ux psz acbias;

% --------------  save radiances here --------------------- %
%{
clear fa;     fa = s.fa; 
clear fc;     fc = s.fc';
clear fa2c; fa2c = s.fa2c;
clear ra rc ra2c;
ra = s.ra(:,s.ing);
rc = s.rc(:,s.ing);
ra2c = s.ra2c(:,s.ing);
  whos fa fc fa2c ra rc ra2c;
savfn = '/home/chepplew/data/sno/airs_cris/radiances/2013Q1_sno_airs_cris_jpl_lw.mat';
save(savfn,'fa','fc','fa2c','ra','rc','ra2c','-v7.3');
%}

% --------------    subset night or day ---------------------
idn = find(s.csolz);                                                    % include all scenes for LW and MW bands
  s.idn = find(s.csolz > 95 & s.asolz > 95);
  s.idd = find(s.csolz < 90 & s.asolz < 90);
  fprintf(1,'Found %d day, %d night\n',numel(s.idd),numel(s.idn));
% --------------    subset NH or SH ---------------------
idn = find(s.clat);                                                    % include all scenes for LW and MW bands
  s.ish = find(s.clat < 0);
  s.inh = find(s.clat > 0);
  fprintf(1,'Found %d north, %d south\n',numel(s.inh),numel(s.ish));
% ----------------------------------------------------------------------- %


%{
% sanity plots
cbt = real(rad2bt(s.fc(s.ichns), s.crad));
dbt = real(rad2bt(s.fd(s.ichns), s.drad));
figure(10);clf;plot(s.fc(s.ichns), nanmean(cbt,2),'-', s.fc(s.ichns), nanmean(dbt,2),'-');
ich = 402   % (899.375)
btbins = [190:1:330];  btcens = [190.5:1:329.5];
pdf_c402 = histcounts(cbt(ich,:), btbins);
pdf_d402 = histcounts(dbt(ich,:), btbins);
figure(10);clf;plot(btcens, pdf_c402,'.-', btcens,pdf_d402,'.-');grid on;
dbins = [-20:1:20];   dcens = [-19.5:1:19.5];
pdf_bias = histcounts(dbt(ich,:) - cbt(402,:), dbins);
figure(10);clf;plot(dcens, pdf_bias, '.-');grid on;


% ------------- Get source data - start with preparation. ------------- %
clear x cc ccs;
unix(['cd ' dp '; find . -noleaf -maxdepth 1 -type f -name ''sno_airs_cris_2013*.mat'' -printf ''%P\n'' > /tmp/fn.txt;']);
fh = fopen('/tmp/fn.txt');
x  = fgetl(fh);
i  = 1;
while ischar(x)
   cc{i} = x;
   i = i + 1;
   x = fgetl(fh);
end
fclose(fh);
cc  = cellstr(cc);
ccs = sort(cc);
fprintf(1,'Found %d total SNO files\n',numel(ccs));


%}
