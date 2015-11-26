function s = read_sno_airs_iasi_jpl_mat(sdate1, sdate2, xchns)
%
% function read_sno_airs_iasi_jpl_mat() reads data from the JPL SNO AIRS, IASI mat files
%   and the radiances for a selected number of channels, specified by AIRS channel, 
%   and for the specified year and months. Unlike the sister function 
%  'load_sno_airs_iasi_jpl_mat' it does not calculate statistics during load.
%
% Synopsis: read_sno_airs_iasi_jpl_mat('date1','date2',[chan1...chan10]);
%           date1: first month as string: 'YYYY/MM/DD'
%           date2: last month as string: 'YYYY/MM/DD'
%           [chan1, ...]: numeric list of AIRS channels to load (max 10).
%           eg [403 499 737  884  905  998  1021 1297]
%           eg [759 902 1300 1608 1652 1789 1827 2203]
%
% Output: s. A structure containing all accummulated fields.
%
% Notes: If the requested channel is on the bad AIRS channel list or has been
%        modified by AIRS L1C (clean and fill) then the next nearest good channel
%        is substituted.
%
% Dependencies: i) AIRS good channel list; ii) nominal AIRS and IASI frequency grids.
%    iii) fixed path and file name syntax for SNO files: /asl/s1/chepplew/projects/sno/...
%       airs_iasi/JPL/sno_airs_iasi_YYYYMMDD.mat
%
% Notes: i) No QA is applied. ii) time separation of SNO pairs from file is positive-only
%    so is recomputed here. iii) the SNO files are monthly sets.
%
% VERSION: Nov 2015: changed dir() to unix(find...)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /home/chepplew/gitLib/asl_sno/run

addpath /home/chepplew/gitLib/asl_sno/source
addpath /asl/matlab2012/aslutil/                  % drdbt.m
addpath /asl/packages/ccast/source                % seq_match.m
addpath /asl/matlib/plotutils                     % aslprint.m

% Specify if SNO from JPL or ASL
SRC = 'JPL';

% check dates are entered correctly:
if (length(sdate1) ~= 10 || length(sdate2) ~= 10) fprintf(1,'Error in dates\n'); exit; end
dnum1 = datenum(sdate1,'yyyy/mm/dd'); dnum2 = datenum(sdate2,'yyyy/mm/dd');
if (dnum1 < datenum('2007/08/01','yyyy/mm/dd') ) fprintf(1,'First date is too early\n'); exit; end
if (dnum2 > datenum('2014/03/31','yyyy/mm/dd') ) fprintf(1,'Last date is too late\n'); exit; end

% Process the date strings
posYrs = [2007:2015];
posMns = [1:12];
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
if(length(xchns) > 10 || length(xchns) < 1 ) fprintf(1,'Wrong number channels\n'); end
if(min(xchns) < 1 || max(xchns) > 2378 ) fprintf(1,'Wrong channel numbers used\n'); end

% get list of good AIRS channels (nig) to use, & bad (nib) to avoid
load('/home/chepplew/projects/airs/master_nig_01_2009.mat');   % nig [1 x 1535]
% nig = importdata('/home/strow/Work/Airs/good_chan_list');
junk = ismember([1:2378], nig);  nib = find(junk == 0);  clear junk;

load('/asl/data/iremis/danz/iasi_f.mat');                  % fiasi [8461x1]
load('/asl/data/airs/airs_freq.mat'); fa=freq; clear freq; % fa    [2378x1]
[xa xi] = seq_match(fa, fiasi);

% Screen the channel selection for AIRS bad channels and update if necessary:
aWavs = fa(xchns);
for i=1:numel(aWavs)
  %tmp     = find(fa  > cWavs(i)-0.1,1);
  tmp = xchns(i);
  aref    = find(nig > tmp,1);  achn(i) = nig(aref);
end
aWavs = fa(achn);
for i=1:numel(aWavs)
  ichn(i) = find(fiasi  > aWavs(i)-0.125, 1);          % better matchup by going 0.125 low.
%  dchn(i) = find(fd  > cWavs(i)-0.125, 1);
end
for i=1:numel(aWavs) s.sWavs{i}  = sprintf('%6.2f',aWavs(i)); end

% ************* load up SNO data ********************
dp     = ['/asl/s1/chepplew/projects/sno/airs_iasi/' SRC '/'];      % standard/';
unix(['cd ' dp '; find . -noleaf -type f -name ''sno_airs_iasi_*.mat'' -printf ''%P\n'' > /tmp/fn.txt;']);
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
%fullfile(dp,ccs{i})  Note that i = 1:length(ccs)
  %snoLst = dir(strcat(dp,'sno_airs_cris_*.mat'));
fprintf(1,'Found %d total SNO files\n',numel(ccs));

% subset range by date as requested:
dstart = datenum([nyr1 nmn1 ndy1]);
dlast  = datenum([nyr2 nmn2 ndy2]);
for i=1:numel(ccs)
  junk = ccs{i}(15:22);
  thisdat = datenum( [str2num(junk(1:4)) str2num(junk(5:6)) str2num(junk(7:8))] );
  if(thisdat <= dstart) ifn1 = i; end
  if(thisdat <= dlast) ifn2 = i; end
end
fprintf(1,'Prcessing SNO files from: %s to %s\n',ccs{ifn1}, ccs{ifn2});

% Load requested SNO data
s.td    = [];  s.arad = [;]; s.irad = [;]; s.drad = [;]; s.itim = [];  s.atim = []; 
s.arlat = []; s.arlon = [];  s.dsn  = []; s.irlat = []; s.irlon = []; s.isolz = [];  
s.nSam  = []; s.alnfr = []; s.ilnfr = [];  s.avrd = [;]; s.avra = [;]; s.avri = [;]; 
s.sdra  = [;]; s.sdri = [;]; s.sdrd = [;];s.iifv  = [];  s.Wavs = []; s.iqual = [];
s.asolz = [];
d  = struct;  d.nSam = [];  d.avra = [];  d.avrd = [];
n  = struct;  n.nSam = [];  n.avra = [];  n.avrd = [];
np = struct; np.nSam = []; np.avra = []; np.avrd = [];
sp = struct; sp.nSam = []; sp.avra = []; sp.avrd = [];

for ifnum = ifn1:ifn2
  vars = whos('-file',strcat(dp,ccs{ifnum}));
  if( ismember('i2ra', {vars.name}) )                  % was raDecon
    %%%if(snoLst(ifnum).bytes > 1.0E4)
    if(vars(4).size(1) > 500)                           % check sample size using alat
      %%fprintf(1,'i: %d, alat: %d, ilat: %d\n', ifnum,vars(4).size(1),vars(15).size(1));
      if( vars(4).size(1) == vars(15).size(1) )         % ensure same no. obs
      g = load(strcat(dp,ccs{ifnum}));
      %if(size(g.alat,1) < size(g.alat,2)) fprintf(1,'Unexpected row/column swapped\n'); end
      %if(size(g.atime,1) ~= size(g.itime,1)) fprintf(1,'Unequal samples\n'); continue; end
      s.arad  = [s.arad, g.ra(achn,:)];                % [arad, [ra(achn,:); avaw]]; etc
      s.irad  = [s.irad, g.ri(ichn,:)];                % 1317 chns (12 guard chans)
      s.drad  = [s.drad, g.i2ra(achn,:)];              %
      s.atim  = [s.atim, g.atime'];
      s.itim  = [s.itim, g.itime'];
      s.arlat = [s.arlat,g.alat'];              s.arlon = [s.arlon, g.alon'];
      s.irlat = [s.irlat,g.ilat'];              s.irlon = [s.irlon, g.ilon'];
      s.td    = [s.td,   g.tdiff'];                     %
      s.dsn   = [s.dsn,  g.dist'];
      s.isolz = [s.isolz,g.isolzen'];           s.asolz = [s.asolz, g.asolzen'];
      s.nSam  = [s.nSam,single(size(g.ra,2))];
%      s.alnfr = [s.alnfr, g.alandfrac'];       s.ilnfr = [s.clnfr,g.ilandfrac']; 
%      s.iqual = [s.iqual, g.iqual'];
%      s.iifv  = [s.iifv, g.iifov'];
      s.avra  = [s.avra,nanmean(g.ra,2)];       s.sdra = [s.sdra,nanstd(g.ra,1,2)];  
      s.avri  = [s.avri,nanmean(g.ri,2)];       s.sdri = [s.sdri,nanstd(g.ri,1,2)];
      s.avrd  = [s.avrd,nanmean(g.i2ra,2)];     s.sdrd = [s.sdrd,nanstd(g.i2ra,1,2)];

      idy     = find(g.asolzen < 90);
      int     = find(g.asolzen >= 90);
      d.nSam  = [d.nSam, numel(idy)];
      n.nSam  = [n.nSam, numel(int)];
      d.avra  = [d.avra, nanmean(g.ra(:,idy),2)];
      n.avra  = [n.avra, nanmean(g.ra(:,int),2)];
      d.avrd  = [d.avrd, nanmean(g.i2ra(:,idy),2)];
      n.avrd  = [n.avrd, nanmean(g.i2ra(:,int),2)];

      inp     = find(g.alat > 50);
      ins     = find(g.alat < -50);
      np.nSam = [np.nSam, numel(inp)];
      sp.nSam = [sp.nSam, numel(ins)];
      np.avra = [np.avra, nanmean(g.ra(:,inp),2)];
      sp.avra = [sp.avra, nanmean(g.ra(:,ins),2)];
      np.avrd = [np.avrd, nanmean(g.i2ra(:,inp),2)];
      sp.avrd = [sp.avrd, nanmean(g.i2ra(:,ins),2)];
            
      clear g;
      end
    else
      fprintf('%d insufficient number samples\n',ifnum)
    end
  else
    fprintf('skip %s ',ccs{ifnum}(15:22) );
  end
  fprintf('.');
end    
fprintf('\n');
s.td2  = s.atim - s.itim;              % tdiff from JPL are all real positive :WRONG
s.Wavs = aWavs;                       % used by sno_quantile.m
fprintf(1,'Total number of SNO pairs: %d\n',numel(s.td2));

%{
% plot options
addpath /asl/matlib/plotutils                % aslprint.m
junk = ones(1,size(s.td2,2)); junk = s.td2*86400; junk = s.dsn;
figure(1);clf;simplemap(s.arlat, s.arlon, junk); title('AIRS IASI SNO delay map (secs)');
  % aslprint('./figs/AirsIasi_SNO_delayMap.png');
figure(2);clf;semilogy(fa,avra(:,1),'b',fc,avrc(:,1),'m',fd,avrd(:,1),'c');grid;
figure(1);clf;hist(junk,200),xlim([0 22]);xlabel('Separation (km)'); ylabel('population');
   title('AIRS IASI SNO separation histogram (km)');
  % aslprint('./figs/AirsIasi_SNO_distanceHist.png');
figure(1);clf;subplot(1,2,1);hist(s.arlat,1040);xlim([-78 -70]);
  subplot(1,2,2);hist(s.arlat,1040);xlim([70 78]);
   xlabel('latitude bin');ylabel('population');title('AIRS CRIS stnd SNO All data Distribution');


   % aslprint('AirsCris_AllSno_lat_hist.png')
figure(3);clf;hist(abt(4,:),100); xlabel('B.T. (K)');ylabel('population');
   title('AIRS All Sample SNO 1414 wn population'); % set(gca, 'YScale', 'log')
   % aslprint('AirsIasi_AllSno_1414wn_hist.png')

%}

% ********** SECTION II - Full spectrum all sample stats: ***********
%            ============================================
ratpm=0; ritpm=0; rdtpm=0; ratps=0; ritps=0; rdtps=0;

for i = 1:numel(s.nSam)
  ratpm = ratpm + s.avra(:,i).*s.nSam(i);    ritpm = ritpm + s.avri(:,i).*s.nSam(i);
  rdtpm = rdtpm + s.avrd(:,i).*s.nSam(i);
  ratps = ratps + ( s.sdra(:,i).*s.sdra(:,i) + s.avra(:,i).*s.avra(:,i) )*s.nSam(i);
  ritps = ritps + ( s.sdri(:,i).*s.sdri(:,i) + s.avri(:,i).*s.avri(:,i) )*s.nSam(i);
  rdtps = rdtps + ( s.sdrd(:,i).*s.sdrd(:,i) + s.avrd(:,i).*s.avrd(:,i) )*s.nSam(i);
end
gavra = ratpm/sum(s.nSam);   gavri = ritpm/sum(s.nSam);    gavrd = rdtpm/sum(s.nSam);
gsdra = real(sqrt( ratps/sum(s.nSam) - gavra.*gavra ));
gsdri = real(sqrt( ritps/sum(s.nSam) - gavri.*gavri ));
gsdrd = real(sqrt( rdtps/sum(s.nSam) - gavrd.*gavrd ));
gsera = gsdra/sqrt(sum(s.nSam));   gseri = gsdri/sqrt(sum(s.nSam));
gserd = gsdrd/sqrt(sum(s.nSam));

gavba = real(rad2bt(fa, gavra));    gavbi = real(rad2bt(fiasi, gavri));
gavbd = real(rad2bt(fa, gavrd));
bias  = gavbd-gavba;
gdrsd = 0.5*sqrt( (gsdra.*gsdra + gsdrd.*gsdrd) );
btm = 0.5*(gavba + gavbd);
  mdr = 1E-3*(1./drdbt(fa,btm) );
dbts  = mdr.*gdrsd;
dbte  = dbts/sqrt(sum(s.nSam));

%{
% Plotting Options
 figure(1);clf;h1=subplot(2,1,1);plot(fa,gavba,'b-',fiasi,gavbi,'r-',fa,gavbd,'g-');
   grid;axis([600 2700 210 270]);title('2015.6mo AIRS (b), IASI (r), I2A (g) SNO Mean');
   legend('Airs 1b','Iasi','iasi2airs','location','north'); ylabel('BT (K)');
  h2=subplot(2,1,2);plot(fa(nig),bias(nig),'m-',fa(nig),dbte(nig),'c-');grid;
  axis([600 2700 -0.8 0.8]);xlabel('wn (cm-1)');ylabel('Bias A-I (K)');
   legend('Airs - Iasi2Airs','location','north');
   linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp=get(h1,'position');
   set(h1,'position',[pp(1) pp(2)-pp(4)*0.1 pp(3) pp(4)*1.1])
   pp=get(h2,'position'); set(h2,'position',[pp(1) pp(2) pp(3) pp(4)*1.1]);
   % aslprint('./figs/2015_6mo_AirsIasi_SNO_BTbiasSpectrum.png')
  AI.fa = fa; AI.fi = fi; AI.gavba = gavba; AI.gavbi = gavbi; AI.gavbd = gavbd; 
  AI.dbte = dbte; 
  % save('SNO_AI_201X_BTspect.mat','AI');
%}

% --------------------------------
%          Day night
% --------------------------------
n.ratpm = 0; n.rdtpm = 0; n.mbias = 0;
for i = 1:numel(n.nSam)
  n.ratpm = n.ratpm + n.avra(:,i).*n.nSam(i);
  n.rdtpm = n.rdtpm + n.avrd(:,i).*n.nSam(i);
end
n.ratpm = n.ratpm/sum(n.nSam);
n.rdtpm = n.rdtpm/sum(n.nSam);
n.abtm  = real(rad2bt(fa,n.ratpm));
n.dbtm  = real(rad2bt(fa,n.rdtpm));
n.mbias = n.dbtm - n.abtm;

d.ratpm = 0; d.rdtpm = 0; d.mbias = 0;
for i = 1:numel(d.nSam)
  d.ratpm = d.ratpm + d.avra(:,i).*d.nSam(i);
  d.rdtpm = d.rdtpm + d.avrd(:,i).*d.nSam(i);
end
d.ratpm = d.ratpm/sum(d.nSam);
d.rdtpm = d.rdtpm/sum(d.nSam);
d.abtm  = real(rad2bt(fa,d.ratpm));
d.dbtm  = real(rad2bt(fa,d.rdtpm));
d.mbias = d.dbtm - d.abtm;

fprintf(1,'Number day samples: %d, night: %d\n',sum(d.nSam),sum(n.nSam));

% -----------------------------------
%           north/south polar
% -----------------------------------
np.ratpm = 0; np.rdtpm = 0; np.mbias = 0;
for i = 1:numel(np.nSam)
  np.ratpm = np.ratpm + np.avra(:,i).*np.nSam(i);
  np.rdtpm = np.rdtpm + np.avrd(:,i).*np.nSam(i);
end
np.ratpm = np.ratpm/sum(np.nSam);
np.rdtpm = np.rdtpm/sum(np.nSam);
np.abtm  = real(rad2bt(fa,np.ratpm));
np.dbtm  = real(rad2bt(fa,np.rdtpm));
np.mbias = np.dbtm - np.abtm;

sp.ratpm = 0; sp.rdtpm = 0; sp.mbias = 0;
for i = 1:numel(sp.nSam)
  sp.ratpm = sp.ratpm + sp.avra(:,i).*sp.nSam(i);
  sp.rdtpm = sp.rdtpm + sp.avrd(:,i).*sp.nSam(i);
end
sp.ratpm = sp.ratpm/sum(sp.nSam);
sp.rdtpm = sp.rdtpm/sum(sp.nSam);
sp.abtm  = real(rad2bt(fa,sp.ratpm));
sp.dbtm  = real(rad2bt(fa,sp.rdtpm));
sp.mbias = sp.dbtm - sp.abtm;


%{
bands = [640, 1650; 2170, 2680]
figure(2);clf;h1=subplot(2,1,1);
  plot(fa,n.abtm,'b-',fa,n.dbtm,'g-');grid on;axis([bands(1,:) 210 260]);
  title('AIRS (b) IASI (g) SNO LW');ylabel('BT K');
  h2=subplot(2,1,2);plot(fa,n.mbias,'m-',fa,d.mbias,'c-');
    grid on;axis([bands(1,:) -1 1]);legend('night','day');ylabel('Bias K');xlabel('wn cm-1');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp=get(h1,'position');
  set(h1,'position',[pp(1) pp(2)-pp(4)*0.1 pp(3) pp(4)*1.1])
  pp=get(h2,'position'); set(h2,'position',[pp(1) pp(2) pp(3) pp(4)*1.1]);
  %aslprint('./figs/AI_SNO_JPL_DN_biasBT_LW.png');
%}
end