function s = read_sno_airs_iasi_jpl_mat(sdate1, sdate2, xchns)
%
% function read_sno_airs_iasi_jpl_mat() reads data from the JPL SNO AIRS, IASI mat files
%   and the radiances for a selected number of channels, specified by AIRS channel, 
%   and for the specified year and months. Unlike the sister function 
%  'load_sno_airs_iasi_jpl_mat' it does not calculate statistics during load.
%
% Synopsis: read_sno_airs_iasi_jpl_mat('date1','date2',[chan1...chan10]);
%           sdate1: first month as string: 'YYYY/MM/DD', typically: '2011/01/02';
%           sdate2: last month as string: 'YYYY/MM/DD',  typically: '2014/03/01';
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
addpath /home/chepplew/gitLib/asl_sno/data        % cris_freq*.mat
addpath /asl/matlab2012/aslutil/                  % drdbt.m
addpath /asl/packages/ccast/source                % seq_match.m
addpath /asl/matlib/plotutils                     % aslprint.m
addpath /asl/packages/airs_decon/source           % hamm_app

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

xx=load('cris_freq_2grd.mat');  fcris = xx.vchan; clear xx;    % 1317 chns (12 guard chans)
load('/asl/data/iremis/danz/iasi_f.mat');                  % fiasi [8461x1]
load('/asl/data/airs/airs_freq.mat'); fa=freq; clear freq; % fa    [2378x1]
[xa xi]   = seq_match(sort(fa), fiasi);
%%fa2c      = fcris([1:1037 1158:1305]);
xx=load('/asl/s1/chepplew/projects/sno/airs_iasi/JPL/airs2cris_CAF_20101101.mat'); fa2c=xx.a2cfrq; clear xx;
[xa2c xc] = seq_match(fa2c, fcris);

bands = [640, 1100; 1250, 1650; 2170, 2680];

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
  cchn(i) = find(fcris  > aWavs(i)-0.125, 1);
end
for i=1:numel(aWavs) s.sWavs{i}  = sprintf('%6.2f',aWavs(i)); end


% ---------------------------------------------------
%               load up SNO data 
% ---------------------------------------------------
dp     = ['/asl/s1/chepplew/projects/sno/airs_iasi/' SRC '/'];      % standard/';
unix(['cd ' dp '; find . -maxdepth 1 -noleaf -type f -name ''sno_airs_iasi_*.mat'' -printf ''%P\n'' > /tmp/fn.txt;']);
fh = fopen('/tmp/fn.txt');
x  = fgetl(fh);
i  = 1;
while ischar(x)
   cc{i} = x;
   i = i + 1;
   x = fgetl(fh);
end
fclose(fh);
cc  = unique(cellstr(cc));
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
s = struct;
s.td    = [];  s.arad = [;]; s.irad = [;]; s.a2rc = [;]; s.i2ra = []; s.i2rc = [];
s.itim  = [];  s.atim = []; 
s.arlat = []; s.arlon = [];   s.dsn  = []; s.irlat = []; s.irlon = []; s.isolz = [];  
s.nSam  = []; s.alnfr = [];  s.ilnfr = [];  
s.avra  = [;]; s.avri = [;]; s.avri2c = [;]; s.avra2c = [];
s.sdra  = [;]; s.sdri = [;]; s.sdri2c = [];  s.sdra2c = [];
s.iifv  = [];  s.Wavs = []; s.iqual = [];
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
      s.i2ra  = [s.i2ra, g.i2ra(achn,:)];              %
      s.a2rc  = [s.a2rc, g.ra2c(cchn,:)];
      s.i2rc  = [s.i2rc, g.i2rc(cchn,:)];
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
      s.avra    = [s.avra,nanmean(g.ra,2)];       s.sdra   = [s.sdra,nanstd(g.ra,1,2)];  
      s.avri    = [s.avri,nanmean(g.ri,2)];       s.sdri   = [s.sdri,nanstd(g.ri,1,2)];
      s.avri2c  = [s.avri2c,nanmean(g.i2rc,2)];   s.sdri2c = [s.sdri2c,nanstd(g.i2rc,1,2)];
% ned quality control on ra2c
        junk = real(g.ra2c);
        inZ  = find(junk < -999);  junk(inZ) = NaN;
      s.avra2c  = [s.avra2c,nanmean(junk,2)];   s.sdra2c = [s.sdra2c,nanstd(junk,1,2)];

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
s.Wavs = aWavs;                        % used by sno_quantile.m
fprintf(1,'Total number of SNO pairs: %d\n',numel(s.td2));

%{
% Sanity Check & plot options
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

K=3;
btairs  = real( rad2bt(fa,s.avra(:,K)) );
   junk = single(hamm_app(double(s.avri(:,K))) );
btiasi  = real( rad2bt(fiasi,junk));
bta2c   = real( rad2bt(fa2c,s.avra2c(:,K)) );
   junk = single( hamm_app(double(s.avri2c(:,K))) );
bti2c   = real(rad2bt(fcris,junk));
  whos btairs btiasi bta2c bti2c
bands = [640, 1100; 1200, 1620; 2170, 2400];
figure(2);clf;plot(fa,btairs,'b-',fiasi,btiasi,'g-');grid on;axis([bands(1,:) 210 250]);
  title('2010.12 SNO AIRS (b) IASI (g) mean'); xlabel('wn cm^{-1}');ylabel('BT K');
  %saveas(gcf,'./figs/201011_AI_SNO_AirsIasi_meanBT_MW.png','png');
figure(3);clf;
 h1=subplot(2,1,1);plot(fa2c,bta2c,'b-',fcris,bti2c,'g-');grid on; axis([bands(1,:) 210 250]);
  title('2010.12 Airs.Iasi.SNO A2C (b) I2C (g) BT mean');ylabel('BT K');
 h2=subplot(2,1,2);plot(fa2c,bta2c - bti2c(xc),'m-');grid on;axis([bands(1,:) -1 1]);
  xlabel('wn cm^{-1}');ylabel('BT K');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])  
  %saveas(gcf,'./figs/201011_AI_SNO_A2C_I2C_meanBT_SW.png','png');

%}

% ********** SECTION II - Full spectrum all sample stats: ***********
%            ============================================
s.avhi2c = s.avri2c;  % single( hamm_app(double(s.avri2c)) );              % Don't 4get to hamming the IASI
s.avha2c = s.avra2c;  % single( hamm_app(double(s.avra2c)) );
r=struct;
r.ratpm=0; r.ritpm=0; r.ra2cm=0; r.ri2cm=0;
r.ratps=0; r.ritps=0; r.ra2cs=0; r.ri2cs=0;

for i = 1:numel(s.nSam)
  r.ratpm = r.ratpm + s.avra(:,i).*s.nSam(i);    r.ritpm = r.ritpm + s.avri(:,i).*s.nSam(i);
  r.ra2cm = r.ra2cm + s.avha2c(:,i).*s.nSam(i);  r.ri2cm = r.ri2cm + s.avhi2c(:,i).*s.nSam(i);
  r.ratps = r.ratps + ( s.sdra(:,i).*s.sdra(:,i) + s.avra(:,i).*s.avra(:,i) )*s.nSam(i);
  r.ritps = r.ritps + ( s.sdri(:,i).*s.sdri(:,i) + s.avri(:,i).*s.avri(:,i) )*s.nSam(i);
  r.ra2cs = r.ra2cs + ( s.sdra2c(:,i).*s.sdra2c(:,i) + s.avha2c(:,i).*s.avha2c(:,i) )*s.nSam(i);
  r.ri2cs = r.ri2cs + ( s.sdri2c(:,i).*s.sdri2c(:,i) + s.avhi2c(:,i).*s.avhi2c(:,i) )*s.nSam(i);
end
r.gavra    = r.ratpm/sum(s.nSam);            r.gavri   = r.ritpm/sum(s.nSam);    
r.gavra2c  = r.ra2cm/sum(s.nSam);            r.gavri2c = r.ri2cm/sum(s.nSam);
r.gsdra    = real(sqrt( r.ratps/sum(s.nSam) - r.gavra.*r.gavra ));
r.gsdri    = real(sqrt( r.ritps/sum(s.nSam) - r.gavri.*r.gavri ));
r.gsdra2c  = real(sqrt( r.ra2cs/sum(s.nSam) - r.gavra2c.*r.gavra2c ));
r.gsdri2c  = real(sqrt( r.ri2cs/sum(s.nSam) - r.gavri2c.*r.gavri2c ));
r.gsera    = r.gsdra/sqrt(sum(s.nSam));      r.gseri   = r.gsdri/sqrt(sum(s.nSam));
r.gsera2c  = r.gsdra2c/sqrt(sum(s.nSam));    r.gseri2c = r.gsdri2c/sqrt(sum(s.nSam));

r.gavba    = real(rad2bt(fa, r.gavra));      r.gavbi   = real(rad2bt(fiasi, r.gavri));
r.gavba2c  = real(rad2bt(fa2c, r.gavra2c));  r.gavbi2c = real(rad2bt(fcris, r.gavri2c));
r.bias     = r.gavba2c - r.gavbi2c;
r.gdrsd    = 0.5*sqrt( (r.gsdra.*r.gsdra + r.gsdrd.*r.gsdrd) );
btm = 0.5*(gavba + gavbd);
  mdr = 1E-3*(1./drdbt(fa,btm) );
dbts  = mdr.*gdrsd;
dbte  = dbts/sqrt(sum(s.nSam));

%{
save('AI_sno_A2C_I2C_plot_data.mat','s','r','fa','fcris','fiasi','fa2c','bands','xc');
% Sanity Check & Plotting Options
 figure(1);clf;h1=subplot(2,1,1);
   plot(fa,r.gavba,'b-',fiasi,r.gavbi,'r-',fa2c,r.gavba2c,'g-',fcris,r.gavbi2c,'c-');axis([bands(1,:) 210 255])
   grid on;title('2011 AI SNO AIRS (b), IASI (r), A2C (g) I2C (c) Mean');
   ylabel('BT (K)'); %%%legend('Airs 1b','Iasi','iasi2airs','location','north'); 
  h2=subplot(2,1,2);plot(fa2c,r.gavba2c - r.gavbi2c(xc),'.-');grid on;axis([bands(1,:) -0.8 0.8]);
   xlabel('wn (cm-1)');ylabel('Bias K');legend('A2C - I2C','location','north');
   linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp1=get(h1,'position');
   set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
   pp2=get(h2,'position'); set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1]);
   % aslprint('./figs/AI_SNO_A2C_I2C_BTbias_LW.png')

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
