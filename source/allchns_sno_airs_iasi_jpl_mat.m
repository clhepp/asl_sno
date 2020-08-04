function [s wmstats] = allchns_sno_airs_iasi_jpl_mat(sdate, igrp)

% Load and prep stats for all channels from a SMALL SNO set

% INPUT: sdate calendar date to start collection
%        e.g. sdate='2011/01/02';
%        igrp - the group number for the channels in multiples of 200.
%           selected to cover the common grid channels (1185)
%           valid values: 1, 2, 3.
%
% OUTPUT: wmstats - see notes at bottom for details
%
% VERSION: 1: tuned to accommodate channels in one of the CrIS equivalent bands
%         of SNO data.
%         Replaced dir() with unix(find...) which generates temporary file.
%
% NOTES: allow 10 GB on compute node.
%
% C Hepplewhite. November 2015.
% CLH. Mar2017. Updated paths to SNO data.
%               common refinemnets with other two 'allchns_sno_X_Y.m' scripts.

%%clearvars s a wmstats g -except sdate igrp;

cd /home/chepplew/gitLib/asl_sno/run
addpath /home/chepplew/gitLib/asl_sno/source
addpath /home/chepplew/gitLib/asl_sno/data        % cris_freq_*.mat
addpath /home/chepplew/myLib/matlib/aslutil       % rad2bt.m
addpath /home/strow/Git/breno_matlab/Math         % Math_bin.m
addpath /asl/matlab2012/aslutil/                  % drdbt.m
addpath /asl/packages/ccast/source                % seq_match.m
addpath /asl/packages/airs_decon/source           % hamm_app
addpath /asl/matlib/plotutils                     % aslprint.m
addpath /home/chepplew/myLib/matlib/math          % remove_6sigma

s = struct;

% set up the dates to process:
if (length(sdate) ~= 10) fprintf(1,'Error in date\n'); exit; end
syr = sdate(1:4);   smn = sdate(6:7);    sdy = sdate(9:10);
nyr = str2num(syr); nmn = str2num(smn);  ndy = str2num(sdy);

% load the frequency grids:
xx=load('/asl/data/airs/airs_freq.mat'); s.fa=xx.freq; clear xx;
load('/asl/data/iremis/danz/iasi_f.mat'); s.fiasi=fiasi;           % fiasi [8461x1]
xx=load('cris_freq_2grd.mat');  s.fc = xx.vchan'; clear xx;       % 1317 chns (12 guard chans)

% subset between the IASI [8461] and AIRS grids [2378]
[xa xd] = seq_match(sort(s.fa), s.fiasi);

% set up the channels to process (applies to the common grid: 2355 total):
nLW = 1262;  nMW = 602;  nSW = 491;
if(igrp < 1 || igrp > 3) fprintf(1,'igrp out of range (1 to 3)\n'); exit; end
if(igrp == 1) s.achns = [1:1262];    cband = 'LW'; nchns = nLW; end          % applies to the IASI->AIRS
if(igrp == 2) s.achns = [1263:1864]; cband = 'MW'; nchns = nMW; end
if(igrp == 3) s.achns = [1865:2355]; cband = 'SW'; nchns = nSW; end
fprintf(1,'Doing channels %d to %d (freq %8.3f to %8.3f)\n',s.achns(1),s.achns(end),...
   s.fa(s.achns(1)), s.fa(s.achns(end)) );

% for plotting wavenumber scale:
bands = [640, 900; 900, 1320; 1300 2560];                % fcris(ichns(1):ichns(end));
bands = [640, 1100; 1200, 1760; 2150 2550];              % fcris(ichns(1):ichns(end));
bands = [640, 1140; 1210, 1615; 2182 2550];              % 

% Get quantile profiler and set to prf.
%load('/home/strow/Matlab/Sno/prob_vector.mat');  p = p(1:200:end);
%load('/home/chepplew/projects/sno/prob_vector61.mat');               % p [1x61]
% Alternate profiler - softer tails than p.
%junk = 0.0:0.1:10; yp = sigmf(junk,[2 4]); clear junk; 
junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
% Choose which profiler to use (goes in prf)
s.prf  = yp;

% Prep the source directory and files before load.
clear x cc ccs;
dp = '/home/chepplew/data/sno/airs_iasi/JPL/'; % standard/';
unix(['cd ' dp '; find . -noleaf -maxdepth 1 -type f -name ''sno_airs_iasi_*.mat'' -printf ''%P\n'' > /tmp/fn.txt;']);
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

% pre-count the number of SNO pairs available
a = dir([dp 'sno_airs_iasi_*.mat']);
ch = 402; clear fs;
for i=1:1:length(a);
  vars = whos('-file',strcat(dp,a(i).name));
  if( ismember('i2ra', {vars.name}) )                   % AIRS->CrIS is present
    alat = load([dp a(i).name],'alat');
    fs(i) = length(alat.alat);
    fprintf(1,'.')
  end
end
na = sum(fs);
disp(['total SNO pairs to load: ' num2str(na)]);

% Confirm the date range before loading
dstart = datenum([nyr nmn ndy]);
dlast  = datenum([2014 03 01]);
ifn = 1;
for i=1:numel(ccs)
  junk = ccs{i}(15:22);
  thisdat = datenum( [str2num(junk(1:4)) str2num(junk(5:6)) str2num(junk(7:8))] );
  if(thisdat <= dstart) ifn1 = i; end
  if(thisdat <= dlast)  ifn2 = i; end
end
fprintf(1,'Processing SNO files from %s to %s\n',ccs{ifn1}, ccs{ifn2});

% Now load the data
clear g;
s.td    = [];   s.dsn  = [];  s.arad = [];  s.drad = [];
s.atim  = [];   s.alat = [];  s.alon = []; s.asolz = [];
s.itim  = [];   s.ilat = [];  s.ilon = []; s.isolz = [];  s.iqual = [];  
for ifnum = ifn1:ifn2
  vars = whos('-file',strcat(dp,ccs{ifnum}));
    %if(snoLst(ifnum).bytes > 1.0E4)
      if( vars(4).size(1) == vars(15).size(1) )         % ensure same no. obs
      g = load(strcat(dp,ccs{ifnum}));
%      if(ismember('i2ra',{vars.name})) 
       s.arad  = [s.arad, g.ra(s.achns,:)];               % AIRS
       s.drad  = [s.drad, g.i2ra(s.achns,:)];             % IASI->AIRS
%      end
      s.iqual = [s.iqual, g.iqual'];
      s.alat  = [s.alat,  g.alat'];
      s.alon  = [s.alon,  g.alon'];
      s.asolz = [s.asolz, g.asolzen'];
      s.ilat  = [s.ilat,  g.ilat'];
      s.ilon  = [s.ilon,  g.ilon'];
      s.isolz = [s.isolz, g.isolzen'];
      s.td    = [s.td,    g.tdiff'];
      s.dsn   = [s.dsn,   g.dist'];
      end
    %end
  fprintf(1,'.');
end
clear g;
fprintf(1,'\n');
[nx ny] = size(s.arad); [mx my] = size(s.drad);
fprintf(1,'number of AI SNO pairs: %d\t of AI SNO pairs: %d\n', ny, my);

% quality control
inq = find(s.iqual ~= 0);
s.arad(inq) = NaN;  s.drad(inq) = NaN;                 % need to match pairs.
fprintf(1,'Found %d bad IASI radiances\n',numel(inq));

% Remove outliers
disp(['Removing outliers']);
aibias = s.arad - s.drad;
clear gx;
for i=1:length(nchns)
   n  = single(remove_6sigma(aibias(i,:)));
   nn = single(remove_6sigma(aibias(i,n)));
   gx(i).n = n(nn);
end

% Now find unique set of bad SNO samples
ux = [];
[~, psz] = size(aibias);
for i=1:length(nchns)
   ux = [ux setdiff(1:psz,gx(i).n)];
end
ux  = unique(ux);
s.ing = single(setdiff(1:psz,ux));
disp(['  ' num2str(numel(ux)) ' outliers removed']);
clear gx n nn ux psz aibias;

% --------------    subset night or day ---------------------
idn = find(s.asolz);                                                    % include all scenes for LW and MW bands
  s.idn = find(s.asolz > 95 & s.isolz > 95);
  s.idd = find(s.asolz < 90 & s.isolz < 90);
  fprintf(1,'Found %d day, %d night\n',numel(s.idd),numel(s.idn));
% --------------    subset NH or SH ---------------------
idn = find(s.alat);                                                    % include all scenes for LW and MW bands
  s.ish = find(s.alat < 0);
  s.inh = find(s.alat > 0);
  fprintf(1,'Found %d north, %d south\n',numel(s.inh),numel(s.ish));

%{
% Choose which subset to use
clear day nit nor sou;
idx = ':';
if(LDAY) clear idx; idx = intersect(ing, idd);  end
if(LNIT) clear idx; idx = intersect(ing, idn);  end
if(LNOR) clear idx; idx = intersect(ing, inh);  end
if(LSOU) clear idx; idx = intersect(ing, ish);  end
%}

% Convert rad to BT, do global stats (also used in next section)  
% remove bad data first:
clear abt dbt adm dbm radstd btstd stderr;
idx = s.ing;
ny  = numel(idx);
abt = single( real(rad2bt(s.fa(s.achns), s.arad(:,idx))) );
dbt = single( real(rad2bt(s.fa(s.achns), s.drad(:,idx))) );
abm = nanmean(abt,2);
dbm = nanmean(dbt,2);

radstd   = nanstd( s.drad(:,idx) - s.arad(:,idx), 0, 2 );
adbm     = 0.5*( nanmean(dbt,2) + nanmean(abt,2) );
  mdr    = 1E-3*( 1./drdbt(s.fa(s.achns),adbm) );
btstd    = mdr.*radstd;  
stderr   = btstd./sqrt(ny)';
whos idx abt dbt arad drad abm dbm radstd btstd stderr

%{
figure;plot(s.fa(s.achns), abm,'b-',s.fa(s.achns),dbm,'g-');grid on;
figure;plot(s.fa(s.achns), abm - dbm, 'b-', s.fa(s.achns),stderr,'g-');grid on;
  axis([bands(igrp,:) -2 2]);
%}
% -------------------------- quantiles --------------------------- %
disp('working on quantiles...');
 
% create the scene bins for each channel (AIRS / I2A)
clear qaBins qxBins qdBins qx qa qd qsBins;
qxBins.B = quantile(abt,s.prf,2);
qdBins.B = quantile(dbt,s.prf,2);
qx       = cell2mat(struct2cell(qxBins));
qd       = cell2mat(struct2cell(qdBins));
qsBins   = (qx + qd)/2.0;                        % x:AIRS & d:IASI-to-AIRS

% populate the scene bins for each channel (jj)
clear binsz btbias radstd adbm btstd btser bias250;
for jj = 1:numel(s.achns)
  sbins = qsBins(jj,:);
  clear dbin dbinStd dbinN dbinInd abin abinStd abinN abinInd ubinInd;
  [dbin dbinStd dbinN dbinInd] = Math_bin(dbt(jj,:),dbt(jj,:),sbins); 
  [abin abinStd abinN abinInd] = Math_bin(abt(jj,:),abt(jj,:),sbins);

  for i = 1:length(dbin)                                                       
    ubinInd(i,:) = {union(dbinInd{i},abinInd{i})};                  
  end

  for i = 1:length(dbin)
    binsz(jj,i)   = length(ubinInd{i});
    btbias(jj,i)  = nanmean( dbt(jj,ubinInd{i}) - abt(jj,ubinInd{i}) );
    radstd(jj,i)  = nanstd( s.arad(jj,ubinInd{i}) - s.drad(jj,ubinInd{i}) );
    adbm(i)    = 0.5*( nanmean(dbt(jj,ubinInd{i})) + nanmean(abt(jj,ubinInd{i})) );
      mdr      = 1E-3*( 1./drdbt(s.fa(jj),adbm(i)) );
    btstd(jj,i)   = mdr.*radstd(jj,i);  
    sterr(jj,i)   = btstd(jj,i)./sqrt(binsz(jj,i));
    %%bias250(jj,i) = btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  fprintf(1,'.');
end
fprintf(1,'\n');


% parameter fitting section
% -------------------------
wmstats = struct; blo =[]; bhi = [];
wmstats.abm    = abm;
wmstats.dbm    = dbm;
wmstats.wn     = s.fa(s.achns);
wmstats.stderr = stderr;

for jj = 1:numel(s.achns)
  clear junk;
  wmstats.bias{jj}  = btbias(jj,:);
  wmstats.sterr{jj} = sterr(jj,:);
  wmstats.binqa{jj} = qsBins(jj,1:end-1);
  wmstats.binsz{jj} = single(binsz(jj,:));

  % weighted mean & stats section
  % -----------------------------
  % 1. full range
  % -------------
  qlo(jj,1)  = 1;
  qhi(jj,1)  = 200;
  jtot = sum(binsz(jj,:));
  for i = 1:length(dbin)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  wmstats.b(jj,1)  = nansum(junk);  
  wmstats.mx(jj,1) = nanmax(btbias(jj,:));  
  wmstats.mn(jj,1) = nanmin(btbias(jj,:));
  wmstats.lo(jj,1) = qsBins(jj,1);
  wmstats.hi(jj,1) = qsBins(jj,end);
  if( ~isnan(wmstats.mx(jj,1)) )
    wmstats.ph(jj,1) = qsBins(jj, find(btbias(jj,:) == wmstats.mx(jj,1),1) );
    wmstats.pl(jj,1) = qsBins(jj, find(btbias(jj,:) == wmstats.mn(jj,1),1) );
  end
  clear junk;
  % -------------------------------
  % 2. range set by min sample size
  % -------------------------------
  inband = find(binsz(jj,:) > 500);
  if(numel(inband) < 2) fprintf(1,'ichn: %d\t binsz too small\n',jj); continue; end
  blo(jj,1) = min(inband); 
  bhi(jj,1) = max(inband);                                      % was min(max(inband), find(qsBins(jj,:) > 297,1));
  jtot = sum(binsz(jj,blo(jj,1):bhi(jj,1)));
  for i = blo(jj,1):bhi(jj,1)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  wmstats.b(jj,2)  = nansum(junk);  clear junk;
  wmstats.mx(jj,2) = nanmax(btbias(jj,blo(jj,1):bhi(jj,1)));  
  wmstats.mn(jj,2) = nanmin(btbias(jj,blo(jj,1):bhi(jj,1)));
  wmstats.lo(jj,2) = qsBins(jj,blo(jj,1));
  wmstats.hi(jj,2) = qsBins(jj,bhi(jj,1));
  wmstats.ph(jj,2) = qsBins(jj, find(btbias(jj,:) == nanmax(btbias(jj,blo(jj,1):bhi(jj,1))),1) ); 
  wmstats.pl(jj,2) = qsBins(jj, find(btbias(jj,:) == nanmin(btbias(jj,blo(jj,1):bhi(jj,1))),1) );
  % ---------------------------
  % 3. range set by max std.err
  % ---------------------------
  inband = find(sterr(jj,:) < 0.04 & binsz(jj,:) > 500);        % highly tuned by trial n error
  if(numel(inband) < 2) fprintf(1,'ichn: %d\t Std Err too large\n',jj); continue; end
  blo(jj,2) = min(inband);                                      % prob OK for all wns. 
  bhi(jj,2) = max(inband);
  jtot = sum(binsz(jj,blo(jj,2):bhi(jj,2)));
  for i = blo(jj,2):bhi(jj,2)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  wmstats.b(jj,3)  = nansum(junk);  clear junk;
  wmstats.mx(jj,3) = nanmax(btbias(jj,blo(jj,2):bhi(jj,2)));  
  wmstats.mn(jj,3) = nanmin(btbias(jj,blo(jj,2):bhi(jj,2)));
  wmstats.lo(jj,3) = qsBins(jj,blo(jj,2));
  wmstats.hi(jj,3) = qsBins(jj,bhi(jj,2));
  wmstats.ph(jj,3) = qsBins(jj, find(btbias(jj,:) == nanmax(btbias(jj,blo(jj,2):bhi(jj,2))),1) );
  wmstats.pl(jj,3) = qsBins(jj, find(btbias(jj,:) == nanmin(btbias(jj,blo(jj,2):bhi(jj,2))),1) );
    
end

% save file
% ---------
savfn = ['AI_SNO_2013_wmstats_' cband '.mat'];
fprintf(1,'Saving: %s\n',savfn);
save(savfn,'wmstats');

% display summary of results for selected channel
% -----------------------------------------------
ich = find(s.fa>900,1);   %  759;
sfnam = fieldnames(wmstats);
disp([wmstats.wn(ich)]);
for i = 9:numel(sfnam)
  disp([wmstats.(sfnam{i})(ich,:)]);
end

%{
%    ---------- Plotting Section -------------------
addpath /asl/matlib/plotutils              % aslprint.m
phome = '/home/chepplew/projects/sno/sno_paper_2016/figs/';


% ----- Plot Quantile Stats  -------
bincen = qsBins(ich,1:end-1);
figure(2);clf;h1=subplot(2,1,1);
  plot(bincen,wmstats.bias{ich},'.-',bincen,wmstats.sterr{ich},'r-',...
       bincen,-wmstats.sterr{ich},'r-');
  grid on;axis([200 290 -1 1]);ylabel('Bias (K)');
  %title('AIRS IASI SNO bias 900 wn vs Scene');
  h2=subplot(2,1,2);semilogy(bincen,wmstats.binsz{jj});grid on;axis([200 290 100 20000]);
  xlabel('Scene Temperature (K)');ylabel('population');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])  
  % aslprint([phome 'AI_jplSNO_bias_std_900wn_vScene.png']);



% ----------------
bincen = qsBins(ich,1:end-1);
figure(1);clf;
  subplot(2,1,1);plot(bincen,btbias(ich,:),'k.-',bincen,btser(jj,:),'r-');grid on;
    axis([190 300 -1 1]);
  subplot(2,1,2);semilogy(bincen,binsz(jj,:));axis([190 300 50 20000]);grid on;
figure(2);clf;subplot(2,1,1)
  plot(bincen(blo(jj,1):bhi(jj,1)),btbias(jj,blo(jj,1):bhi(jj,1)));grid on;
 subplot(2,1,2);plot(bincen(blo(jj,2):bhi(jj,2)),btbias(jj,blo(jj,2):bhi(jj,2)));grid on;

% -------------- Global Stats ---------------------
figure(3);clf;h1=subplot(2,1,1);plot(wmstats.wn,wmstats.abm,'b-',wmstats.wn,wmstats.dbm,'g-');
  grid on;ylabel('BT K');xlim([bands(igrp,:)])
  title('Airs (g) IASI (b) 2011-14 SNO mean BT');
  h2=subplot(2,1,2);plot(wmstats.wn,wmstats.dbm - wmstats.abm,'-');
  grid on;xlabel('wavenumber');ylabel('BT Bias K');axis([bands(igrp,:) -1 1]);
  %%ha=findobj(gcf,'type','axes');set(ha(1),'ylim',[-Inf Inf]);
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp=get(h1,'position');
  set(h1,'position',[pp(1) pp(2)-pp(4)*0.1 pp(3) pp(4)*1.1])
  pp=get(h2,'position'); set(h2,'position',[pp(1) pp(2) pp(3) pp(4)*1.1]);
  aslprint(['./figs/AI_SNO_2011x_BTSpectrum_chns' sprintf('%d',igrp) '.png']);

%}


%{
there is one structure: wmstats
There are 200 channels. 
wmstats = 
      cbm: [200x1 single]
      dbm: [200x1 single]
       wn: [1x200 double]
     bias: {1x200 cell}
    btser: {1x200 cell}
     bins: {1x200 cell}
        b: [200x3 single]
       mx: [200x3 single]
       mn: [200x3 single]
       lo: [200x3 single]
       hi: [200x3 single]
       ph: [200x3 single]
       pl: [200x3 single]


wn    = wavenumber of the channel (on the common grid). 
cbm   = mean CrIS BT spectrum (for the 200 channels in group)
dbm   = mean AIRS-to-CrIS BT spectrum (ditto)
bias  = the BT bias A-C for the 200 scene bins.
btser = the standard error for the bias for the 200 scene bins.
bins  = the scene bins used for the bias calcs.
b     = 3 x weighted average values per channel:(full range; 
      restricted range by sample size, restricted range by std.err)
mx: 3 x maximum BT bias (K) for repsective channel and range.
mn: 3 x minimum BT bias (ditto)
lo: 3 x actual min scene temperature in range used
hi: 3 x actual max scene temp in range used
ph: 3 x the scene temp at which the max BT bias occurs
pl: 3 x the scene temp at which the min BT bias occurs

%}
