function [s wmstats] = allchns_sno_iasi_cris_jpl_mat(sdate, igrp)

% Load and prep stats for all channels from a SMALL SNO set

% INPUT: sdate calendar date to start collection
%        e.g. sdate='2012/04/01';
%        igrp - the group number for the channels in multiples of 200.
%           selected to cover the common grid channels (1185)
%           valid range 1:6.
%
% OUTPUT: (intermediate)
%        s.irad  IASI radiances for selected channels
%        s.crad  CrIS radiances for selected channels  
%        s.drad  IASI->CrIS radiances for selected channels 
%        s.csolz Cris Solar Zenith angle
%	 s.isolz IASI Solar Zenith Angle
%        s.iqual IASI quality flag
%	 s.td    SNO time difference between each sensor obs.
%	 s.dsn   SNO distance between each sensor obs.
%	 s.clat  CrIS Obs Latitude      s.clon CrIS Obs Longitude
%	 s.ilat  IASI Obs Latitude      s.ilon IASI Obs Longitude
%        s.ing   good channels retained (after removing 6-sigma outliers)
%
%        wmstats: structure with fields of data:
%        wmstats.cbm            CrIS mean spectral BT of the global sample.
%        wmstats.dbm            A2C  mean spectral BT of the global sample.
%        wmstats.wn             wavenumbers in group. [nv size]
%        wmstats.stderr         standard error of the mean bias of the global sample.
%        wmstats.bias{jj}       bias in quantile bin [nq x nv cells]
%        wmstats.btser{jj}      std error of mean bias in quantile bin [nq x nv cells]
%        wmstats.binqa{jj}      Quantile BT bins used [nq x nv cells] 
%        wmstats.binsz{jj}      Sample population in each quantile bin [nq x nv cells]
%        wmstats.b:             weighted mean bias from all quantile bins [713 x 3 single]
%        wmstats.mx:            Maximum BT bias found in each quantile bin [713 x 3 single]
%        wmstats.mn:            Minimum BT bias found in each quantile bin [713 x 3 single]
%        wmstats.ph:            the scene temp at which the max BT bias occurs [713 x 3 single]
%        wmstats.pl:            the scene temp at which the min BT bias occurs[713 x 3 single]
%        wmstats.blo:           low-end BT of quantile bin corresp. to qlo [713 x 3 single]
%        wmstats.bhi:           high-end BT of quantile bin corresp. to qhi.[713 x 3 single]
%        wmstats.qlo:           low-end quantile bin number [713 x 3 double]
%        wmstats.qhi:           high-end quantile bin number [713 x 3 double]
%
% VERSION: 2: gulp a group of channels in one go loading up to 30 mat files
%         of SNO data.
%         Replaced dir() with unix(find...) which generates temporary file.
%
% NOTES: allow 10 GB on compute node.
%
% C Hepplewhite. November 2015.
% CLH. Mar2017. Updated path to sno data.
%               changed grouping to match LW MW and SW sub-bands

%%clearvars s a wmstats g -except sdate igrp;

cd /home/chepplew/gitLib/asl_sno/run
addpath /home/chepplew/gitLib/asl_sno/source
addpath /home/chepplew/gitLib/asl_sno/data        % cris_freq_*.mat
addpath /home/chepplew/myLib/matlib/aslutil       % rad2bt.m
addpath /home/strow/Git/breno_matlab/Math         % Math_bin.m
addpath /asl/matlab2012/aslutil/                  % drdbt.m
addpath /asl/packages/ccast/source                % seq_match.m
addpath /asl/packages/airs_decon/source           % hamm_app
addpath /home/chepplew/myLib/matlib/math          % remove_6sigma

s = struct;

LDAY=false;
LNIT=false;
LNOR=false;
LSOU=false;

% set up the dates to process:
if (length(sdate) ~= 10) fprintf(1,'Error in date\n'); exit; end
syr = sdate(1:4);   smn = sdate(6:7);    sdy = sdate(9:10);
nyr = str2num(syr); nmn = str2num(smn);  ndy = str2num(sdy);

% set up the channels to process (applies to the 1317 channel common grid w/2 guard chans):
%%nLW = 713; nMW = 441; nSW = 163;
nLW = 717; nMW = 437; nSW = 163;
if(igrp < 1 || igrp > 3) fprintf(1,'igrp out of range (1 to 3)\n'); exit; end
if(igrp == 1) s.ichns = [1:nLW];          cband = 'LW'; nchns = nLW;  end
if(igrp == 2) s.ichns = [nLW+1: nLW+nMW]; cband = 'MW'; nchns = nMW;  end
if(igrp == 3) s.ichns = [nLW+nMW+1:1185]; cband = 'SW'; nchns = nSW;  end
fprintf(1,'Doing channels %d to %d\n',s.ichns(1),s.ichns(end));

% for plotting at the end:
bands = [640, 900; 900, 1320; 1300 2560];                % fcris(ichns(1):ichns(end));
bands = [640, 1095; 1210, 1760; 2150 2550];                % fcris(ichns(1):ichns(end));

% load the frequency grids:
xx=load('cris_freq_2grd.mat');  s.fc = xx.vchan; clear xx;     % 1317 chns (12 guard chans)
load('/asl/data/iremis/danz/iasi_f.mat');  s.fi = fiasi;       % fiasi [8461x1]
[xc xi] = seq_match(sort(s.fc), s.fi);
%s.fd = f_cris([1:1037 1158:1305]);                            % 1185 prior knowledge from Howards decon routine
%[xd, xc] = seq_match(s.fd, s.fc);                             % track indexes for both grids

% Get quantile profiler and set to prf.
%load('/home/strow/Matlab/Sno/prob_vector.mat');  p = p(1:200:end);
%load('/home/chepplew/projects/sno/prob_vector61.mat');               % p [1x61]
% Alternate profiler - softer tails than p.
%junk = 0.0:0.1:10; yp = sigmf(junk,[2 4]); clear junk; 
junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
% Choose which profiler to use (goes in prf)
s.prf  = yp;


clear x cc ccs;
%dp = '/asl/s1/chepplew/projects/sno/iasi_cris/JPL/'; % standard/';
dp = '/home/chepplew/data/sno/iasi_cris/JPL/';     % standard
unix(['cd ' dp '; find . -noleaf -maxdepth 1 -type f -name ''sno_iasi_cris*.mat'' -printf ''%P\n'' > /tmp/fn.txt;']);
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
  %snoLst = dir(strcat(dp,'sno_iasi_cris_*.mat'));
fprintf(1,'Found %d total SNO files\n',numel(ccs));

a = dir([dp 'sno_iasi_cris201*.mat']);
ch = 402; clear fs;
for i=1:length(a);
  vars = whos('-file',strcat(dp,a(i).name));
  if( ismember('i2rc', {vars.name}) )                   % AIRS->CrIS is present
    clat = load([dp a(i).name],'clat');
    fs(i) = length(clat.clat);
    fprintf(1,'.')
  end
end
na = sum(fs);
disp(['total SNO pairs to load: ' num2str(na)]);

dstart = datenum([nyr nmn ndy]);
for i=1:numel(ccs)
  junk = ccs{i}(14:21);                                  % sno_iasi_cris20120401.mat
  thisdat = datenum( [str2num(junk(1:4)) str2num(junk(5:6)) str2num(junk(7:8))] );
  if(thisdat <= dstart) ifn1 = 1; end
  %%if(thisdat <= dlast)  ifn2 = i; end
end
if(numel(ccs) > 36) 
  ifn2 = ifn1 + numel(ccs);                      % only 36 SNO files from 2001/01/01.
  else ifn1 = 1; ifn2 = numel(ccs);
end
fprintf(1,'Processing SNO files from %s to %s\n',ccs{ifn1}, ccs{ifn2});

clear g;
s.td    = [];   s.crad  = [];  s.irad = [];  s.drad = [];  s.itim = []; s.ctim = []; 
s.clat  = [];   s.clon  = [];  s.dsn  = [];  s.ilat = [];  s.ilon = []; s.csolz = [];
s.iqual = [];   s.isolz = [];
%s.alnfr = [];  s.clnfr = []; s.cifv  = [];
%a.nSam  = [];   a.avrd = [;]; a.avra = [;]; a.avrc = [;]; a.sdra = [;]; a.sdrc = [;]; 
%a.sdrd  = [;];

for ifnum = ifn1:ifn2
  vars = whos('-file',strcat(dp,ccs{ifnum}));
  if( ismember('i2rc', {vars.name}) )                 % IASI->CrIS is present
    if(vars(5).size(1) > 10)
      if( vars(5).size(1) == vars(17).size(1) )       % ensure same no. obs
        fprintf(1,'proc: %s\n', ccs{ifnum})
        g  = load(strcat(dp,ccs{ifnum}));
        rc = single( hamm_app(double(g.rc(s.ichns,:))) );
	rd = single( hamm_app(double(g.i2rc(s.ichns,:))) );
        s.crad  = [s.crad,  rc];              % 
        s.drad  = [s.drad,  rd];            %
        s.irad  = [s.irad,  g.ri(xi(s.ichns),:)];
        s.csolz = [s.csolz, g.csolzen'];
	s.isolz = [s.isolz, g.isolzen'];
        s.iqual = [s.iqual, g.iqual'];
	s.td    = [s.td,    g.tdiff'];
	s.dsn   = [s.dsn,   g.dist'];
	s.clat  = [s.clat,  g.clat'];      s.clon = [s.clon,  g.clon'];
	s.ilat  = [s.ilat,  g.ilat'];      s.ilon = [s.ilon,  g.ilon'];
      else
        fprintf(1,'unequal skip: %d\n',ifnum);
      end
    else
      fprintf(1,'small skip: %d\n', ifnum);
    end
  else
    fprintf(1,'no i2rc: %d\n',ifnum);
  end
  clear vars;
  %%fprintf(1,'.');
end
clear g;
fprintf(1,'\n');
[nx ny] = size(s.drad);
fprintf(1,'number of SNO pairs: %d\n', ny);

% quality control 1: IASI qual flag
inq = find(s.iqual ~= 0);
fprintf(1,'Found %d bad i2ra radiances\n',numel(inq));
s.drad(:,inq) = NaN;  s.crad(:,inq) = NaN;  s.irad(:,inq) = NaN;   % need to match pairs.

% Quality control 2: 6-sigma outliers
disp(['Removing outliers']);
cdbias = s.crad - s.drad;
clear gx;
for i=1:length(s.ichns)
   n  = remove_6sigma(cdbias(i,:));
   nn = remove_6sigma(cdbias(i,n));
   gx(i).n = n(nn);
end

% Now find unique set of bad SNO samples
ux = [];
[~, psz] = size(cdbias);
for i=1:length(s.ichns)
   ux = [ux setdiff(1:psz,gx(i).n)];
end
ux  = unique(ux);
s.ing = single(setdiff(1:psz,ux));
disp(['  ' num2str(numel(ux)) ' outliers removed']);
clear gx n nn ux psz;

% --------------    subset night or day ---------------------
idn = find(s.csolz);                                                    % include all scenes for LW and MW bands
  s.idn = find(s.csolz > 95 & s.isolz > 95);
  s.idd = find(s.csolz < 90 & s.isolz < 90);
  fprintf(1,'Found %d day, %d night\n',numel(s.idd),numel(s.idn));
% --------------    subset NH or SH ---------------------
ihh = find(s.clat);                                                    % include all scenes for LW and MW bands
  s.ish = find(s.clat < 0);
  s.inh = find(s.clat > 0);
  fprintf(1,'Found %d north, %d south\n',numel(s.inh),numel(s.ish));
% Choose which subset to use
clear day nit nor sou;
idx = s.ing;
if(LDAY) clear idx; idx = intersect(s.ing, s.idd);  end
if(LNIT) clear idx; idx = intersect(s.ing, s.idn);  end
if(LNOR) clear idx; idx = intersect(s.ing, s.inh);  end
if(LSOU) clear idx; idx = intersect(s.ing, s.ish);  end

 
% convert rad to BT and get global stats:
[nx ny] = size(s.drad(:,idx));
clear cbt ibt dbt cbm ibm dbm btstd btser;
ibt     = real( rad2bt(fiasi(xi(s.ichns))',s.irad(:,idx)) );
cbt     = real( rad2bt(s.fc(s.ichns), s.crad(:,idx)) );
dbt     = real( rad2bt(s.fc(s.ichns), s.drad(:,idx)) );
cbm     = nanmean(cbt,2); 
dbm     = nanmean(dbt,2); 
ibm     = nanmean(ibt,2);

radstd  = nanstd( s.drad(:,idx) - s.crad(:,idx),0,2 );
cdbm    = 0.5*( nanmean(dbt,2) + nanmean(cbt,2) );
  mdr   = 1E-3*( 1./drdbt(s.fc(s.ichns),cdbm') );
btstd   = mdr.*radstd';  
stderr  = btstd'./sqrt(ny);
whos cbt dbt ibt crad s cbm dbm ibm btstd stderr

if(LDAY) day.cbm = cbm; day.dbm = dbm; day.btser = btser; day.btbias = btbias; 
  day.cbt = cbt; day.dbt = dbt;  end
if(LNIT) nit.cbm = cbm; nit.dbm = dbm; nit.btser = btser; nit.btbias = btbias; 
  nit.cbt = cbt; nit.dbt = dbt;  end
if(LNOR) nor.cbm = cbm; nor.dbm = dbm; nor.btser = btser; nor.btbias = btbias; 
  nor.cbt = cbt; nor.dbt = dbt; end
if(LSOU) sou.cbm = cbm; sou.dbm = dbm; sou.btser = btser; sou.btbias = btbias; 
  sou.cbt = cbt; sou.dbt = dbt; end


% -------------------------- quantiles --------------------------- %
disp('working on quantiles...');
% create the scene bins for each channel
clear qcBins qdBins qc qd qsBins;
qcBins.B = quantile(cbt,s.prf,2);
qdBins.B = quantile(dbt,s.prf,2);
qc       = cell2mat(struct2cell(qcBins));
qd       = cell2mat(struct2cell(qdBins));
qsBins   = (qc + qd)/2.0;                        % c:CrIS & d:IASI-to-CrIS

% populate the scene bins for each channel (jj)
clear binsz btbias radstd cdbm btstd btser bias250;
for jj = 1:numel(s.ichns)
  sbins = qsBins(jj,:);
  clear dbin dbinStd dbinN dbinInd xbin xbinStd xbinN xbinInd ubinInd;
  [dbin dbinStd dbinN dbinInd] = Math_bin(dbt(jj,:),dbt(jj,:),sbins); 
  [cbin cbinStd cbinN cbinInd] = Math_bin(cbt(jj,:),cbt(jj,:),sbins);

  for i = 1:length(dbin)                                                       
    ubinInd(i,:) = {union(dbinInd{i},cbinInd{i})};                  
  end

  for i = 1:length(dbin)
    binsz(jj,i)   = length(ubinInd{i});
    btbias(jj,i)  = nanmean( dbt(jj,ubinInd{i}) - cbt(jj,ubinInd{i}) );
    radstd(jj,i)  = nanstd( s.drad(jj,idx(ubinInd{i})) - s.crad(jj,idx(ubinInd{i})) );
    cdbm(i)    = 0.5*( nanmean(dbt(jj,ubinInd{i})) + nanmean(cbt(jj,ubinInd{i})) );
      mdr      = 1E-3*( 1./drdbt(s.fc(jj),cdbm(i)) );
    btstd(jj,i)   = mdr.*radstd(jj,i);  
    btser(jj,i)   = btstd(jj,i)./sqrt(binsz(jj,i));
    %%bias250(jj,i) = btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  jtot  = sum(binsz(jj,:));
  jmdr  = 1E-3*( 1./drdbt(s.fc(jj),cbm(jj)) );
  jbtse = jmdr.* nanstd(s.drad(jj,idx) - s.crad(jj,idx),1,2) / sqrt(jtot);
  fprintf(1,'.');
end
fprintf(1,'\n');

%{
% Sanity check
jj=400;
figure(2);clf;plot(qsBins(jj,1:end-1),btbias(jj,:),'b.-');grid on;
%}

% parameter fitting section
% -------------------------
wmstats = struct; blo =[]; bhi = [];
wmstats.cbm = cbm;
wmstats.dbm = dbm;
wmstats.wn  = s.fc(s.ichns);
wmstats.stderr = stderr;

for jj = 1:numel(s.ichns)
  clear junk;
  wmstats.bias{jj}   = btbias(jj,:);
  wmstats.btser{jj}  = btser(jj,:);
  wmstats.binqa{jj}  = qsBins(jj,1:end-1);
  wmstats.binsz{jj}  = binsz(jj,:);

  % range select by bin size and hot scenes
  inband = find(binsz(jj,:) > 500);
  if(numel(inband) < 2) fprintf(1,'ichn: %d\t binsz too small\n',jj); continue; end
  blo(jj,1) = min(inband); 
  bhi(jj,1) = max(inband);                                      % was min(max(inband), find(qsBins(jj,:) > 297,1));
  % range select by std err - deal with non-monotonic var at tails
  inband = find(btser(jj,:) < 0.06 & binsz(jj,:) > 500);        % highly tuned by trial n error
  if(numel(inband) < 2) fprintf(1,'ichn: %d\t Std Err too large\n',jj); continue; end
  blo(jj,2) = min(inband);                                      % prob OK for all wns. 
  bhi(jj,2) = max(inband);
  % weighted mean & stats section
  % -----------------------------
  % full range
  % ----------
  jtot = sum(binsz(jj,:));
  for i = 1:length(dbin)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  wmstats.b(jj,1)  = nansum(junk);  
  wmstats.mx(jj,1) = nanmax(btbias(jj,:));  
  wmstats.mn(jj,1) = nanmin(btbias(jj,:));
  wmstats.lo(jj,1) = qsBins(jj,1);
  wmstats.hi(jj,1) = qsBins(jj,end);
  wmstats.ph(jj,1) = qsBins(jj, find(btbias(jj,:) == nanmax(btbias(jj,:)),1) );
  wmstats.pl(jj,1) = qsBins(jj, find(btbias(jj,:) == nanmin(btbias(jj,:)),1) );
  clear junk;
  % ---------------
  % range set by min sample size
  % ---------------
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
  % --------------
  % range set by max std.err
  % ---------------
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
savfn = ['sno_IC_wmstats_2012x_chns400x' sprintf('%02d',igrp) '.mat'];
%{
fprintf(1,'Saving: %s\n',savfn);
save(savfn,'wmstats');

% display summary of results for selected channel
% -----------------------------------------------
jj = 200;
sfnam = fieldnames(wmstats);
disp([wmstats.wn(jj)]);
for i = 7:numel(sfnam)
  disp([wmstats.(sfnam{i})(jj,:)]);
end
%}
%{ 
% ------------------ Plot Options ------------------------------
j = igrp;    ich = 402;
figure(j);clf;h1=subplot(2,1,1);plot(s.fc(s.ichns),day.cbm,'-',fcris(ichns),day.dbm,'-',...
  fcris(ichns),nit.cbm,'-',fcris(ichns),nit.dbm,'-');ylabel('BT (K)');
  grid on;axis([bands(igrp,:) 205 260]);title('IC SNO Day,Night SH Mean, Bias, std error');
  h2=subplot(2,1,2);plot(fcris(ichns),day.btbias,'-',fcris(ichns), day.btser,'-',...
  fcris(ichns),nit.btbias,'-',fcris(ichns), nit.btser,'-');xlabel('wavenumber cm^{-1}');
  grid on;axis([bands(igrp,:),-0.1 0.3]);ylabel('dBT (K)');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp1=get(h1,'position');
  set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position'); set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1]);
figure(j);clf;plot(fcris(ichns),cbm - dbm, '-');grid on;axis([bands(igrp,:) -0.6 0.6]);
figure(j);clf;[hax,hl1,hl2]=plotyy(fcris(ichns),cbm - dbm,fcris(ichns), btser);
  grid on;xlim(hax(:),[bands(igrp,:)]);ylim(hax(1),[-0.2 0.6]);ylim(hax(2),[0 0.05]);
  hax(1).YTick = [-0.2:0.1:0.5];
  xlabel('wavenumber cm^{-1}');ylabel(hax(1),'Mean Bias K');ya2=ylabel(hax(2),'Std. error K');
  set(ya2, 'Units', 'Normalized', 'Position', [1.05, 0.7, 0]);
  title('2012-14 IASI CrIS SNO MW Bias'); legend('bias','std. error','Location','north');
  % aslprint('../figs/IC_jplSNO_Bias_stdErr_MW.png',1);
% for 900 wn channel:
  btbins   = [180:1:300];
  c404_pdf = histcounts(cbt(ich,:),btbins);    d404_pdf = histcounts(dbt(ich,:),btbins);
  btbinx   = (btbins(1:end-1)+btbins(2:end))/2;
  figure(3);clf;hold on;plot(btbinx,c404_pdf,'+-');plot(btbinx,d404_pdf,'o-');grid on;
  aslprint('../figs/IC_jplSNO_900wn_population_vs_scene.pdf',1)
%  day/night
  c404d_pdf = histcounts(day.cbt(ich,:),btbins); d404d_pdf = histcounts(day.dbt(ich,:),btbins);
  c404n_pdf = histcounts(nit.cbt(ich,:),btbins); d404n_pdf = histcounts(nit.dbt(ich,:),btbins);
  figure(2);clf;h1=subplot(2,1,1);hold on;plot(btbinx,c404d_pdf,'+-');plot(btbinx,d404d_pdf,'o-');
    grid on;legend('CrIS day','IASI day');title('IC SNO 900wn sample pop vs scene (day night)');
    h2=subplot(2,1,2);hold on;plot(btbinx,c404n_pdf,'+-');plot(btbinx,d404n_pdf,'o-');grid on;
    xlabel('Scene BT (K)');ylabel('population');legend('CrIS night','IASI night');
    % aslprint('../figs/IC_jplSNO_900wn_population_vs_scene_day_night.png');
% north (boreal)/south (austral) hemisphere
  c404b_pdf = histcounts(nor.cbt(ich,:),btbins); d404b_pdf = histcounts(nor.dbt(ich,:),btbins);
  c404a_pdf = histcounts(sou.cbt(ich,:),btbins); d404a_pdf = histcounts(sou.dbt(ich,:),btbins);
  figure(4);clf;h1=subplot(2,1,1);hold on;plot(btbinx,c404b_pdf,'+-');plot(btbinx,d404b_pdf,'o-');
    grid on;legend('CrIS NH','IASI NH');title('IC SNO 900wn sample pop vs scene (north south)');
    h2=subplot(2,1,2);hold on;plot(btbinx,c404a_pdf,'+-');plot(btbinx,d404a_pdf,'o-');grid on;
    xlabel('Scene BT (K)');ylabel('population');legend('CrIS SH','IASI SH');
    % aslprint('../figs/IC_jplSNO_900wn_population_vs_scene_NH_SH.png');

% ----------------
addpath /asl/matlib/plotutils              % aslprint.m
jj = 200;
bincen = wmstats.binqa{jj};                 % qsBins(jj,1:end-1);
figure(1);clf;
  h1=subplot(2,1,1);plot(bincen,btbias(jj,:),'k.-',bincen,btser(jj,:),'r-');grid on;
    axis([180 300 -2 1]);title('IC SNO Night Scene Dep. ch:2487wn');
  h2=subplot(2,1,2);semilogy(bincen,wmstats.binsz{jj});axis([180 300 50 20000]);grid on;xlabel('Scene K');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp=get(h1,'position');
  set(h1,'position',[pp(1) pp(2)-pp(4)*0.1 pp(3) pp(4)*1.1])
  pp=get(h2,'position'); set(h2,'position',[pp(1) pp(2) pp(3) pp(4)*1.1]);
figure(2);clf;subplot(2,1,1)
  plot(bincen(blo(jj,1):bhi(jj,1)),btbias(jj,blo(jj,1):bhi(jj,1)));grid on;
 subplot(2,1,2);plot(bincen(blo(jj,2):bhi(jj,2)),btbias(jj,blo(jj,2):bhi(jj,2)));grid on;

% Global sample stats:
figure(3);clf;h1=subplot(2,1,1);plot(wmstats.wn,wmstats.cbm,'b-',wmstats.wn,wmstats.dbm,'g-');
  grid on;title('CrIS (b) IASI (g) 2012-14 SNO mean BT');ylabel('BT K');xlim(bands(igrp,:));
  h2=subplot(2,1,2);plot(wmstats.wn,wmstats.cbm - wmstats.dbm,'b-',...
                         wmstats.wn,10.*wmstats.stderr,'g-');
  grid on;xlabel('wavenumber');ylabel('BT Bias K');axis([bands(igrp,:) -0.1 0.5]);
  legend('CrIS - I2C','Location','North');
  %%ha=findobj(gcf,'type','axes');set(ha(1),'ylim',[-Inf Inf]);
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp1=get(h1,'position');
  set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position'); set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1]);
  %aslprint(['./figs/IC_SNO_2012x_I2C_C_BTSpectrum_chns' sprintf('%d',igrp) '.png']);

%}
