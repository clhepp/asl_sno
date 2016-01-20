function wmstats = allchns_sno_iasi_cris_jpl_mat(sdate, igrp)

% Load and prep stats for all channels from a SMALL SNO set

% INPUT: sdate calendar date to start collection
%        e.g. sdate='2012/04/01';
%        igrp - the group number for the channels in multiples of 200.
%           selected to cover the common grid channels (1185)
%           valid range 1:6.
%
% OUTPUT: wmstats - see notes at bottom for details
%
% VERSION: 1: tuned to accommodate 200 channels in one gulp loading 30 days
%         of SNO data.
%         Replaced dir() with unix(find...) which generates temporary file.
%
% NOTES: allow 10 GB on compute node.
%
% C Hepplewhite. November 2015.

%%clearvars s a wmstats g -except sdate igrp;

cd /home/chepplew/gitLib/asl_sno/run
addpath /home/chepplew/gitLib/asl_sno/source
addpath /home/chepplew/gitLib/asl_sno/data        % cris_freq_*.mat
addpath /home/chepplew/myLib/matlib/aslutil       % rad2bt.m
addpath /home/strow/Git/breno_matlab/Math         % Math_bin.m
addpath /asl/matlab2012/aslutil/                  % drdbt.m
addpath /asl/packages/ccast/source                % seq_match.m
addpath /asl/packages/airs_decon/source           % hamm_app

% set up the dates to process:
if (length(sdate) ~= 10) fprintf(1,'Error in date\n'); exit; end
syr = sdate(1:4);   smn = sdate(6:7);    sdy = sdate(9:10);
nyr = str2num(syr); nmn = str2num(smn);  ndy = str2num(sdy);

% set up the channels to process (applies to the common grid):
if (igrp < 1 || igrp > 3) fprintf(1,'igrp out of range (1 to 3)\n'); exit; end
ichns = [(igrp-1)*400 + 1:igrp*400];                     % applies to the IASI->AIRS
if(igrp == 3) ichns = [801:1317]; end
fprintf(1,'Doing channels %d to %d\n',ichns(1),ichns(end));

% for plotting at the end:
bands = [640, 900; 900, 1320; 1300 2560];                % fcris(ichns(1):ichns(end));

% load the frequency grids:
xx=load('cris_freq_2grd.mat');  fcris = xx.vchan; clear xx;    % 1317 chns (12 guard chans)
load('/asl/data/iremis/danz/iasi_f.mat');                      % fiasi [8461x1]
[xc xi] = seq_match(sort(fcris), fiasi);

% Get quantile profiler and set to prf.
%load('/home/strow/Matlab/Sno/prob_vector.mat');  p = p(1:200:end);
%load('/home/chepplew/projects/sno/prob_vector61.mat');               % p [1x61]
% Alternate profiler - softer tails than p.
%junk = 0.0:0.1:10; yp = sigmf(junk,[2 4]); clear junk; 
junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
% Choose which profiler to use (goes in prf)
prf  = yp;


clear x cc ccs;
dp = '/asl/s1/chepplew/projects/sno/iasi_cris/JPL/'; % standard/';
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

dstart = datenum([nyr nmn ndy]);
for i=1:numel(ccs)
  junk = ccs{i}(14:21);                                  % sno_iasi_cris20120401.mat
  thisdat = datenum( [str2num(junk(1:4)) str2num(junk(5:6)) str2num(junk(7:8))] );
  if(thisdat <= dstart) ifn1 = 1; end
  %%if(thisdat <= dlast)  ifn2 = i; end
end
if(numel(ccs) > 36) ifn2 = ifn1 + numel(ccs);                      % only 36 SNO files from 2001/01/01.
else ifn1 = 1; ifn2 = numel(ccs);
end
%fprintf(1,'Processing SNO files from: %s to %s\n',snoLst(ifn1).name, snoLst(ifn2).name);
fprintf(1,'Processing SNO files from %s to %s\n',ccs{ifn1}, ccs{ifn2});

clear g;
s = struct;
s.td    = [];   s.crad  = [;]; s.irad = [;]; s.drad = [;]; s.itim = [];  s.ctim = []; 
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
        g = load(strcat(dp,ccs{ifnum}));
        s.irad  = [s.irad,  g.ri(xi(ichns),:)];
        s.crad  = [s.crad,  g.rc(ichns,:)];              % 
        s.drad  = [s.drad,  g.i2rc(ichns,:)];            %
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

% quality control
inq = find(s.iqual ~= 0);
fprintf(1,'Found %d bad i2ra radiances\n',numel(inq));
s.drad(inq) = NaN;  s.crad(inq) = NaN;  s.irad(inq) = NaN;   % need to match pairs.

% subset night (for SW band :- igrp=3)
idn = ':';                                                    % include all scenes for LW and MW bands
if(igrp == 3 )
  idn = find(s.csolz > 95 & s.isolz > 95);
  fprintf(1,'Found %d night\n',numel(idn));
end
 
% convert Obs to BT
clear cbt ibt dbt cbm ibm dbm;
  junk = single( hamm_app(double(s.irad(:,idn))) );
ibt    = real(rad2bt(fiasi(xi(ichns))',junk) );
  junk = single( hamm_app(double(s.crad(:,idn))) );
cbt    = real( rad2bt(fcris(ichns),junk) );
  junk = single( hamm_app(double(s.drad(:,idn))) );
dbt   =  real( rad2bt(fcris(ichns), junk) );
crad  = s.crad(:,idn);  
cbm   = nanmean(cbt,2); 
dbm   = nanmean(dbt,2); 
ibm   = nanmean(ibt,2);       
whos cbt dbt ibt crad s cbm dbm ibm

%{ 
% Sanity check
figure(1);clf;plot(fcris(ichns),cbm,'b-',fcris(ichns),dbm,'g-');grid on;axis([bands(igrp,:) 215 255]);
figure(1);clf;plot(fcris(ichns),cbm - dbm, 'm.-');grid on;axis([bands(igrp,:) -0.6 0.6]);
%}

% create the scene bins for each channel
clear qcBins qdBins qc qd qsBins;
qcBins.B = quantile(cbt,prf,2);
qdBins.B = quantile(dbt,prf,2);
qc       = cell2mat(struct2cell(qcBins));
qd       = cell2mat(struct2cell(qdBins));
qsBins   = (qc + qd)/2.0;                        % c:CrIS & d:IASI-to-CrIS

% populate the scene bins for each channel (jj)
clear binsz btbias radstd cdbm btstd btser bias250;
for jj = 1:numel(ichns)
  sbins = qsBins(jj,:);
  clear dbin dbinStd dbinN dbinInd abin abinStd abinN abinInd ubinInd;
  [dbin dbinStd dbinN dbinInd] = Math_bin(dbt(jj,:),dbt(jj,:),sbins); 
  [cbin cbinStd cbinN cbinInd] = Math_bin(cbt(jj,:),cbt(jj,:),sbins);

  for i = 1:length(dbin)                                                       
    ubinInd(i,:) = {union(dbinInd{i},cbinInd{i})};                  
  end

  for i = 1:length(dbin)
    binsz(jj,i)   = length(ubinInd{i});
    btbias(jj,i)  = nanmean( dbt(jj,ubinInd{i}) - cbt(jj,ubinInd{i}) );
    radstd(jj,i)  = nanstd( s.drad(jj,ubinInd{i}) - s.crad(jj,ubinInd{i}) );
    cdbm(i)    = 0.5*( nanmean(dbt(jj,ubinInd{i})) + nanmean(cbt(jj,ubinInd{i})) );
      mdr      = 1E-3*( 1./drdbt(fcris(jj),cdbm(i)) );
    btstd(jj,i)   = mdr.*radstd(jj,i);  
    btser(jj,i)   = btstd(jj,i)./sqrt(binsz(jj,i));
    %%bias250(jj,i) = btbias(jj,i)./drd250(jj);                 % option hard wired
  end
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
wmstats.wn  = fcris(ichns);

for jj = 1:numel(ichns)
  clear junk;
  wmstats.bias{jj}  = btbias(jj,:);
  wmstats.btser{jj} = btser(jj,:);
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

% Plotting section
% ----------------
addpath /asl/matlib/plotutils              % aslprint.m
jj = 200;
bincen = wmstats.bins{jj};                 % qsBins(jj,1:end-1);
figure(1);clf;
  h1=subplot(2,1,1);plot(bincen,btbias(jj,:),'k.-',bincen,btser(jj,:),'r-');grid on;
    axis([150 280 -20 1]);title('IC SNO Night Scene Dep. ch:2487wn');
  h2=subplot(2,1,2);semilogy(bincen,binsz(jj,:));axis([150 280 50 20000]);grid on;xlabel('Scene K');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp=get(h1,'position');
  set(h1,'position',[pp(1) pp(2)-pp(4)*0.1 pp(3) pp(4)*1.1])
  pp=get(h2,'position'); set(h2,'position',[pp(1) pp(2) pp(3) pp(4)*1.1]);
figure(2);clf;subplot(2,1,1)
  plot(bincen(blo(jj,1):bhi(jj,1)),btbias(jj,blo(jj,1):bhi(jj,1)));grid on;
 subplot(2,1,2);plot(bincen(blo(jj,2):bhi(jj,2)),btbias(jj,blo(jj,2):bhi(jj,2)));grid on;

figure(3);clf;h1=subplot(2,1,1);plot(wmstats.wn,wmstats.cbm,'b-',wmstats.wn,wmstats.dbm,'g-');
  grid on;title('CrIS (b) IASI (g) 2012-14 SNO mean BT');ylabel('BT K');xlim(bands(igrp,:));
  h2=subplot(2,1,2);plot(wmstats.wn,wmstats.cbm - wmstats.dbm,'m.-');
  grid on;xlabel('wavenumber');ylabel('BT Bias K');axis([bands(igrp,:) -0.7 0.7]);
  legend('CrIS - I2C','Location','North');
  %%ha=findobj(gcf,'type','axes');set(ha(1),'ylim',[-Inf Inf]);
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp1=get(h1,'position');
  set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position'); set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1]);
  %aslprint(['./figs/IC_SNO_2012x_I2C_C_BTSpectrum_chns' sprintf('%d',igrp) '.png']);

%{
there is one structure: wmstats
There are 400 channels. 
wmstats = 
      cbm: [400x1 single]
      dbm: [400x1 single]
       wn: [1x400 double]
     bias: {1x400 cell}
    btser: {1x400 cell}
     bins: {1x400 cell}
        b: [400x3 single]
       mx: [400x3 single]
       mn: [400x3 single]
       lo: [400x3 single]
       hi: [400x3 single]
       ph: [400x3 single]
       pl: [400x3 single]


wn     = wavenumber of the channel (on the common grid). 
cbm    = mean CrIS BT spectrum (for the 400 channels in group)
dbm    = mean IASI-to-CrIS BT spectrum (ditto)
bias   = the BT bias I-C for the 400 scene bins.
btser  = the standard error for the bias for the 400 scene bins.
binqa  = the scene bins used for the bias calcs from the quantile profiler (in K).
binsz  = the population in each bin of the quantile profiler for each channel.
b      = 3 x weighted average values per channel:(full range; 
      restricted range by sample size, restricted range by std.err)
mx: 3 x maximum BT bias (K) for repsective channel and range.
mn: 3 x minimum BT bias (ditto)
lo: 3 x actual min scene temperature in range used
hi: 3 x actual max scene temp in range used
ph: 3 x the scene temp at which the max BT bias occurs
pl: 3 x the scene temp at which the min BT bias occurs

%}
