function [s wmstats] = allchns_sno_airs_cris_jpl_mat(sdate, igrp)

% Load and prep stats for all channels from a 'SMALL' SNO set (see notes below)
%
% INPUT: sdate calendar date to start collection
%        e.g. sdate='2013/01/02';
%        igrp - the group number for the channels corresponding to CrIS bands.
%           selected to cover the common grid channels (1185)
%           valid range 1:3.
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
%
% VERSION: 1: Gulps subset channels loading 2013 of SNO data.
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
addpath /home/strow/Git/breno_matlab/Math         % Math_bin.m
addpath /asl/matlab2012/aslutil/                  % drdbt.m
addpath /asl/packages/ccast/source                % seq_match.m
addpath /asl/packages/airs_decon/source           % hamm_app
addpath /home/chepplew/myLib/matlib/math          % remove_6sigma

s = struct;

% set up the dates to process:
if (length(sdate) ~= 10) fprintf(1,'Error in date\n'); exit; end
syr = sdate(1:4);   smn = sdate(6:7);    sdy = sdate(9:10);
nyr = str2num(syr); nmn = str2num(smn);  ndy = str2num(sdy);

% load the frequency grids:
xx=load('/asl/data/airs/airs_freq.mat'); s.fa=xx.freq; clear xx;
xx=load('cris_freq_nogrd.mat'); f_cris=xx.vchan; clear xx;    % 1305 chns (no guard chans)
xx=load('cris_freq_2grd.mat');  s.fc = xx.vchan; clear xx;    % 1317 chns (12 guard chans)
s.fd = f_cris([1:1037 1158:1305]);                            % 1185 prior knowledge from Howards decon routine
[xd, xc]    = seq_match(s.fd, s.fc);                          % track indexes for both grids
[ax, ac]    = seq_match(sort(s.fa), s.fc);


% set up the channels to process (applies to the 1185 channel common grid):
%%nLW = 713; nMW = 441; nSW = 163;
nLW = 713; nMW = 324; nSW = 148;
if(igrp < 1 || igrp > 3) fprintf(1,'igrp out of range (1 to 3)\n'); exit; end
if(igrp == 1) s.ichns = [1:nLW];          cband = 'LW'; nchns = nLW;  
  s.achns = [find(s.fa >= s.fc(s.ichns(1)),1):find(s.fa >= s.fc(s.ichns(end)),2)]; end
if(igrp == 2) s.ichns = [nLW+1: nLW+nMW]; cband = 'MW'; nchns = nMW;  
  [s.achns,~] = seq_match(sort(s.fa), s.fc); end
if(igrp == 3) s.ichns = [nLW+nMW+1:1185]; cband = 'SW'; nchns = nSW;  
  [s.achns,~] = seq_match(sort(s.fa), s.fc); end
fprintf(1,'Doing CrIS channels %d to %d\n',s.ichns(1),s.ichns(end));

% for plotting at the end:
bands = [640, 900; 900, 1320; 1300 2560];                  % fcris(ichns(1):ichns(end));
bands = [640, 1100; 1200, 1760; 2150 2550];                % fcris(ichns(1):ichns(end));
bands = [650, 1095; 1210, 1615; 2182 2550];

% Get quantile profiler and set to prf.
%load('/home/strow/Matlab/Sno/prob_vector.mat');  yp = p(1:200:end);
load('/home/chepplew/projects/sno/prob_vector61.mat');               % p [1x61]
% Alternate profiler - softer tails than p.
junk = 0.0:0.1:10; yp = sigmf(junk,[2 4]); clear junk; 
junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
% Choose which profiler to use (goes in prf)
s.prf  = yp;

% ------------- Get source data - start with preparation. ------------- %
dp = '/home/chepplew/data/sno/airs_cris/JPL/standard/';
%dp = '/umbc/lustre/strow/asl/sno/airs_cris/JPL/standard/';
a = dir([dp 'sno_airs_cris_2013*.mat']);
ch = 402; clear fs;
for i=1:1:length(a);
  vars = whos('-file',strcat(dp,a(i).name));
  if( ismember('raDecon', {vars.name}) )                   % AIRS->CrIS is present
    alat = load([dp a(i).name],'alat');
    fs(i) = length(alat.alat);
    fprintf(1,'.')
  end
end
na = sum(fs);
disp(['total SNO pairs to load: ' num2str(na)]);

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

dstart = datenum([nyr nmn ndy]);
dlast  = datenum([2013 12 31]);
ifn1   = 1;
for i=1:numel(ccs)
  junk = ccs{i}(15:22);
  thisdat = datenum( [str2num(junk(1:4)) str2num(junk(5:6)) str2num(junk(7:8))] );
  if(thisdat <= dstart) ifn1 = i; end
  if(thisdat <= dlast)  ifn2 = i; end
end
s.fnames = ccs;
fprintf(1,'Processing SNO files from %s to %s\n',ccs{ifn1}, ccs{ifn2});

% ------------- load data into memory --------------------- %
clear g;
n1 = 1;
s.arad  = single(zeros(length(s.achns),na)); 
s.crad  = single(zeros(nchns,na)); 
s.drad  = single(zeros(nchns,na)); 
s.ctim  = [na]; s.atim = [na]; s.cifov = single(zeros(1,na));
s.alat  = single(zeros(1,na)); s.alon  = single(zeros(1,na));  
s.clat  = single(zeros(1,na)); s.clon  = single(zeros(1,na)); 
s.csolz = single(zeros(1,na)); s.asolz = single(zeros(1,na)); 
s.cland = single(zeros(1,na)); s.aland = single(zeros(1,na));
s.td    = single(zeros(1,na)); s.dsn   = single(zeros(1,na));
for ifnum = ifn1:1:ifn2
  n2 = n1 + fs(ifnum)-1;
  vars = whos('-file',strcat(dp,ccs{ifnum}));
  if( ismember('raDecon', {vars.name}) )                   % AIRS->CrIS is present
      g = load(strcat(dp,ccs{ifnum}));
      if( size(g.ra,2) == size(g.raDecon,2) )               % can get incomplete SNO pairs
        s.arad(:,n1:n2) = g.ra(s.achns,:);                 % [arad, [ra(achn,:); avaw]]; etc
        rc = single(hamm_app(double(g.rc(xc(s.ichns),:))));
	rd = single(hamm_app(double(g.raDecon(xd(s.ichns),:))));
	s.crad(:,n1:n2) = rc;
	s.drad(:,n1:n2) = rd;
%        s.crad  = [s.crad,  g.rc(xc(s.ichns),:)];           % 1317 chns (12 guard chans)
%        s.drad  = [s.drad,  g.raDecon(xd(s.ichns),:)];      %
%        s.csolz = [s.csolz, g.csolzen'];
%        s.asolz = [s.asolz, g.asolzen'];
%        s.td    = [s.td,    g.tdiff'];
        s.csolz(n1:n2) = g.csolzen;
	s.asolz(n1:n2) = g.asolzen;
        s.td(n1:n2)    = g.tdiff;
%        s.dsn   = [s.dsn,   g.dist'];
        s.dsn(n1:n2)   = g.dist;
	s.clat(n1:n2)  = g.clat;      s.clon(n1:n2) = g.clon;
	s.alat(n1:n2)  = g.alat;      s.alon(n1:n2) = g.alon;
	s.ctim(n1:n2)  = g.ctime;     s.atim(n1:n2) = g.atime;
	s.cifov(n1:n2) = g.cifov;
	s.cland(n1:n2) = g.clandfrac; s.aland(n1:n2) = g.alandfrac;
%        s.clat  = [s.clat,  g.clat'];      s.clon = [s.clon,  g.clon'];
%        s.alat  = [s.alat,  g.alat'];      s.alon = [s.alon,  g.alon'];
%        s.ctim  = [s.ctim,  g.ctime'];     s.atim = [s.atim,  g.atime'];
      else
        disp(['skip: ra .ne. raDecon. file: ' num2str(ifnum)]);
      end
  end
  fprintf(1,'%03d ',ifnum);
  n1 = n1 + fs(ifnum);
end
fprintf(1,'\n');
clear g rc rd;
[nx ny] = size(s.crad);
fprintf(1,'number of SNO pairs: %d\n', ny);

% ------------- Quality Control ---------------------
disp(['Removing outliers']);
acbias = single(s.crad - s.drad);
clear gx;
for i=1:length(s.ichns)
   n  = single(remove_6sigma(acbias(i,:)));
   nn = single(remove_6sigma(acbias(i,n)));
   gx(i).n = n(nn);
end
% Now find unique set of bad SNO samples
ux = [];
[~, psz] = size(acbias);
for i=1:length(s.ichns)
   ux = [ux setdiff(1:psz,gx(i).n)];
end
ux  = unique(ux);
s.ing = single(setdiff(1:psz,ux));
disp(['  ' num2str(numel(ux)) ' outliers removed']);
clear gx n nn ux psz acbias;

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
%% ----------------------------------------------------------------------- %%

% Convert rad to BT, do global stats (also used in next section)  
% remove bad data first:
idx = s.ing;
cbt = single( real(rad2bt(s.fc(xc(s.ichns)), s.crad(:,idx))) );
dbt = single( real(rad2bt(s.fd(xd(s.ichns)), s.drad(:,idx))) );
cbm = nanmean(cbt,2);
dbm = nanmean(dbt,2);

radstd   = nanstd( s.drad(:,idx) - s.crad(:,idx), 0, 2 );
cdbm     = 0.5*( nanmean(dbt,2) + nanmean(cbt,2) );
  mdr    = 1E-3*( 1./drdbt(s.fc(xc(s.ichns)),cdbm') );
btstd    = mdr.*radstd';  
stderr   = btstd'./sqrt(ny);
whos idx cbt dbt crad drad cbm dbm radstd btstd stderr

%{
figure;plot(s.fd(s.ichns), cbm,'b-',s.fd(s.ichns),dbm,'g-');grid on;
figure;plot(s.fd(s.ichns), cbm - dbm, 'b-', s.fd(s.ichns),stderr,'g-');grid on;
%}
% -------------------------- quantiles --------------------------- %
disp('working on quantiles...');

% create the scene bins for each channel
clear qsBins qxBins qdBins qx qd;
qxBins.B = quantile(cbt,s.prf,2);
qdBins.B = quantile(dbt,s.prf,2);
qx       = cell2mat(struct2cell(qxBins));
qd       = cell2mat(struct2cell(qdBins));
qsBins   = (qx + qd)/2.0;                        % x:AIRS & d:IASI-to-AIRS

% populate the scene bins for each channel (jj)
% NB the s.crad, s.drad have not been subset.
clear binsz btbias radstd cdbm btstd btser bias250 num_cbin num_dbin;
for jj = 1:numel(s.ichns)
  sbins = qsBins(jj,:);
  clear dbin dbinStd dbinN dbinInd xbin xbinStd xbinN xbinInd ubinInd;
  [dbin dbinStd dbinN dbinInd] = Math_bin(dbt(jj,:),cbt(jj,:)-dbt(jj,:),sbins); 
  [cbin cbinStd cbinN cbinInd] = Math_bin(cbt(jj,:),cbt(jj,:)-dbt(jj,:),sbins);

  % diagnostics: record the separate number samples in each bin:
  num_cbin(jj,:) = cbinN;
  num_dbin(jj,:) = dbinN;
  
  for i = 1:length(dbin)                                                       
    ubinInd(i,:) = {union(dbinInd{i},cbinInd{i})};                  
  end

  for i = 1:length(dbin)
    binsz(jj,i)   = length(ubinInd{i});
    btbias(jj,i)  = nanmean( dbt(jj,ubinInd{i}) - cbt(jj,ubinInd{i}) );
    radstd(jj,i)  = nanstd( s.drad(jj,idx(ubinInd{i})) - s.crad(jj,idx(ubinInd{i})) );
    cdbm(i)    = 0.5*( nanmean(dbt(jj,ubinInd{i})) + nanmean(cbt(jj,ubinInd{i})) );
      mdr      = 1E-3*( 1./drdbt(s.fd(jj),cdbm(i)) );
    btstd(jj,i)   = mdr.*radstd(jj,i);  
    btser(jj,i)   = btstd(jj,i)./sqrt(binsz(jj,i));
    %%bias250(jj,i) = btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  jtot  = sum(binsz(jj,:));
  jmdr  = 1E-3*( 1./drdbt(s.fd(jj),cbm(jj)) );
  jbtse = jmdr.* nanstd(s.drad(jj,idx) - s.crad(jj,idx),1,2) / sqrt(jtot);
  fprintf(1,'.');
end
fprintf(1,'\n');

% parameter fitting section
% -------------------------
wmstats        = struct;
wmstats.cbm    = cbm;
wmstats.dbm    = dbm;
wmstats.wn     = s.fd(xd(s.ichns));
wmstats.stderr = stderr;

for jj = 1:numel(s.ichns)
  clear junk;
  wmstats.bias{jj}  = btbias(jj,:);
  wmstats.btser{jj} = btser(jj,:);
  wmstats.binqa{jj} = qsBins(jj,1:end-1);
  wmstats.binsz{jj} = single(binsz(jj,:));

  % weighted mean & stats section
  % -----------------------------
  % 1. full range
  % ----------
  qlo(jj,1)  = 1;
  qhi(jj,1)  = 200;
  jtot = sum(binsz(jj,:));
  for i = 1:length(dbin)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  wmstats.b(jj,1)   = nansum(junk);  
  wmstats.mx(jj,1)  = nanmax(btbias(jj,:));  
  wmstats.mn(jj,1)  = nanmin(btbias(jj,:));
  wmstats.blo(jj,1) = qsBins(jj,1);
  wmstats.bhi(jj,1) = qsBins(jj,end);
  wmstats.ph(jj,1)  = qsBins(jj, find(btbias(jj,:) == nanmax(btbias(jj,:)),1) );
  wmstats.pl(jj,1)  = qsBins(jj, find(btbias(jj,:) == nanmin(btbias(jj,:)),1) );
  clear junk;
  % -------------------------------
  % 2. range set by min sample size
  % -------------------------------
  % range select by bin size > 500 samples
  inband    = find(binsz(jj,:) > 500);
  qlo(jj,2) = min(inband); 
  qhi(jj,2) = max(inband);                                      % was min(max(inband), find(qsBins(jj,:) > 297,1));
  jtot = sum(binsz(jj,qlo(jj,2):qhi(jj,2)));
  for i = qlo(jj,2):qhi(jj,2)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  wmstats.b(jj,2)   = nansum(junk);  clear junk;
  wmstats.mx(jj,2)  = nanmax(btbias(jj,qlo(jj,2):qhi(jj,2)));  
  wmstats.mn(jj,2)  = nanmin(btbias(jj,qlo(jj,2):qhi(jj,2)));
  wmstats.blo(jj,2) = qsBins(jj,qlo(jj,2));
  wmstats.bhi(jj,2) = qsBins(jj,qhi(jj,2));
  wmstats.ph(jj,2)  = qsBins(jj,find(btbias(jj,:) == nanmax(btbias(jj,qlo(jj,2):qhi(jj,2))),1) ); 
  wmstats.pl(jj,2)  = qsBins(jj,find(btbias(jj,:) == nanmin(btbias(jj,qlo(jj,2):qhi(jj,2))),1) );
  % ------------------------------------
  % range set by max std.err & bin size.
  % ------------------------------------
  inband = find(btser(jj,:) < 0.04 & binsz(jj,:) > 500);        % highly tuned by trial n error
  if(numel(inband) < 2) fprintf(1,'ichn: %d\t Std Err too large\n',jj); continue; end
  qlo(jj,3) = min(inband);                                      % prob OK for all wns. 
  qhi(jj,3) = max(inband);
  jtot = sum(binsz(jj,qlo(jj,3):qhi(jj,3)));
  for i = qlo(jj,3):qhi(jj,3)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  wmstats.b(jj,3)   = nansum(junk);  clear junk;
  wmstats.mx(jj,3)  = nanmax(btbias(jj,qlo(jj,3):qhi(jj,3)));  
  wmstats.mn(jj,3)  = nanmin(btbias(jj,qlo(jj,3):qhi(jj,3)));
  wmstats.blo(jj,3) = qsBins(jj,qlo(jj,3));
  wmstats.bhi(jj,3) = qsBins(jj,qhi(jj,3));
  wmstats.ph(jj,3)  = qsBins(jj,find(btbias(jj,:) == nanmax(btbias(jj,qlo(jj,3):qhi(jj,3))),1) );
  wmstats.pl(jj,3)  = qsBins(jj,find(btbias(jj,:) == nanmin(btbias(jj,qlo(jj,3):qhi(jj,3))),1) );
    
end
wmstats.qlo = qlo;
wmstats.qhi = qhi;

% save file
% ---------
savfn = ['AC_jplSNO_2013_qnstats_' cband '.mat'];
fprintf(1,'Saving: %s\n',savfn);

%{
save(['/home/chepplew/data/sno/airs_cris/' savfn],'s','wmstats','-v7.3');

% display summary of results for selected channel
% -----------------------------------------------
find(s.fc(s.ichns) > 900,1);
jj = 402;
sfnam = fieldnames(wmstats);
disp([wmstats.wn(jj)]);
% display stats for selected wavenumber. (the three methods should be similar)
for i = 9:length(sfnam)
  %disp([sfnam{i}  wmstats.(sfnam{i})(jj,:)]);
  fprintf(1,'%s    \t%8.4f\t%8.4f\t%8.4f\n', sfnam{i}, wmstats.(sfnam{i})(jj,:) );
end

%
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

%}
