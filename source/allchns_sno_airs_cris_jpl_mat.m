function wmstats = allchns_sno_airs_cris_jpl_mat(sdate, igrp)

% Load and prep stats for all channels from a SMALL SNO set

% INPUT: sdate calendar date to start collection
%        e.g. sdate='2013/01/02';
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

clearvars s a wmstats g -except sdate igrp;

cd /home/chepplew/gitLib/asl_sno/run
addpath /home/chepplew/gitLib/asl_sno/source
addpath /home/chepplew/gitLib/asl_sno/data       % cris_freq_*.mat
addpath /home/chepplew/myLib/matlib/aslutil       % rad2bt.m
addpath /home/strow/Git/breno_matlab/Math         % Math_bin.m
addpath /asl/matlab2012/aslutil/                  % drdbt.m
addpath /asl/packages/ccast/source                % seq_match.m

% set up the dates to process:
if (length(sdate) ~= 10) fprintf(1,'Error in date\n'); exit; end
syr = sdate(1:4);   smn = sdate(6:7);    sdy = sdate(9:10);
nyr = str2num(syr); nmn = str2num(smn);  ndy = str2num(sdy);

% set up the channels to process (applies to the common grid):
if (igrp < 1 || igrp > 6) fprintf(1,'igrp out of range (1 to 6)\n'); exit; end
ichns = [(igrp-1)*200 + 1:igrp*200];                        % applies to the AIRS to CRIS
if(igrp == 6) ichns = [1001:1185]; end

% load the frequency grids:
xx=load('/asl/data/airs/airs_freq.mat'); fa=xx.freq; clear xx;
xx=load('cris_freq_nogrd.mat'); f_cris=xx.vchan; clear xx;  % 1305 chns (no guard chans)
xx=load('cris_freq_2grd.mat');  fc = xx.vchan; clear xx;    % 1317 chns (12 guard chans)
fd = f_cris([1:1037 1158:1305]);                            % prior knowledge from Howards decon routine
[xd, xc] = seq_match(fd, fc);                               % track indexes for both grids

% Get quantile profiler and set to prf.
%load('/home/strow/Matlab/Sno/prob_vector.mat');  p = p(1:200:end);
load('/home/chepplew/projects/sno/prob_vector61.mat');               % p [1x61]
% Alternate profiler - softer tails than p.
junk = 0.0:0.1:10; yp = sigmf(junk,[2 4]); clear junk; 
junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
% Choose which profiler to use (goes in prf)
prf  = yp;


clear x cc ccs;
dp = '/asl/s1/chepplew/projects/sno/airs_cris/v10_0_0/standard/'; % standard/';
unix(['cd ' dp '; find . -noleaf -type f -name ''sno_airs_cris_2013*.mat'' -printf ''%P\n'' > /tmp/fn.txt;']);
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

dstart = datenum([nyr nmn ndy]);
for i=1:numel(ccs)
  %junk = snoLst(i).name(15:22);
  junk = ccs{i}(15:22);
  thisdat = datenum( [str2num(junk(1:4)) str2num(junk(5:6)) str2num(junk(7:8))] );
  if(thisdat <= dstart) ifn1 = i; end
  %%if(thisdat <= dlast)  ifn2 = i; end
end
ifn2 = ifn1 + 60;                        % only 10 SNO days
%fprintf(1,'Processing SNO files from: %s to %s\n',snoLst(ifn1).name, snoLst(ifn2).name);
fprintf(1,'Processing SNO files from %s to %s\n',ccs{ifn1}, ccs{ifn2});

clear g;
s.td    = [];   s.arad = [;]; s.crad = [;]; s.drad = [;]; s.ctim = [];  s.atim = []; 
s.arlat = [];  s.arlon = [];  s.dsn  = []; s.crlat = []; s.crlon = []; s.csolz = [];  
%s.alnfr = [];  s.clnfr = []; s.cifv  = [];
%a.nSam  = [];   a.avrd = [;]; a.avra = [;]; a.avrc = [;]; a.sdra = [;]; a.sdrc = [;]; 
%a.sdrd  = [;];

for ifnum = ifn1:ifn2
  %vars = whos('-file',strcat(dp,snoLst(ifnum).name));
  vars = whos('-file',strcat(dp,ccs{ifnum}));
  if( ismember('raDecon', {vars.name}) )              % AIRS->CrIS is present
    %if(snoLst(ifnum).bytes > 1.0E4)
      %g = load(strcat(dp,snoLst(ifnum).name));
      g = load(strcat(dp,ccs{ifnum}));
      %%s.arad  = [s.arad, g.ra(achns,:)];                % [arad, [ra(achn,:); avaw]]; etc
      s.crad  = [s.crad, g.rc(xc(ichns),:)];                % 1317 chns (12 guard chans)
      s.drad  = [s.drad, g.raDecon(xd(ichns),:)];           %

    %end
  end
  fprintf(1,'.');
end
fprintf(1,'\n');
[nx ny] = size(s.crad);
fprintf(1,'number of SNO pairs: %d\n', ny);
  
% convert Obs to BT
%%abt  = real(rad2bt(fa(achns),s.arad));
cbt  = real(rad2bt(fc(xc(ichns)),s.crad)); 
dbt  = real(rad2bt(fd(xd(ichns)),s.drad)); 
crad = s.crad;  
cbm  = nanmean(cbt,2); 
dbm  = nanmean(dbt,2);         
whos cbt dbt crad s cbm dbm


% create the scene bins for each channel
clear qaBins qxBins qdBins qx qa qd qsBins;
qxBins.B = quantile(cbt,prf,2);
qdBins.B = quantile(dbt,prf,2);
qx       = cell2mat(struct2cell(qxBins));
qd       = cell2mat(struct2cell(qdBins));
qsBins   = (qx + qd)/2.0;                        % x:AIRS & d:IASI-to-AIRS

% populate the scene bins for each channel (jj)
clear binsz btbias radstd cdbm btstd btser bias250;
for jj = 1:numel(ichns)
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
    radstd(jj,i)  = nanstd( s.drad(jj,ubinInd{i}) - crad(jj,ubinInd{i}) );
    cdbm(i)    = 0.5*( nanmean(dbt(jj,ubinInd{i})) + nanmean(cbt(jj,ubinInd{i})) );
      mdr      = 1E-3*( 1./drdbt(fd(jj),cdbm(i)) );
    btstd(jj,i)   = mdr.*radstd(jj,i);  
    btser(jj,i)   = btstd(jj,i)./sqrt(binsz(jj,i));
    %%bias250(jj,i) = btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  fprintf(1,'.');
end
fprintf(1,'\n');

% parameter fitting section
% -------------------------
wmstats = struct; blo =[]; bhi = [];
wmstats.cbm = cbm;
wmstats.dbm = dbm;
wmstats.wn  = fd(xd(ichns));

for jj = 1:numel(ichns)
  clear junk;
  wmstats.bias{jj}  = btbias(jj,:);
  wmstats.btser{jj} = btser(jj,:);
  wmstats.bins{jj}  = qsBins(jj,1:end-1);

  % range select by bin size and hot scenes
  inband = find(binsz(jj,:) > 500);
  blo(jj,1) = min(inband); 
  bhi(jj,1) = max(inband);                                      % was min(max(inband), find(qsBins(jj,:) > 297,1));
  % range select by std err - deal with non-monotonic var at tails
  inband = find(btser(jj,:) < 0.04 & binsz(jj,:) > 500);        % highly tuned by trial n error
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
savfn = ['sno_AC_wmstats_2013_chns200x' sprintf('%02d',igrp) '.mat'];
fprintf(1,'Saving: %s\n',savfn);
save(savfn,'wmstats');

% display summary of results for selected channel
% -----------------------------------------------
jj = 90;
sfnam = fieldnames(wmstats);
disp([wmstats.wn(jj)]);
for i = 7:numel(sfnam)
  disp([wmstats.(sfnam{i})(jj,:)]);
end

% Plotting section
% ----------------
jj = 90;
bincen = qsBins(jj,1:end-1);
figure(1);clf;
  subplot(2,1,1);plot(bincen,btbias(jj,:),'k.-',bincen,btser(jj,:),'r-');grid on;
    axis([200 300 -2 2]);
  subplot(2,1,2);semilogy(bincen,binsz(jj,:));axis([200 300 50 20000]);grid on;
figure(2);clf;subplot(2,1,1)
  plot(bincen(blo(jj,1):bhi(jj,1)),btbias(jj,blo(jj,1):bhi(jj,1)));grid on;
 subplot(2,1,2);plot(bincen(blo(jj,2):bhi(jj,2)),btbias(jj,blo(jj,2):bhi(jj,2)));grid on;

figure(3);clf;subplot(2,1,1);plot(wmstats.wn,wmstats.cbm,'b-',wmstats.wn,wmstats.dbm,'g-');grid on;
  subplot(2,1,2);plot(wmstats.wn,wmstats.cbm - wmstats.dbm,'m.-');grid on;
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
