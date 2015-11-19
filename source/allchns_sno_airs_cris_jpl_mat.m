function [s a] = allchns_sno_airs_cris_jpl_mat(sdate)

% Load and prep stats for all channels from a SMALL SNO set

% INPUT: sdate
%        e.g. '2013/01/01';




if (length(sdate) ~= 10) fprintf(1,'Error in date\n'); exit; end
syr = sdate(1:4);   smn = sdate(6:7);    sdy = sdate(9:10);
nyr = str2num(syr); nmn = str2num(smn);  ndy = str2num(sdy);


xx=load('/asl/data/airs/airs_freq.mat'); fa=xx.freq; clear xx;
xx=load('cris_freq_nogrd.mat'); f_cris=xx.vchan; clear xx;  % 1305 chns (no guard chans)
xx=load('cris_freq_2grd.mat');  fc = xx.vchan; clear xx;    % 1317 chns (12 guard chans)
fd = f_cris([1:1037 1158:1305]);                            % prior knowledge from Howards decon routine

dp     = '/asl/s1/chepplew/projects/sno/airs_cris/v10_0_0/standard/'; % standard/';
snoLst = dir(strcat(dp,'sno_airs_cris_*.mat'));
fprintf(1,'Found %d total SNO files\n',numel(snoLst));

dstart = datenum([nyr nmn ndy]);
for i=1:numel(snoLst)
  junk = snoLst(i).name(15:22);
  thisdat = datenum( [str2num(junk(1:4)) str2num(junk(5:6)) str2num(junk(7:8))] );
  if(thisdat <= dstart) ifn1 = i; end
  %%if(thisdat <= dlast)  ifn2 = i; end
end
ifn2 = ifn1 + 30;                        % only 10 SNO days

fprintf(1,'Processing SNO files from: %s to %s\n',snoLst(ifn1).name, snoLst(ifn2).name);

s.td   = [];   s.arad = [;]; s.crad = [;]; s.drad = [;]; s.ctim = [];  s.atim = []; 
s.arlat = []; s.arlon = [];  s.dsn  = []; s.crlat = []; s.crlon = []; s.csolz = [];  
s.alnfr = []; s.clnfr = [];  s.cifv  = [];
a.nSam = [];  a.avrd = [;]; a.avra = [;]; a.avrc = [;]; a.sdra = [;];  a.sdrc = [;]; 
a.sdrd = [;];

for ifnum = ifn1:ifn2
  vars = whos('-file',strcat(dp,snoLst(ifnum).name));
  if( ismember('raDecon', {vars.name}) )              % AIRS->CrIS is present
    if(snoLst(ifnum).bytes > 1.0E4)
      g = load(strcat(dp,snoLst(ifnum).name));
      s.arad  = [s.arad, g.ra(1:200,:)];                % [arad, [ra(achn,:); avaw]]; etc
      s.crad  = [s.crad, g.rc(1:200,:)];                % 1317 chns (12 guard chans)
      s.drad  = [s.drad, g.raDecon(1:200,:)];           %


    end
  end
  fprintf(1,'.');
end
fprintf(1,'\n');

%  

  abt  = real(rad2bt(fa(1:200),s.arad));
  xbt  = real(rad2bt(fc(1:200),s.crad)); 
  dbt  = real(rad2bt(fd(1:200),s.drad)); 
  xrad = s.crad;            
  whos abt xbt dbt xrad
%


clear qaBins qxBins qdBins qx qa qd qsBins;
qxBins.B = quantile(xbt,prf,2);
qdBins.B = quantile(dbt,prf,2);
qx       = cell2mat(struct2cell(qxBins));
qd       = cell2mat(struct2cell(qdBins));
qsBins   = (qx + qd)/2.0;                        % x:AIRS & d:IASI-to-AIRS


clear binsz btbias radstd dbm btstd btser bias250;
for jj = 1:200
  sbins = btbins; sbins = qsBins(jj,:);
  clear dbin dbinStd dbinN dbinInd xbin xbinStd xbinN xbinInd ubinInd;
  [dbin dbinStd dbinN dbinInd] = Math_bin(dbt(jj,:),dbt(jj,:),sbins); 
  [xbin xbinStd xbinN xbinInd] = Math_bin(xbt(jj,:),xbt(jj,:),sbins);

  for i = 1:length(dbin)                                                       
    ubinInd(i,:) = {union(dbinInd{i},xbinInd{i})};                  
  end

  for i = 1:length(dbin)
    binsz(jj,i)   = length(ubinInd{i});
    btbias(jj,i)  = nanmean( dbt(jj,ubinInd{i}) - xbt(jj,ubinInd{i}) );
    radstd(jj,i)  = nanstd( s.drad(jj,ubinInd{i}) - xrad(jj,ubinInd{i}) );
    dbm(i)     = 0.5*( nanmean(dbt(jj,ubinInd{i})) + nanmean(xbt(jj,ubinInd{i})) );
      mdr      = 1E-3*( 1./drdbt(fd(jj),dbm(i)) );
    btstd(jj,i)   = mdr.*radstd(jj,i);  
    btser(jj,i)   = btstd(jj,i)./sqrt(binsz(jj,i));
    %%bias250(jj,i) = btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  fprintf(1,'.');
end

% parameter fitting section
% -------------------------
wmstats = struct;
for jj = 1:200
  clear junk;
  %%%wmstats.wn(jj) = sd.Wavs(jj);

  % range select by bin size and hot scenes
  inband = find(binsz(jj,:) > 100);
  blo(jj,1) = min(inband); 
  bhi(jj,1) = max(inband);                                      % was min(max(inband), find(qsBins(jj,:) > 297,1));
  % range select by std err - deal with non-monotonic var at tails
  inband = find(btser(jj,:) < 2 & binsz(jj,:) > 50);        % highly tuned by trial n error
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
  wmstats.ph(jj,1) = qsBins(jj, find(btbias(jj,:) == nanmax(btbias(jj,:))) ); % scene Temp of hi bias
  wmstats.pl(jj,1) = qsBins(jj, find(btbias(jj,:) == nanmin(btbias(jj,:))) ); % scene Temp of lo bias
  clear junk;

end

% Plotting section
% ----------------
jj = 90;
bincen = qsBins(jj,1:end-1);
figure(2);clf;plot(bincen,btbias(jj,:),'k.-',bincen,btser(jj,:),'r-');grid on;
figure(2);clf;plot(bincen(blo(jj,1):bhi(jj,1)),btbias(jj,blo(jj,1):bhi(jj,1)));grid on;
figure(2);clf;plot(bincen(blo(jj,1):bhi(jj,1)),btbias(jj,blo(jj,1):bhi(jj,1)));grid on;

