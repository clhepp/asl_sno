 function wmstats = sno_quantile(sd,ans)

% function sno_quantila(sd,ans);
%
% Synopsis: sd: structure of arrays retuend from the reader.
%           ans: for IASI or CRIS paired with AIRS
%                possible values are: 'CRIS'; 'IASI';
%
%
%
%
%
%
%
%
%

addpath /home/chepplew/myLib/matlib/aslutil       % rad2bt.m
addpath /home/strow/Git/breno_matlab/Math         % Math_bin.m
addpath /asl/matlab2012/aslutil/                  % drdbt.m

% Check and process input string 'ans' for pairing with AIRS
if(~strcmp(ans, 'CRIS') && ~strcmp(ans, 'IASI')) fprintf(1, 'ERROR: wrong ans\n'); end
if(strcmp(ans, 'IASI')) IASI = 1;  CRIS = 0; end
if(strcmp(ans, 'CRIS')) CRIS = 1;  IASI = 0; end

% Get quantile profiler and set to prf.
%load('/home/strow/Matlab/Sno/prob_vector.mat');  p = p(1:200:end);
load('/home/chepplew/projects/sno/prob_vector61.mat');               % p [1x61]
% Alternate profiler - softer tails than p.
junk = 0.0:0.1:10; yp = sigmf(junk,[2 4]); clear junk; 
junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
% Choose which profiler to use (goes in prf)
prf  = yp;

% for bias referenced to 250K black body
rad250 = bt2rad(sd.Wavs,250);             % reference for scaling to 250 K
drd250 = 1000.0*drdbt(sd.Wavs,250);       % convert to K.mW-1.m-2.sr-1.cm.


% For IASI have to apply basic QA to remove zero valued radiances
% test for zero in the long wave channel only but apply to all.
if(IASI)
  inaZ = find(sd.arad(1,:) <= 0.0); iniZ = find(sd.irad(1,:) <= 0.0); indZ = find(sd.drad(1,:) <= 0.0);
    fprintf('%d, %d, %d\n',numel(inaZ),numel(iniZ),numel(indZ));
  inaB = find(sd.arad > 199); iniB = find(sd.irad > 199); indB = find(sd.drad > 199);
    fprintf('%d, %d, %d\n',numel(inaB),numel(iniB),numel(indB));
  inQ = find(sd.iqual(1,:) == 2);
  % Note: Standard set: only iniZ are significant
  %       Tropical subset: all three need to be set:
  arad_orig = sd.arad; irad_orig = sd.irad; drad_orig = sd.drad;
  for k = 1:numel(sd.Wavs)
    %arad(inaZ) = NaN; irad(inaZ) = NaN; drad(inaZ) = NaN;
    sd.arad(k,iniZ) = NaN; sd.irad(k,iniZ) = NaN; sd.drad(k,iniZ) = NaN;
    sd.arad(k,inQ) = NaN;  sd.irad(k,inQ) = NaN;  sd.drad(k,inQ) = NaN;
  end
end
% Basic QA for sd.crad and sd.drad
if(CRIS)
  indZ = find(sd.drad(1,:) <= 0.001); incZ = find(sd.crad(1,:) <= 0.001);
    fprintf(1,'%d, %d\n', numel(indZ), numel(incZ));
  for k = 1:numel(sd.Wavs)
    sd.drad(k,indZ) = NaN; sd.crad(k,indZ) = NaN;
    sd.drad(k,incZ) = NaN; sd.crad(k,incZ) = NaN;
  end
end 


if(CRIS)
  % subset if needed 
  inD  = find(sd.csolz < 90);    inN = find(sd.csolz >= 90);
  abt  = real(rad2bt(sd.Wavs,sd.arad));
  xbt  = real(rad2bt(sd.Wavs,sd.crad)); 
  dbt  = real(rad2bt(sd.Wavs,sd.drad)); 
  xrad = sd.crad;            
end
if(IASI) 
  xbt  = real(rad2bt(sd.Wavs,sd.arad)); 
  ibt  = real(rad2bt(sd.Wavs,sd.irad));
  dbt  = real(rad2bt(sd.Wavs,sd.drad)); 
  xrad = sd.arad;      
end

%{
binA   = 180;   binZ = 334;    Dbin = 2;
btbins = binA:Dbin:binZ;
bincen = (binA+0.5*Dbin):Dbin:(binZ-0.5*Dbin);
%}

clear qaBins qxBins qdBins qx qa qd qsBins;
qxBins.B = quantile(xbt,prf,2);
qdBins.B = quantile(dbt,prf,2);
qx       = cell2mat(struct2cell(qxBins));
qd       = cell2mat(struct2cell(qdBins));
qsBins   = (qx + qd)/2.0;                        % x:AIRS & d:IASI-to-AIRS

clear binsz btbias radstd dbm btstd btser bias250;
%clear btbias radstd dbm btstd btser binsz bias250;

for jj = 1:numel(sd.Wavs)
  if(exist('btbins')) sbins = btbins; end
  if(exist('qsBins')) sbins = qsBins(jj,:); end
  clear dbin dbinStd dbinN dbinInd xbin xbinStd xbinN xbinInd ubinInd;
  [dbin dbinStd dbinN dbinInd] = Math_bin(dbt(jj,:),dbt(jj,:),sbins); 
  [xbin xbinStd xbinN xbinInd] = Math_bin(xbt(jj,:),xbt(jj,:),sbins);

  for i = 1:length(dbin)                                                       
    ubinInd(i,:) = {union(dbinInd{i},xbinInd{i})};                  
  end

  for i = 1:length(dbin)
    binsz(jj,i)   = length(ubinInd{i});
    btbias(jj,i)  = nanmean( dbt(jj,ubinInd{i}) - xbt(jj,ubinInd{i}) );
    radstd(jj,i)  = nanstd( sd.drad(jj,ubinInd{i}) - xrad(jj,ubinInd{i}) );
    dbm(i)     = 0.5*( nanmean(dbt(jj,ubinInd{i})) + nanmean(xbt(jj,ubinInd{i})) );
      mdr      = 1E-3*( 1./drdbt(sd.Wavs(jj),dbm(i)) );
    btstd(jj,i)   = mdr.*radstd(jj,i);  
    btser(jj,i)   = btstd(jj,i)./sqrt(binsz(jj,i));
    bias250(jj,i) = btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  fprintf(1,'.');
end

%{
k = 1;
figure(2);clf;subplot(2,1,1);plot(qsBins(k,1:end-1),btbias(k,:),'.-');grid on;
 axis([200 320 -0.7 0.7]);
 subplot(2,1,2);semilogy(qsBins(k,1:end-1),binsz(k,:),'.-');grid on;
%}

% parameter fitting section
% -------------------------
wmstats = struct;
for jj = 1:numel(sd.Wavs)
  clear junk;
  wmstats.wn(jj)    = sd.Wavs(jj);
  wmstats.bias{jj}  = btbias(jj,:);
  wmstats.btser{jj} = btser(jj,:);
  wmstats.binbt{jj} = qsBins(jj,1:end-1);


  % range select by bin size and hot scenes
  inband = find(binsz(jj,:) > 1000);
  blo(jj,1) = min(inband); 
  bhi(jj,1) = max(inband);                                      % was min(max(inband), find(qsBins(jj,:) > 297,1));
  % range select by std err - deal with non-monotonic var at tails
  inband = find(btser(jj,:) < 0.04 & binsz(jj,:) > 900);     % highly tuned by trial n error
  blo(jj,2) = min(inband);                                      % prob OK for all wns. 
  bhi(jj,2) = max(inband);
  %%%bhi(2) = min(max(inband), find(qsBins(jj,:) > 297,1));     % depends on wn.

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
  % selection 1
  % -----------
  jtot = sum(binsz(jj,blo(jj,1):bhi(jj,1)));
  for i = blo(jj,1):bhi(jj,1)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  wmstats.b(jj,2)  = nansum(junk);  clear junk;
  wmstats.mx(jj,2) = nanmax(btbias(jj,blo(jj,1):bhi(jj,1)));  
  wmstats.mn(jj,2) = nanmin(btbias(jj,blo(jj,1):bhi(jj,1)));
  wmstats.lo(jj,2) = qsBins(blo(jj,1));
  wmstats.hi(jj,2) = qsBins(bhi(jj,1));
  wmstats.ph(jj,2) = qsBins( find(btbias(jj,:) == nanmax(btbias(jj,blo(jj,1):bhi(jj,1)))) ); % scene Temp of hi bias
  wmstats.pl(jj,2) = qsBins( find(btbias(jj,:) == nanmin(btbias(jj,blo(jj,1):bhi(jj,1)))) ); % scene Temp of lo bias
  % selection 2
  % -----------
  jtot = sum(binsz(jj,blo(jj,2):bhi(jj,2)));
  for i = blo(jj,2):bhi(jj,2)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  wmstats.b(jj,3)  = nansum(junk);  clear junk;
  wmstats.mx(jj,3) = nanmax(btbias(jj,blo(jj,2):bhi(jj,2)));  
  wmstats.mn(jj,3) = nanmin(btbias(jj,blo(jj,2):bhi(jj,2)));
  wmstats.lo(jj,3) = qsBins(blo(jj,2));
  wmstats.hi(jj,3) = qsBins(bhi(jj,2));
  wmstats.ph(jj,3) = qsBins( find(btbias(jj,:) == nanmax(btbias(jj,blo(jj,2):bhi(jj,2)))) ); % scene Temp of hi bias
  wmstats.pl(jj,3) = qsBins( find(btbias(jj,:) == nanmin(btbias(jj,blo(jj,2):bhi(jj,2)))) ); % scene Temp of lo bias

end   % end jj


%{
% Plotting section
% ----------------
jj = 1;
figure(2);clf;plot(bincen(blo(jj,1):bhi(jj,1)),btbias(jj,blo(jj,1):bhi(jj,1)));grid on;


  if(IASI) titStr = sprintf('AIRS IASI SNO %s wn', sd.sWavs{jj}); 
  figNam = sprintf('./figs/AirsIasi_Sno_qBTBias_%swn',sd.sWavs{jj}); end
  if(CRIS) titStr = sprintf('AIRS CrIS SNO %s wn', sd.sWavs{jj}); 
  figNam = sprintf('./figs/AirsCris_Sno_qBTBias_%swn',sd.sWavs{jj}); end
  
  bincen = qsBins(jj,1:end-1); xlims=[floor(min(bincen)) ceil(max(bincen))];
  figure('Visible','On');clf;
    h1=subplot(2,1,1);plot(bincen,btbias(jj,:),'k.-',bincen,btser(jj,:),'r-',bincen,-btser(jj,:),'r-');
    title(titStr);axis([xlims -0.8 0.5]);grid on;ylabel('BT Bias (K)');
    if(IASI) legend('I2A - Airs','Std Err','location','North'); end
    if(CRIS) legend('C2A - Airs','Std Err','location','North'); end
    h2=subplot(2,1,2);semilogy(bincen,binsz(jj,:),'m.-');axis([xlims 10 10^5]);grid;
    legend('Bin population','location','South');xlabel('Scene BT (K)');
    ha=findobj(gcf,'type','axes');set(ha(1),'ylim',[10^2 10^5],'ytick',[10^2 10^3 10^4 10^5]);
    linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp=get(h1,'position');
    set(h1,'position',[pp(1) pp(2)-pp(4)*0.1 pp(3) pp(4)*1.1])
    pp=get(h2,'position'); set(h2,'position',[pp(1) pp(2) pp(3) pp(4)*1.1]);
    %aslprint('./figs/AirsIasi_Sno_qBTBias_901d0wn.png')
    fprintf(1,'saving graphic %d\n', jj)
    %saveas(gcf,[figNam '.fig'],'fig'); saveas(gcf, [figNam '.png'],'png');
  
%end      % end of jj for loop  
%}

end
