function r = stats_sno_airs_cris_asl_mat(s,band, subset)

% stats_sno_airs_cris_asl_mat.m
%
% first run:  [s] = load_sno_airs_cris_asl_mat(sdate1, sdate2, xchns, res, src, vers);
% INPUTS: s (structire from load)
%         band {'LW','MW','SW'} aid to sampling
%         subset: {'night','day', 'tropics','north','south',''}; (empty = no subset)
% OUTPUT: r. Structure of fields:
%            abm, cbm, dbm. Mean sensor BT
%            arm, crm, drm. Mean sensor radiance.
%            pdf_* .        PDF distros
%            q              structure of quantile statistics.
%            fov            structure of FOV subset statistics.
%            *bias*         Various global stats
%            *chns, f?,     Chanel listings for each sensor.
%            *bins, *cens   Binning values for PDFs.
%            + some other variables to track data source.
%

addpath /home/chepplew/gitLib/asl_sno/source
addpath /asl/packages/airs_decon/source             % hamm_app.m
addpath /asl/matlib/aslutil                         % rad2bt.m
addpath /home/chepplew/projects/cris                % cris_freq*.mat
addpath /home/strow/Git/breno_matlab/Math           % Math_bin.m
addpath /home/chepplew/gitLib/airs_deconv/source    % seq_match.m

% Check spectral band
if(~ismember(band,{'LW','MW','SW'})) error('Invalid band'); return; end

% Check subset
subset = lower(subset);
if(~ismember(subset,{'night','day','tropics','north','south',''})) 
  error('Invalid subset'); return; end

% Copy key variables over
fa   = s.fa(s.achns);
fc   = s.fc;      % (s.cchns);
fd   = s.fd(s.dchns);
vers = strrep(s.vers,'_','.');
res  = upper(s.res);

if(strcmp(res,'LOW'))    CR = 'LR'; end
if(strcmp(res,'MEDIUM')) CR = 'MR'; end
if(strcmp(res,'HIGH'))   CR = 'HR'; end

% ----------------------------------------------------------------------
%                   SUBSET by Chosen Type 
% ----------------------------------------------------------------------
disp('Working on Prime subset');
switch subset
  case 'night'
    iis = find(s.asolz > 90 & s.csolz > 90);
    iix = intersect(s.ig, iis);
  case 'day'
    iis = find(s.asolz <= 90 & s.csolz <= 90);
    iix = intersect(s.ig, iis);
  case 'tropics'
    iis = find(s.cLat > -40 & s.cLat < 40);
    iix = intersect(s.ig, iis);
  case 'north'
    iis = find(s.cLat > 0);
    iix = intersect(s.ig, iis);
  case 'south'
    iis = find(s.cLat < 0);
    iix = intersect(s.ig, iis);
  case ''
    iix = s.ig;
end

% Default is to apply hamming 
hamm = 0;

s.ra     = s.ra(:,iix);
if hamm
   disp('Doing Hamming')
   s.rc  = hamm_app(double(s.rc(:,iix)));
   s.rd  = hamm_app(double(s.ra2c(:,iix)));
else
  s.rc   = s.rc(:,iix);
  s.rd   = s.ra2c(:,iix);
end
s.ra2c = [];     % save some space in memory

% --------------- convert to BT ----------------------------
disp('converting to BT')
abt      = real(rad2bt(fa, s.ra));
cbt      = real(rad2bt(fc, s.rc));
dbt      = real(rad2bt(fd, s.rd));

czs = find(cbt == 0);
if(length(czs)) cbt(czs) = NaN; end

%nbr_cbt  = real(rad2bt(fc, s.nbr_rLW(:,:,ig)));     clear junk;
%nbr_abt  = real(rad2bt(fa, s.nbr_ra(:,:,ig)));
%nbr_dbt  = real(rad2bt(fd, s.nbr_rd(:,:,ig)));

% ---------------- Basic Stats ------------------------------

disp('calculating global basic stats');
%btbias     = dbt - cbt;
abm        = nanmean(abt,2);
cbm        = nanmean(cbt,2);
dbm        = nanmean(dbt,2);
bias_btmn  = nanmean(cbt - dbt,2);
bias_sd    = nanstd(cbt - dbt, 0,2);
radstd     = nanstd( s.rd - s.rc,0,2 );
 cdbm      = 0.5*( nanmean(dbt,2) + nanmean(cbt,2) );
 mdr       = 1.0*( 1./drdbt(fd,cdbm) );           % was 1E-3*
btstd      = mdr.*radstd;
btser      = btstd./sqrt(size(cbt,2));

% Work in radiance space (also handles -ve radiances in SW correctly)
arm      = nanmean(s.ra,2);
crm      = nanmean(s.rc,2);
drm      = nanmean(s.rd,2);
abtrm    = rad2bt(fa, arm);
cbtrm    = rad2bt(fc, crm);
dbtrm    = rad2bt(fd, drm);
% bias_btrm = cbtrm - dbtrm

r_mn     = 0.5*( s.rc + s.rd );
bt_mn    = rad2bt(fd, r_mn);
r_diff   = s.rc - s.rd;
bias_bt  = bt_mn - rad2bt(fd, r_mn - r_diff);
bias_mn  = real(nanmean(bias_bt,2));


%nbr_cbm  = nanmean(nbr_cbt,3);
%nbr_abm  = nanmean(nbr_abt,3);
%nbr_dbm  = nanmean(nbr_dbt,3);
%nbr_cbsd = nanstd(nbr_cbt,0,3);
%nbr_absd = nanstd(nbr_abt,0,3);
%nbr_dbsd = nanstd(nbr_dbt,0,3);
% whos *bt *bm bt* bias_* btser nbr_*

% --------------- neighbour stats at SNOs ----------------------------

%nsd_abt  = nanstd(nbr_abt,0,2);
%nsd_cbt  = nanstd(nbr_cbt,0,2);
%nsd_dbt  = nanstd(nbr_dbt,0,2);
%
%nmx_cbt  = max(nbr_cbt,[],2);
%nmn_cbt  = min(nbr_cbt,[],2);
%nex_cbt  = nmx_cbt - nmn_cbt;
%
%inLowSD = find(squeeze(nsd_cbt <= 6));
%inLowEX = find(squeeze(nex_cbt(ich,1,:)) <= 6);
%xbins = [0:1:60];   xcens = [0.5:1:59.5];
%for i=1:length(s.cchns) pdf_nex_cbt(i,:) = histcounts(squeeze(nex_cbt(i,1,:)), xbins); end

% whos nsd_* nex_* pdf_*
 
% ---------------- generate mean vals from neighbours ---------------
%sim_ra   = squeeze(nanmean(s.nbr_ra,2));
%sim_rc   = squeeze(nanmean(s.nbr_rLW,2));
%sim_cbt  = real(rad2bt(s.fc(s.cchn), sim_rc'));
%sim_abt  = real(rad2bt(s.fa(s.achn), sim_ra'));
%

% --------- generate simulated Obs using SNOs & neighbours ---------------
clear sim_*;
%for j=1:length(s.cchns)
%  for i=1:size(s.rc,2) sim_ra(j,i) = 0.8*s.rc(j,i) + 0.1*sum(s.nbr_rLW(j,1:2,i)); end
%end
%sim_abt  = real(rad2bt(s.fc(s.cchns), sim_ra(:,ig)));
%sim_abt  = real(rad2bt(s.fa(s.achn(4)), sim_ra'));
  whos sim_*;

% --------------------- full set PDFs --------------------------
disp('calculating pdfs');
clear pdf
tbtbins  = [190.0: 1.0: 330]; tbtcens = [190.5: 1.0: 329.5];
for i=1:size(cbt,1) 
  [pdf.cbt(i,:),~,pdf.cbin(i,:)] = histcounts(cbt(i,:), tbtbins); end
for i=1:size(abt,1) 
  [pdf.abt(i,:),~,pdf.abin(i,:)] = histcounts(abt(i,:), tbtbins); end
for i=1:size(dbt,1) 
  [pdf.dbt(i,:),~,pdf.dbin(i,:)] = histcounts(dbt(i,:), tbtbins); end
%for i=1:size(nbr_cbt,1) pdf.nbr_cbt(i,:) = histcounts(nbr_cbt(i,:), tbtbins); end
%for i=1:size(nbr_abt,1) pdf.nbr_abt(i,:) = histcounts(nbr_abt(i,:), tbtbins); end
%for i=1:size(nbr_dbt,1) pdf.nbr_dbt(i,:) = histcounts(nbr_dbt(i,:), tbtbins); end
%for i=1:size(sim_rc,1)  pdf.sim_cbt(i,:) = histcounts(sim_cbt(i,:), tbtbins); end
%%for i=1:size(abt,1)     pdf.sim_abt(i,:) = histcounts(sim_abt(i,:), tbtbins); end
%     indexes_of_cbt_with_value_in_bin_70 = find(pdf.cbin(cch,:)==70);
biasbins = [-10:0.05:10];  biascens = [-9.975:0.05:9.975];
for i=1:size(dbt,1) pdf.bias(i,:) = histcounts(dbt(i,:)-cbt(i,:), biasbins); end
%
%for i=1:size(nsd_cbt,1) pdf.nsd_cbt(i,:) = histcounts(nsd_cbt(i,:), biasbins); end
switch band
  case 'LW'
    drad = 1;
    radbins = [10 : drad : 160];
    wvn = 900;
  case 'MW'
    drad = 0.5;
    radbins = [2 : drad : 100];
    wvn = 1231;
  case 'SW'
    drad = 0.04;
    junk = [0.001 : drad : 4];
    %radbins = junk;  % 0.25*(junk.^2);
    wvn = 2456;
    radbins = bt2rad(wvn,tbtbins);
end
radcens = mean([radbins(1:end-1);radbins(2:end)]);
clear *rbtcens
for i=1:length(fa) arbtcens(i,:) = rad2bt(fa(i),radcens); end
for i=1:length(fc) crbtcens(i,:) = rad2bt(fc(i),radcens); end
for i=1:length(fd) drbtcens(i,:) = rad2bt(fd(i),radcens); end

for i=1:size(s.rc,1) pdf.crad(i,:) = histcounts(real(s.rc(i,:)), radbins); end
for i=1:size(s.ra,1) pdf.arad(i,:) = histcounts(s.ra(i,:), radbins); end
for i=1:size(s.rd,1) pdf.drad(i,:) = histcounts(real(s.rd(i,:)), radbins); end

pdf.tbtbins  = tbtbins;
pdf.tbtcens  = tbtcens;
pdf.biasbins = biasbins;
pdf.biascens = biascens;
pdf.radbins  = radbins;
pdf.radcens  = radcens;
% save space
pdf.abin = int16(pdf.abin);
pdf.cbin = int16(pdf.cbin);
pdf.dbin = int16(pdf.dbin);
%{
figure(3);clf; plot(rbtcens, pdf.crad(402,:),'.-')
  hold on;plot(rbtcens, pdf.arad(795,:),'.-', rbtcens,pdf.drad(402,:),'.-')

%}
% pdf
 
% -------------------------- full-set quantiles --------------------------- %
disp('working on quantiles...');

% Get quantile profiler and set to prf.
%load('/home/strow/Matlab/Sno/prob_vector.mat');  yp = p(1:200:end);
load('/home/chepplew/projects/sno/prob_vector61.mat');               % p [1x61]
junk = 0.0:0.1:10; yp = sigmf(junk,[2 4]); clear junk yp; 
junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
junk = [-5:.25:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./4.0;  clear junk y0;
% Choose which profiler to use (goes in prf)
s.prf  = yp;

% create the scene bins for each channel and choose AIRS or simulated AIRS
xrd = real(s.rd);
xrc = s.rc;
%xbt = sim_abt;  xra = sim_ra;

clear bt rd q *binInd num_*
% BT
bt.qcBins   = quantile(cbt,s.prf,2);
bt.qdBins   = quantile(dbt,s.prf,2);
bt.qcd      = (bt.qcBins + bt.qdBins)/2.0;
% Radiance
rd.qcBins   = quantile(xrc,s.prf,2);
rd.qdBins   = quantile(xrd,s.prf,2);
rd.qcd      = (rd.qcBins + rd.qdBins)/2.0;

%
Yedges  = [-10:0.5:10];
dbin_mn = [];
dbin_sd = [];
dbin_n  = [];

%
for jj = 1:numel(fc)
  Bedges  = bt.qcd(jj,:);
  Redges  = rd.qcd(jj,:);
  Bdata   = cbt(jj,:)-dbt(jj,:);
  Rdata   = xrc(jj,:)-xrd(jj,:);
  ibins   = {};

  for ibin=1:numel(bt.qcd(jj,:))-1
    ibins{ibin}=find(dbt(jj,:)>=bt.qcd(jj,ibin) & dbt(jj,:)<bt.qcd(jj,ibin+1) );
  end
  nbins = numel(ibins);

  for ic=1:nbins
    dbin_mn(jj,ic) = nanmean(Bdata(ibins{ic}));
    dbin_sd(jj,ic) = nanstd(Bdata(ibins{ic}));
    dbin_n(jj,ic)  = numel(ibins{ic});
  end

  [dN,~,~,dbinX,dbinY] = histcounts2(dbt(jj,:),cbt(jj,:)-dbt(jj,:),Bedges,Yedges);
  [cN,~,~,cbinX,cbinY] = histcounts2(cbt(jj,:),cbt(jj,:)-dbt(jj,:),Bedges,Yedges);
  [eN,~,~,ebinX,ebinY] = histcounts2(xrd(jj,:),xrc(jj,:)-xrd(jj,:),Redges,Yedges);
  [fN,~,~,fbinX,fbinY] = histcounts2(xrc(jj,:),xrc(jj,:)-xrd(jj,:),Redges,Yedges);

  dindX = [];  cindX = [];
  eindX = [];  findX = [];
  clear ubinInd;
  for i=1:numel(Bedges)
    dindX{i} = find(dbinX == i);
    cindX{i} = find(cbinX == i);
    ubinInd{i} = union(dindX{i},cindX{i});
    eindX{i} = find(ebinX == i);
    findX{i} = find(fbinX == i);
    rbinInd{i} = union(dindX{i},cindX{i});
  end

%{
  [dbin dbinStd dbinN dbinInd] = Math_bin(dbt(jj,:),cbt(jj,:)-dbt(jj,:),bt.qcd(jj,:)); 
  [cbin cbinStd cbinN cbinInd] = Math_bin(cbt(jj,:),cbt(jj,:)-dbt(jj,:),bt.qcd(jj,:));
  num_cbin(jj,:) = cbinN;
  num_dbin(jj,:) = dbinN;
  for i = 1:length(dbin)                                                       
    ubinInd(i,:) = {union(dbinInd{i},cbinInd{i})};                  
  end
  [~, ~, ~, dbinInd] = Math_bin(xrd(jj,:),xrc(jj,:)-xrd(jj,:),rd.qcd(jj,:)); 
  [~, ~, ~, cbinInd] = Math_bin(xrc(jj,:),xrc(jj,:)-xrd(jj,:),rd.qcd(jj,:));
  for i = 1:length(dbin)                                                       
    rbinInd(i,:) = {union(dbinInd{i},cbinInd{i})};                  
  end
%}

  for i = 1:length(dN)
    % radiance
    q.rbinsz(jj,i)  = length(rbinInd{i});
    r_mn       = 0.5*( xrc(jj,rbinInd{i}) + xrd(jj,rbinInd{i}) );
    bt_mn      = rad2bt(fd(jj), r_mn);
    rdelta     = xrd(jj,rbinInd{i}) - xrc(jj,rbinInd{i});
    bias_bt    = bt_mn - rad2bt(fd(jj), r_mn - rdelta);
    q.bias_mn(jj,i)  = real(nanmean(bias_bt,2));
    % Brightness Temperature
    q.bbinsz(jj,i)  = length(ubinInd{i});
    q.btbias(jj,i)  = nanmean( dbt(jj,ubinInd{i}) - cbt(jj,ubinInd{i}) );
    q.radstd(jj,i)  = nanstd( xrd(jj,ubinInd{i}) - xrc(jj,ubinInd{i}) );
    cdbm(i)    = 0.5*( nanmean(dbt(jj,ubinInd{i})) + nanmean(cbt(jj,ubinInd{i})) );
      mdr      = 1.0*( 1./drdbt(fd(jj),cdbm(i)) );          % was 1E-3*
    q.btstd(jj,i)   = mdr.*q.radstd(jj,i);  
    q.btser(jj,i)   = q.btstd(jj,i)./sqrt(q.bbinsz(jj,i));
    %%bias250(jj,i) = q.btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  jtot  = sum(q.bbinsz(jj,:));
  jmdr  = 1E-3*( 1./drdbt(fd(jj),cbm(jj)) );
  jbtse = jmdr.* nanstd(xrd(jj,:) - xrc(jj,:),1,2) / sqrt(jtot);   % s.rc
  fprintf(1,'.');
end
fprintf(1,'\n');
q.rd_qn  = rd.qcd;
q.bt_qn  = bt.qcd;

%  whos q     % binsz btbias btstd btser

% ------------------- Hot Bin Investigation ----------------------
%{
disp('working on hot scenes');
%btbins  = [190.0: 0.2: 330]; btcens = [190.1: 0.2: 329.9];
ich = 12;
  dHot290  = find(dbt(ich,:) > 290 & dbt(ich,:) <= 291);
  cHot290  = find(cbt(ich,:) > 290 & cbt(ich,:) <= 291);
uHot290    = union(dHot290, cHot290);
  dHot304  = find(dbt(ich,:) > 304 & dbt(ich,:) <= 305);
  cHot304  = find(cbt(ich,:) > 304 & cbt(ich,:) <= 305);
uHot304    = union(dHot304, cHot304);
 whos uHot*
nanmean(dbt(ich,uHot290) - cbt(ich,uHot290));
nanmean(dbt(ich,uHot304) - cbt(ich,uHot304));

for i=1:size(cbt,1) pdf.cbt_bin290(i,:) = histcounts(cbt(i,uHot290), btbins); end
for i=1:size(abt,1) pdf.abt_bin290(i,:) = histcounts(abt(i,uHot290), btbins); end
for i=1:size(dbt,1) pdf.dbt_bin290(i,:) = histcounts(dbt(i,uHot290), btbins); end
%
for i=1:size(cbt,1) pdf.cbt_bin304(i,:) = histcounts(cbt(i,uHot304), btbins); end
for i=1:size(abt,1) pdf.abt_bin304(i,:) = histcounts(abt(i,uHot304), btbins); end
for i=1:size(dbt,1) pdf.dbt_bin304(i,:) = histcounts(dbt(i,uHot304), btbins); end
%
%for i=1:size(cbt,1) pdf.nbr_cbt_bin290(i,:) = histcounts(nbr_cbt(i,:,uHot290), btbins); end
%for i=1:size(abt,1) pdf.nbr_abt_bin290(i,:) = histcounts(nbr_abt(i,:,uHot290), btbins); end
%for i=1:size(dbt,1) pdf.nbr_dbt_bin290(i,:) = histcounts(nbr_dbt(i,:,uHot290), btbins); end
%for i=1:size(cbt,1) pdf.nbr_cbt_bin304(i,:) = histcounts(nbr_cbt(i,:,uHot304), btbins); end
%for i=1:size(abt,1) pdf.nbr_abt_bin304(i,:) = histcounts(nbr_abt(i,:,uHot304), btbins); end
%for i=1:size(dbt,1) pdf.nbr_dbt_bin304(i,:) = histcounts(nbr_dbt(i,:,uHot304), btbins); end
%}
% ----------------------------------------------------------------------
%                   SUBSET by CrIS FOV 
% ----------------------------------------------------------------------
disp('working on CrIS FOV subsets');
clear xFOV nFOV fov;
for i = 1:9 xFOV{i} = find(s.cFov(iix) == i); end
for i = 1:9 nFOV(i) = numel(xFOV{i}); end
for i = 1:9
  fov(i).r_mn    = 0.5*( s.rc(:,xFOV{i}) + s.rd(:,xFOV{i}));
  fov(i).r_diff  = s.rc(:,xFOV{i}) - s.rd(:,xFOV{i});
  fov(i).bias_mn = real(nanmean(rad2bt(fd, fov(i).r_mn) - ...
                       rad2bt(fd, fov(i).r_mn - fov(i).r_diff), 2));
  fov(i).cbt     = real( rad2bt(fc, s.rc(:, xFOV{i})) );
  fov(i).dbt     = real( rad2bt(fd, s.rd(:, xFOV{i})) );
end
%clear rc_tmp rd_tmp;
%

for i=1:9
  for j=1:length(s.cchns) 
    fov(i).pdf_crad(j,:) = histcounts(real(s.rc(j,xFOV{i})), radbins); 
  end
  for j=1:length(s.dchns) 
    fov(i).pdf_drad(j,:) = histcounts(real(s.rd(j,xFOV{i})), radbins); 
  end
end

for i = 1:9
  radstd   = nanstd(r_diff(:,xFOV{i}),0,2 );
   cdbm    = 0.5*( nanmean(dbt(:,xFOV{i}),2) + nanmean(cbt(:,xFOV{i}),2) );
   mdr     = 1E-3*( 1./drdbt(fd,cdbm) );
  btstd    = mdr.*radstd;
  fov(i).btser = btstd./sqrt(numel(xFOV{i}));
end

%  -------- Quantile Analysis ---------
for i=1:9
  qn_c = quantile(fov(i).cbt,s.prf,2);
  qn_d = quantile(fov(i).dbt,s.prf,2);
  fov(i).qn = qn_c; % (qn_c + qn_d)/2.0;
end
disp('computing quantiles for each FOV - will take a while!')
for i=1:9 
  fov(i).qbinsz = []; fov(i).qbtbias = []; fov(i).qbtstd = []; fov(i).qbtser = []; 
end
for i=1:9
  xrd = s.rd(:,xFOV{i});
  xrc = s.rc(:,xFOV{i});
for j = 1:size(cbt,1)
  sbins = fov(i).qn(j,:);
  [dbin dbinStd dbinN dbinInd] = ...
      Math_bin(fov(i).dbt(j,:),fov(i).cbt(j,:) - fov(i).dbt(j,:),sbins); 
  [cbin cbinStd cbinN cbinInd] = ...
      Math_bin(fov(i).cbt(j,:),fov(i).cbt(j,:) - fov(i).dbt(j,:),sbins);
  for k = 1:length(dbin)                                                       
    ubinInd(k,:) = {union(dbinInd{k},cbinInd{k})};                  
    fov(i).qbinsz(j,k)  = length(ubinInd{k});
    fov(i).qbtbias(j,k) = nanmean(fov(i).dbt(j,ubinInd{k}) - fov(i).cbt(j,ubinInd{k}) );

    radstd  = nanstd( xrd(j,ubinInd{k}) - xrc(j,ubinInd{k}) );
    cdbm    = 0.5*( nanmean(fov(i).dbt(j,ubinInd{k})) + ...
                    nanmean(fov(i).cbt(j,ubinInd{k})) );
      mdr      = 1E-3*( 1./drdbt(fd(j),cdbm) );
    fov(i).qbtstd(j,k)   = mdr.*radstd;  
    fov(i).qbtser(j,k)   = fov(i).qbtstd(j,k)./sqrt(fov(i).qbinsz(j,k));

  end
end
fprintf(1,'.')
end

% ----------------- choose which variables to return ------------
r.src   = s.src;
r.res   = s.res;
r.vers  = vers;
r.band  = band;
r.sdate = s.sdate;      r.edate = s.edate;
r.nsam  = size(dbt,2);
r.fa    = fa;           r.fc = fc;             r.fd = fd;
r.achns = s.achns;      r.cchns = s.cchns;     r.dchns = s.dchns;
r.cbt   = single(cbt);
r.btbias   = single(dbt - cbt);
r.abm      = abtrm;     r.cbm = cbtrm;         r.dbm = dbtrm;
r.arm      = arm;       r.crm = crm;           r.drm = drm;
%r.abtm    = abm;       r. cbtm = cbm;         r.dbtm = dbm;
r.bias_mn  = bias_mn;   r.bias_sd = bias_sd;   r.btser = btser;  r.btstd = btstd;
r.fov      = fov;
r.pdf      = pdf;
r.q        = q;
r.fov      = fov;

disp('completed and return structure r');

%{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end here for now %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])  


%}
