function r = stats_sno_iasi_cris_asl_mat(s, band, subset)

% BIAS is CrIS minus IASI2CRIS
% first run:  [s] = load_sno_iasi_cris_asl_mat(sdate, edate, xchns, src);
% currently no neighbour FOVs
% Manual edits required for which missions (IASAI-1,2 CrIS-1,2) and resolution.
% INPUTS: s (structire from load)
%         band {'LW','MW','SW'} aid to sampling
%         subset: {'night','day', 'tropics','north','south',''}; (empty = no subset)
%
% 

addpath /home/chepplew/gitLib/asl_sno/source
addpath /asl/packages/airs_decon/source             % hamm_app.m
addpath /asl/matlib/aslutil                         % rad2bt.m
addpath /home/strow/Git/breno_matlab/Math           % Math_bin.m
addpath /home/chepplew/gitLib/airs_deconv/source    % seq_match.m

% Check spectral band
if(~ismember(band,{'LW','MW','SW'})) error('Invalid band'); return; end

% Check subset
subset = lower(subset);
if(~ismember(subset,{'night','day','tropics','north','south',''}))
  error('Invalid subset'); return; end

% plot options
% set(gcf,'Resize','off');

fi   = s.fi(s.ichns);
fc   = s.fc;   % (s.cchns);
fd   = s.fd(s.dchns);
%vers = strrep(s.vers,'_','.');
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
    iis = find( s.csolz > 90);   % s.isolz > 90 &
    iix = intersect(s.iok, iis); 
  case 'day'
    iis = find(s.csolz <= 90);    % s.isolz <= 90 & 
    iix = intersect(s.iok, iis);
  case 'tropics'
    iis = find(s.clat > -40 & s.clat < 40);
    iix = intersect(s.iok, iis);
  case 'north'
    iis = find(s.clat > 0);
    iix = intersect(s.iok, iis);
  case 'south'
    iis = find(s.clat < 0);
    iix = intersect(s.iok, iis);
  case ''
    iix = s.iok;
end

% Default: apply hamming
hamm  = 0;

s.ri     = s.ri(:,iix);
s.rd     = s.rd(:,iix);
if hamm
   s.rc  = hamm_app(double(s.rc(:,iix)));
   %s.rd  = hamm_app(double(s.rd));
else
  s.rc   = s.rc(:,iix);
  %s.rd   = s.rd(:,s.iok);
end

% --------------- convert to BT ----------------------------
csz = find(s.rc == 0);
isz = find(s.ri == 0);
if(length(csz)) s.rc(csz) = NaN; end
if(length(isz)) s.ri(isz) = NaN; end

dbt      = real(rad2bt(fd, s.rd));
cbt      = real(rad2bt(fc, s.rc));
ibt      = real(rad2bt(fi, s.ri));

%nbr_cbt  = real(rad2bt(s.fc(s.cchns), s.nbr_rLW(:,:,s.iok)));
%nbr_ibt  = real(rad2bt(s.fi(s.ichns), s.nbr_ri(:,:,s.iok)));
%nbr_dbt  = real(rad2bt(s.fd(s.dchns), s.nbr_rd(:,:,s.iok)));
% ---------------- Basic Stats ------------------------------
btbias     = cbt - dbt;
ibm        = nanmean(ibt,2);
cbm        = nanmean(cbt,2);
dbm        = nanmean(dbt,2);
bias_btmn  = nanmean(cbt - dbt,2);
bias_sd    = nanstd(cbt - dbt, 0,2);

radstd   = nanstd( s.rc - s.rd,0,2 ); 
 cdbm    = 0.5*( nanmean(dbt,2) + nanmean(cbt,2) );
 mdr     = 1E-3*( 1./drdbt(s.fd(s.dchns),cdbm) );    % was 1E-3*()
btstd    = mdr.*radstd;  
btser    = btstd./sqrt(size(cbt,2));
%nbr_cbm  = nanmean(nbr_cbt,3);
%nbr_ibm  = nanmean(nbr_ibt,3);
%nbr_dbm  = nanmean(nbr_dbt,3);
%nbr_cbsd = nanstd(nbr_cbt,0,3);                     % global std.dev
%nbr_absd = nanstd(nbr_abt,0,3);
%nbr_dbsd = nanstd(nbr_dbt,0,3);

% Work in radiance space (also handles -ve radiances in SW correctly)
irm      = nanmean(s.ri,2);
crm      = nanmean(s.rc,2);
drm      = nanmean(s.rd,2);
ibtrm    = rad2bt(fi, irm);
cbtrm    = rad2bt(fc, crm);
dbtrm    = real(rad2bt(fd, drm));

r_mn     = 0.5*( s.rc + s.rd );
bt_mn    = rad2bt(fd, r_mn);
r_diff   = s.rc - s.rd;
bias_bt  = bt_mn - rad2bt(fd, r_mn - r_diff);
bias_mn  = real(nanmean(bias_bt,2));

  whos *bt* bias* *bm bias_* nbr_*

% --------------- neighbour stats at SNOs ----------------------------

%nsd_abt  = nanstd(nbr_abt,0,2);
%nsd_cbt  = nanstd(nbr_cbt,0,2);
%nsd_dbt  = nanstd(nbr_dbt,0,2);
%
%nmx_cbt  = max(nbr_cbt,[],2);
%nmn_cbt  = min(nbr_cbt,[],2);
%nex_cbt  = nmx_cbt - nmn_cbt;
%
%inLowSD = find(squeeze(nsd_cbt <= 6);
%inLowEX = find(squeeze(nex_cbt(ich,1,:)) <= 6);
%xbins = [0:1:60];   xcens = [0.5:1:59.5];
%for i=1:length(s.cchns) pdf_nex_cbt(i,:) = histcounts(squeeze(nex_cbt(i,1,:)), xbins); end

% whos nsd_* nex_cbt pdf_*;
 
% ---------------- generate mean vals from neighbours ---------------
%sim_ra   = squeeze(nanmean(s.nbr_ra,2));
%sim_rc   = squeeze(nanmean(s.nbr_rLW,2));
%sim_cbt  = real(rad2bt(s.fc(s.cchn), sim_rc'));
%sim_abt  = real(rad2bt(s.fa(s.achn), sim_ra'));
%

% --------- generate simulated Obs using SNOs & neighbours ---------------
clear sim_rc sim_ra;
%for j=1:length(s.cchns)
%  for i=1:size(s.rc,2) sim_ra(j,i) = 0.8*s.rc(j,i) + 0.1*sum(s.nbr_rLW(j,1:2,i)); end
%end
%sim_abt  = real(rad2bt(s.fc(s.cchns), sim_ra(:,ig)));
%sim_abt  = real(rad2bt(s.fa(s.achn(4)), sim_ra'));
%  whos sim_rc sim_ra sim_cbt sim_abt;

% --------------------- full set PDFs --------------------------
btbins  = [190.0: 1.0: 330]; btcens = [190.5: 1.0: 329.5];
for i=1:size(cbt,1) pdf_cbt(i,:) = histcounts(cbt(i,:), btbins); end
for i=1:size(ibt,1) pdf_ibt(i,:) = histcounts(ibt(i,:), btbins); end
for i=1:size(dbt,1) pdf_dbt(i,:) = histcounts(dbt(i,:), btbins); end
%for i=1:size(nbr_cbt,1) pdf_nbr_cbt(i,:) = histcounts(nbr_cbt(i,:), btbins); end
%for i=1:size(nbr_abt,1) pdf_nbr_abt(i,:) = histcounts(nbr_abt(i,:), btbins); end
%for i=1:size(nbr_dbt,1) pdf_nbr_dbt(i,:) = histcounts(nbr_dbt(i,:), btbins); end
%for i=1:size(sim_rc,1)  pdf_sim_cbt(i,:) = histcounts(sim_cbt(i,:), btbins); end
%for i=1:size(abt,1)     pdf_sim_abt(i,:) = histcounts(sim_abt(i,:), btbins); end
%
biasbins = [-20:0.2:20];  biascens = [-19.9:0.2:19.9];
for i=1:size(dbt,1) pdf_bias(i,:) = histcounts(dbt(i,:)-cbt(i,:), biasbins); end
%
%for i=1:size(nsd_cbt,1) pdf_nsd_cbt(i,:) = histcounts(nsd_cbt(i,:), biasbins); end

switch band
  case 'LW'
    dtemp = 1;
    radbins = [10 : dtemp : 160];
    wvn = 900;
  case 'MW'
    dtemp = 0.5;
    radbins = [2 : dtemp : 100];
    wvn = 1231;
  case 'SW'
    dtemp = 0.005;
    junk = [0.0 : dtemp : 1];
    radbins = junk.^2;
    wvn = 2360;
end
radcens = [0.5*(radbins(2)+radbins(1)) :dtemp :0.5*(radbins(end-1)+radbins(end))];
rbtbins = rad2bt(wvn,radbins);
rbtcens = rad2bt(wvn,radcens);

for i=1:size(s.rc,1) pdf_crad(i,:) = histcounts(real(s.rc(i,:)), radbins); end
for i=1:size(s.ri,1) pdf_irad(i,:) = histcounts(real(s.ri(i,:)), radbins); end
for i=1:size(s.rd,1) pdf_drad(i,:) = histcounts(real(s.rd(i,:)), radbins); end

 whos pdf*
 
% -------------------------- quantiles --------------------------- %
disp('working on quantiles...');

% Get quantile profiler and set to prf.
%load('/home/strow/Matlab/Sno/prob_vector.mat');  yp = p(1:200:end);
load('/home/chepplew/projects/sno/prob_vector61.mat');               % p [1x61]
junk = 0.0:0.1:10; yp = sigmf(junk,[2 4]); clear junk yp; 
junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
% Choose which profiler to use (goes in prf)
s.prf  = yp;

% create the scene bins for each channel and choose AIRS or simulated AIRS
%xbt = dbt;      xra = s.rd;
%xbt = sim_abt;  xra = sim_ra;

clear qcBins qdBins qsBins qc qd;
qcBins   = quantile(cbt,s.prf,2);
qdBins   = quantile(dbt,s.prf,2);
%qc       = cell2mat(struct2cell(qcBins));
%qd       = cell2mat(struct2cell(qdBins));
qcd      = (qcBins + qdBins)/2.0;

for jj = 1:numel(s.cchns)
  sbins = qcd(jj,:);
  [dbin dbinStd dbinN dbinInd] = Math_bin(dbt(jj,:),cbt(jj,:)-dbt(jj,:),sbins); 
  [cbin cbinStd cbinN cbinInd] = Math_bin(dbt(jj,:),cbt(jj,:)-dbt(jj,:),sbins);

  % diagnostics: record the separate number samples in each bin:
  num_cbin(jj,:) = cbinN;
  num_dbin(jj,:) = dbinN;
  
  for i = 1:length(dbin)                                                       
    ubinInd(i,:) = {union(dbinInd{i},cbinInd{i})};                  
  end

  for i = 1:length(dbin)
    q.binsz(jj,i)   = length(ubinInd{i});
    q.btbias(jj,i)  = nanmean( dbt(jj,ubinInd{i}) - cbt(jj,ubinInd{i}) );
    q.radstd(jj,i)  = nanstd( s.rc(jj,ubinInd{i}) - s.rd(jj,ubinInd{i}) );
    cdbm(i)    = 0.5*( nanmean(dbt(jj,ubinInd{i})) + nanmean(cbt(jj,ubinInd{i})) );
      mdr      = 1E-3*( 1./drdbt(s.fd(s.dchns(jj)),cdbm(i)) );
    q.btstd(jj,i)   = mdr.*q.radstd(jj,i);  
    q.btser(jj,i)   = q.btstd(jj,i)./sqrt(q.binsz(jj,i));
    %%bias250(jj,i) = btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  jtot  = sum(q.binsz(jj,:));
  jmdr  = 1E-3*( 1./drdbt(s.fd(s.dchns(jj)),cbm(jj)) );
  jbtse = jmdr.* nanstd(s.rc(jj,:) - s.rc(jj,:),1,2) / sqrt(jtot);   % s.rc
  fprintf(1,'.');
end
fprintf(1,'\n');
q.qn = qcd;
  whos q

% ----------------------------------------------------------------------
%                   SUBSET by CrIS FOV 
% ----------------------------------------------------------------------
disp('working on CrIS FOV subsets');
clear xFOV nFOV fov;
for i = 1:9 xFOV{i} = find(s.cfov(iix) == i); end
for i = 1:9 nFOV(i) = numel(xFOV{i}); end
for i = 1:9
  fov(i).r_mn    = 0.5*( s.rc(:,xFOV{i}) + s.rd(:,xFOV{i}));
  fov(i).r_diff  = s.rc(:,xFOV{i}) - s.rd(:,xFOV{i});
  fov(i).bias_mn = real(nanmean(rad2bt(fd, fov(i).r_mn) - ...
                       rad2bt(fd, fov(i).r_mn - fov(i).r_diff), 2));
  fov(i).cbt     = real( rad2bt(fc, s.rc(:, xFOV{i})) );
  fov(i).dbt     = real( rad2bt(s.fd(s.dchns), s.rd(:, xFOV{i})) );
end
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

%for i=1:9
%  for j=1:length(s.cchns) fov(i).pdf(j,:) = histcounts(fov(i).cbt(j,:), btbins); end
%end
%
%for i=1:9
%  fov(i).mbias = nanmean(fov(i).cbt - fov(i).dbt, 2);
%end

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

      radstd  = nanstd( xrd(j,ubinInd{k}) - xrc(j,ubinInd{k}) );   % s.rc
      cdbm    = 0.5*( nanmean(fov(i).dbt(j,ubinInd{k})) + ...
                      nanmean(fov(i).cbt(j,ubinInd{k})) );
        mdr   = 1E-3*( 1./drdbt(fd(j),cdbm) );
      fov(i).qbtstd(j,k)   = mdr.*radstd;  
      fov(i).qbtser(j,k)   = fov(i).qbtstd(j,k)./sqrt(fov(i).qbinsz(j,k));

    end
  end
  fprintf(1,'.')
end

% ----------------- choose which variables to return ------------
r.src   = s.src;
r.res   = s.res;
r.vers  = s.vers;
r.band     = band;
r.sdate    = s.sdate;     r.edate = s.edate;
r.nsam     = size(dbt,2);
r.fi       = s.fi;        r.fc = s.fc;              r.fd = s.fd;
r.ichns    = s.ichns;     r.cchns = s.cchns;        r.dchns = s.dchns;
r.cbt      = cbt;
r.btbias   = single(btbias);
r.ibm = ibm;              r.cbm = cbm;              r.dbm = dbm;
r.bias_mn  = bias_mn;     r.bias_sd = bias_sd;      r.btser = btser;  %r.btstd = btstd;
r.irm      = irm;         r.crm     = crm;          r.drm   = drm;

r.fov      = fov;
r.pdf_ibt  = pdf_ibt;     r.pdf_cbt  = pdf_cbt;     r.pdf_dbt  = pdf_dbt;
r.pdf_bias = pdf_bias;
r.pdf_crad = pdf_crad;    r.pdf_irad = pdf_irad;    r.pdf_drad = pdf_drad;
r.radbins  = radbins;     r.radcens = radcens;
r.btbins   = rbtbins;     r.btcens  = rbtcens;
r.biasbins = biasbins;    r.biascens = biascens;
r.q        = q;

%


disp('completed and return structure r');



%{
% ------------------- Hot Bin Investigation ----------------------
btbins  = [190.0: 0.2: 330]; btcens = [190.1: 0.2: 329.9];
ich = 12;
  dHot290  = find(dbt(ich,:) > 290 & dbt(ich,:) <= 291);
  cHot290  = find(cbt(ich,:) > 290 & cbt(ich,:) <= 291);
uHot290    = union(dHot290, cHot290);
  dHot304  = find(dbt(ich,:) > 304 & dbt(ich,:) <= 305);
  cHot304  = find(cbt(ich,:) > 304 & cbt(ich,:) <= 305);
uHot304    = union(dHot304, cHot304);
 whos uHot290 uHot302
nanmean(dbt(ich,uHot290) - cbt(ich,uHot290))
nanmean(dbt(ich,uHot304) - cbt(ich,uHot304))

clear pdf_*bin*
for i=1:size(cbt,1) pdf_cbt_bin290(i,:) = histcounts(cbt(i,uHot290), btbins); end
for i=1:size(abt,1) pdf_abt_bin290(i,:) = histcounts(abt(i,uHot290), btbins); end
for i=1:size(dbt,1) pdf_dbt_bin290(i,:) = histcounts(dbt(i,uHot290), btbins); end
%
for i=1:size(cbt,1) pdf_cbt_bin304(i,:) = histcounts(cbt(i,uHot304), btbins); end
for i=1:size(abt,1) pdf_abt_bin304(i,:) = histcounts(abt(i,uHot304), btbins); end
for i=1:size(dbt,1) pdf_dbt_bin304(i,:) = histcounts(dbt(i,uHot304), btbins); end
%
for i=1:size(cbt,1) pdf_nbr_cbt_bin290(i,:) = histcounts(nbr_cbt(i,:,uHot290), btbins); end
for i=1:size(abt,1) pdf_nbr_abt_bin290(i,:) = histcounts(nbr_abt(i,:,uHot290), btbins); end
for i=1:size(dbt,1) pdf_nbr_dbt_bin290(i,:) = histcounts(nbr_dbt(i,:,uHot290), btbins); end
for i=1:size(cbt,1) pdf_nbr_cbt_bin304(i,:) = histcounts(nbr_cbt(i,:,uHot304), btbins); end
for i=1:size(abt,1) pdf_nbr_abt_bin304(i,:) = histcounts(nbr_abt(i,:,uHot304), btbins); end
for i=1:size(dbt,1) pdf_nbr_dbt_bin304(i,:) = histcounts(nbr_dbt(i,:,uHot304), btbins); end

figure(5);clf;semilogy(btcens, pdf_dbt_bin290(ich,:),'.-', btcens, pdf_cbt_bin290(ich,:),'.-');
  xlim([275 315]);grid on;hold on;
  plot(btcens, pdf_dbt_bin304(ich,:),'.-', btcens, pdf_cbt_bin304(ich,:),'.-');
  legend('AIRS','CrIS','AIRS','CrIS')

figure(6);clf;semilogy(btcens, pdf_nbr_cbt_bin304(ich,:),'.-');xlim([275 315]);grid on;hold on;
  semilogy(btcens, pdf_nbr_cbt_bin290(ich,:),'.-')
%}
%{
% ----------------------------------------------------------------
%                     PLOTTING SECTION 
% ----------------------------------------------------------------
%
figure(1);clf;plot(s.fi(s.ichns),ibm,'-',s.fc(s.cchns),cbm,'-',s.fd(s.dchns),dbm,'-');
  grid on;legend('IASI','CrIS','I2C','Location','southEast');
  xlabel('wavenumber cm^{-1}');ylabel('B.T. (K)');title('');
  %saveas(gcf,[phome '2017_iasi2_cris1_mean_bt_lw_spectrum.png'],'png')
% ------------ Maps ----------------------
figure(2);clf;simplemap(s.clat, s.clon, s.tdiff*24*60);title('Delay IASI-CrIS mins');
figure(2);clf;simplemap(s.clat, s.clon, s.dist); title('Separation deg');
ich = 404;   % find(s.fc(s.cchns) >900,1)
figure(2);clf;simplemap(s.clat(iix), s.clon(iix), cbt(ich,:)');title('CrIS BT (K)');
figure(2);clf;simplemap(s.clat(iix), s.clon(iix), (dbt(ich,:) - cbt(ich,:))');
  title('900wn Bias IASI - CrIS (K)');

% ------------ Histograms -----------------
ich = find(s.fi(s.ichns)>900,1); cch = find(s.fc(s.cchns)>900,1); 
dch = find(s.fd(s.dchns)>900,1);
pc_diff_pdf = 100.*(pdf_cbt - pdf_abt)./pdf_cbt;
figure(3);clf;plot(btcens,pdf_cbt(cch,:),'.-', btcens,pdf_dbt(dch,:),'.-',...
      btcens,pdf_ibt(ich,:),'-'); grid on;xlim([190 310]);
   xlabel('Scene BT bin (K)');ylabel('Number in bin');legend('CrIS','IASI2CrIS','IASI');
   title('2017 ASL I2:C1 SNO 900 cm^{-1} channel')
  %saveas(gcf,[phome '2017_i2c1_sno_900wn_pdf.png'],'png')
   
figure(2);clf;
 h1=subplot(2,1,1);plot(btcens,pdf_cbt(cch,:),'.-', btcens,pdf_dbt(dch,:),'.-'); grid on;
 title('Obs count by 900cm-1 temp. bin');axis([200 330 0 Inf]);
 h2=subplot(2,1,2);plot(btcens,pdf_cbt(cch,:),'.-', btcens,pdf_dbt(dch,:),'.-'); grid on;
   xlim([200 320]);legend('CrIS','AIRS');xlabel('Tb (K)');
 h2=subplot(2,1,2);plot(btcens, pc_diff_pdf(ich,:),'.-');grid on;ylabel('% diff CrIS-AIRS');

figure(2);clf;plot(btcens,pdf_sim_cbt(ich,:),'.-', btcens,pdf_dbt(ich,:),'.-'); grid on;

figure(2);clf;plot(btcens, 0.25*pdf_nbr_cbt(ich,:),'.-', btcens, pdf_cbt(ich,:),'.-');
figure(2);clf;plot(btcens, 0.50*pdf_nbr_abt(ich,:),'.-', btcens, pdf_abt(ich,:),'.-');
figure(2);clf;plot(btcens, 0.50*pdf_nbr_dbt(ich,:),'.-', btcens, pdf_dbt(ich,:),'.-');

figure(3);clf;plot(biascens, pdf_bias(ich,:), '.-');grid on;

% ------------ fractional PDF differences --------------
pc_diff_nbr = 100.*(pdf_nbr_cbt - pdf_nbr_abt)./pdf_nbr_cbt;
figure(2);clf;
 h1=subplot(2,1,1);plot(btcens,pdf_nbr_cbt(4,:),'.-', btcens,pdf_nbr_abt(4,:),'.-'); grid on;
 title('Nbr. Obs count by 900cm-1 temp. bin');axis([200 330 0 Inf]);
 %h2=subplot(2,1,2);plot(btcens,pdf_cbt(4,:),'.-', btcens,pdf_abt(4,:),'.-'); grid on;xlim([300 320]);
 %legend('CrIS','AIRS');xlabel('Tb (K)');
 h2=subplot(2,1,2);plot(btcens, pc_diff_nbr(4,:),'.-');grid on;ylabel('% diff CrIS-AIRS');
 xlim([200 330])

figure(2);clf;plot(btcens, (pdf_sim_cbt(4,:) - pdf_cbt(4,:))./pdf_cbt(4,:),'.-');hold on;
   plot(btcens, (pdf_sim_abt(4,:) - pdf_abt(4,:))./pdf_abt(4,:),'.-');                      
   xlim([200 330]); grid on; legend('CrIS','AIRS');
   title('Fraction difference AIRS and CrIS neighbors vs SNO'); 

% ------------ Basic Stats -------------------------
figure(1);clf;plot(s.fc(s.cchns),bias_mn,'-', s.fc(s.cchns),btser,'-');
  axis([645 1100 -0.8 0.8]);grid on;legend('CrIS - IASI','std.err.');
  xlabel('wavenumber cm^{-1}');ylabel('CrIS minus IASI (K)');
  title('2017 I2:C1 SNO mean bias LW'); 
  %saveas(gcf,[phome '2017_i2c1_sno_mean_bias_stderr_lw.png'],'png');
  
% choose hot subset (or no subset here)
 idx = ':';  % idx = uHot305;
 junk = squeeze(nbr_cbt(4,idx,:));
 dmn4 = dmn_btLW(4,idx);
 clear diff4;  diff4 = zeros(length(dmn4),4);
 for i=1:length(dmn4) diff4(i,:) = junk(i,:) - dmn4(i); end
clear junk; junk = reshape(diff4,[],4*length(dmn4));
xbins = [-40:0.2:40];  xcens = [-39.9:0.2:39.9];
pdf_diff4 = histcounts(junk, xbins);
pdf_nbr1  = histcounts(nbr_cbt(4,:,:),  btbins);
pdf_nbr2  = histcounts(nbr_cbt(4,uHot290,:),  btbins);
pdf_nbr3  = histcounts(nbr_cbt(4,uHot305,:),  btbins);
pdf_sno   = histcounts(cbt(4,:),  btbins);
pdf_sno2  = histcounts(cbt(4,uHot290), btbins);
pdf_sno3  = histcounts(cbt(4,uHot305), btbins);
 
clear pdf_*

% -------------------------- quantiles --------------------------- %
ich = 12;   % or 16
figure(4);clf;h1=subplot(2,1,1);plot(qsBins(ich,1:end-1), btbias(ich,:), '.-');grid on;
  axis([190 330 -1 0.3]);title('2013 ASL AC SNO 902 wn Bias'); ylabel('AIRS - CrIS K')
  h2=subplot(2,1,2);semilogy(qsBins(ich,1:end-1), cbinN,'-');grid on;xlim([190 330]);
  xlabel('Scene BT at 902 wn (K)')

  % saveas(gcf, [phome  '2013_aslSNO_900wn_bias_vs_scene.png'], 'png');
  
% ----------------- Hot Scene Investigation -----------------------
figure(5);clf;simplemap(s.cLat(ig(uHot290)), s.cLon(ig(uHot290)), btbias(uHot290));
  figure(2);plot(dbtcens, pdf_hot290,'+-');grid on;  disp( num2str(nanmean(dbt(uHot290))) )
     disp( num2str(nanstd(dxm_btLW(4,uHot290))) )
figure(1);clf;simplemap(s.cLat(uHot305), s.cLon(uHot305), dbt(uHot305)');
  figure(2);plot(dbtcens, pdf_hot305,'+-');grid on;  disp( num2str(nanmean(dbt(uHot305))) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end here for now %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])  


%btm   = 0.5*(dbtm + cbtm);
%mdr   = 1E-3*(1./drdbt(fd,btm') );
%drse  = sqrt((a.gddrs.^2 + a.gcdrs.^2))/sqrt(sum(a.nSam));
%dbse  = mdr.*drse';
    bias_mn  = nanmean(cbt - dbt,2);
%}
