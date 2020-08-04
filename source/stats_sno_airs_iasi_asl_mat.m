function r = stats_sno_airs_iasi_asl_mat(s,band)

% stats_sno_airs_iasi_asl_mat.m
%
% first run:  [s] = load_sno_airs_iasi_asl_mat(sdate1, sdate2, xchns, src);


% Hardwire spectral band
if(~ismember(band,{'LW','MW','SW'})) error('Invalid band'); return; end

fa   = s.fa(s.achns);
fi   = s.fi(s.ichns);
fd   = s.fd(s.dchns);
vers = '';

phome = '/home/chepplew/projects/sno/airs_iasi/figs/';
hamm  = 0;

% --------------- convert to BT ----------------------------
disp('converting to BT')
abt      = real(rad2bt(fa, s.ra(:,s.iok)));
if hamm
   s.rd    = hamm_app(double(s.ri2a(:,s.iok)));
   dbt     = real(rad2bt(fd, s.rd));
else
  ibt      = real(rad2bt(fi, s.ri(:,s.iok)));
  s.rd     = s.ri2a(:,s.iok);
  dbt      = real(rad2bt(fd, s.rd));
end
%nbr.cbt  = real(rad2bt(fc, s.nbr_rLW(:,:,ig)));
%nbr.abt  = real(rad2bt(fa, s.nbr_ra(:,:,ig)));
%nbr.dbt  = real(rad2bt(fd, s.nbr_rd(:,:,ig)));
% ---------------- Basic Stats ------------------------------
izx = find(ibt == 0);
ibt(izx) = NaN;

disp('calculating global basic stats');
btbias   = abt - dbt;
abm      = nanmean(abt,2);
ibm      = nanmean(ibt,2);
dbm      = nanmean(dbt,2);
bbias_mn = nanmean(abt - dbt,2);
bias_sd  = nanstd(abt - dbt, 0,2);
radstd   = nanstd( s.ra(:,s.iok) - s.rd,0,2 );
 adbm    = 0.5*( nanmean(abt,2) + nanmean(dbt,2) );
 mdr     = 1E0*( 1./drdbt(fd,adbm) );                 % was 1E-3*()
btstd    = mdr.*radstd;
btser    = btstd./sqrt(size(abt,2));

% Work in radiance space (for IASI > 2400 cm-1 radiance becomes -ve due to noise)
arm      = nanmean(s.ra(:,s.iok),2);
irm      = nanmean(s.ri(:,s.iok),2);
drm      = nanmean(s.rd,2);
r_mn     = 0.5*( s.ra(:,s.iok) + s.rd);
bt_mn    = rad2bt(fd, r_mn);
r_diff   = s.ra(:,s.iok) - s.rd;
bias_bt  = bt_mn - rad2bt(fd, r_mn + r_diff);
rbias_mn = nanmean(bias_bt,2);


%nbr.ibm  = nanmean(nbr.ibt,3);
%nbr.abm  = nanmean(nbr.abt,3);
%nbr.dbm  = nanmean(nbr.dbt,3);
%nbr.ibsd = nanstd(nbr.ibt,0,3);                     % global std.dev
%nbr.absd = nanstd(nbr.abt,0,3);
%nbr.dbsd = nanstd(nbr.dbt,0,3);

% --------------- neighbour stats at SNOs ----------------------------

%nbr.ansd  = nanstd(nbr.abt,0,2);
%nbr.insd  = nanstd(nbr.ibt,0,2);
%nbr.dnsd  = nanstd(nbr.dbt,0,2);
%
%nbr.nmx  = max(nbr.abt,[],2);
%nbr.nmn  = min(nbr.abt,[],2);
%nbr.nex  = nbr.nmx - nbr.nmn;
%
%inLowSD = find(squeeze(nsd_cbt <= 6));
%inLowEX = find(squeeze(nex_cbt(ich,1,:)) <= 6);
%xbins = [0:1:60];   xcens = [0.5:1:59.5];
%for i=1:length(s.cchns) pdf_nex_cbt(i,:) = histcounts(squeeze(nex_cbt(i,1,:)), xbins); end

% whos nsd_abt nsd_cbt nsd_dbt nex_cbt pdf_nex_cbt pdf_nsd_cbt;
 
% ---------------- generate mean vals from neighbours ---------------
%sim.ra   = squeeze(nanmean(s.nbr_ra,2));
%sim.rc   = squeeze(nanmean(s.nbr_rLW,2));
%sim.cbt  = real(rad2bt(s.fc(s.cchn), sim.rc'));
%sim.abt  = real(rad2bt(s.fa(s.achn), sim.ra'));
%

% --------- generate simulated Obs using SNOs & neighbours ---------------
clear sim_rc sim_ra;
%for j=1:length(s.cchns)
%  for i=1:size(s.rc,2) sim_ra(j,i) = 0.8*s.rc(j,i) + 0.1*sum(s.nbr_rLW(j,1:2,i)); end
%end
%sim.abt  = real(rad2bt(s.fc(s.cchns), sim.ra(:,ig)));
%sim.abt  = real(rad2bt(s.fa(s.achn(4)), sim.ra'));

% --------------------- full set PDFs --------------------------
disp('calculating pdfs');
clear pdf
tbtbins  = [190.0: 2.0: 330]; tbtcens = [191.0: 2.0: 329.0];
for i=1:size(ibt,1) 
   [pdf.ibt(i,:),~,pdf.ibin(i,:)] = histcounts(ibt(i,:), tbtbins); end
for i=1:size(abt,1) 
   [pdf.abt(i,:),~,pdf.abin(i,:)] = histcounts(abt(i,:), tbtbins); end
for i=1:size(dbt,1) 
   [pdf.dbt(i,:),~,pdf.dbin(i,:)] = histcounts(dbt(i,:), tbtbins); end
%for i=1:size(nbr.cbt,1) pdf.nbr_cbt(i,:) = histcounts(nbr.cbt(i,:), btbins); end
%for i=1:size(nbr.abt,1) pdf.nbr_abt(i,:) = histcounts(nbr.abt(i,:), btbins); end
%for i=1:size(nbr.dbt,1) pdf.nbr_dbt(i,:) = histcounts(nbr.dbt(i,:), btbins); end
%for i=1:size(sim.rc,1)  pdf.sim_cbt(i,:) = histcounts(sim.cbt(i,:), btbins); end
%%for i=1:size(abt,1)     pdf.sim_abt(i,:) = histcounts(sim.abt(i,:), btbins); end
%
biasbins = [-10:0.05:10];  biascens = [-9.975:0.05:9.975];
for i=1:size(dbt,1) pdf.bbias(i,:) = histcounts(abt(i,:)-dbt(i,:), biasbins); end
%
%for i=1:size(nsd_cbt,1) pdf_nsd_cbt(i,:) = histcounts(nsd_cbt(i,:), biasbins); end
% whos pdf_*

switch band
  case 'LW'
    drad = 1;
    %radbins = [10 : drad : 160];
    wvn = 900;
    radbins = bt2rad(wvn,tbtbins);
  case 'MW'
    drad = 0.5;
    %radbins = [2 : drad : 100];
    wvn = 1231;
    radbins = bt2rad(wvn,tbtbins);
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
for i=1:length(fi) irbtcens(i,:) = rad2bt(fi(i),radcens); end
for i=1:length(fd) drbtcens(i,:) = rad2bt(fd(i),radcens); end

for i=1:size(s.ri,1) pdf.irad(i,:) = histcounts(s.ri(i,s.iok), radbins); end
for i=1:size(s.ra,1) pdf.arad(i,:) = histcounts(s.ra(i,s.iok), radbins); end
for i=1:size(s.rd,1) pdf.drad(i,:) = histcounts(real(s.rd(i,:)), radbins); end

pdf.tbtbins  = tbtbins;
pdf.tbtcens  = tbtcens;
pdf.biasbins = biasbins;
pdf.biascens = biascens;
pdf.radbins  = radbins;
pdf.radcens  = radcens;
% save space
pdf.abin = int16(pdf.abin);
pdf.ibin = int16(pdf.ibin);
pdf.dbin = int16(pdf.dbin);
 
% -------------------------- quantiles --------------------------- %
disp('working on quantiles...');

% Get quantile profiler and set to prf.
%load('/home/strow/Matlab/Sno/prob_vector.mat');  yp = p(1:200:end);
load('/home/chepplew/projects/sno/prob_vector61.mat');               % p [1x61]
junk = 0.0:0.1:10; yp = sigmf(junk,[2 4]); clear junk yp; 
junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
junk = [-5:.25:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./4.0; clear junk y0;
% Choose which profiler to use (goes in prf)
s.prf  = yp;

% create the scene bins for each channel and choose AIRS or simulated AIRS
xrd = real(s.rd);
xra = s.ra(:,s.iok);

clear bt rd q*;
bt.qaBins   = quantile(abt,s.prf,2);
bt.qdBins   = quantile(dbt,s.prf,2);
bt.qad      = (bt.qaBins + bt.qdBins)/2.0;
% radiance
rd.qaBins   = quantile(xra,s.prf,2);
rd.qdBins   = quantile(xrd,s.prf,2);
rd.qad      = (rd.qaBins + rd.qdBins)/2.0;

for jj = 1:numel(s.achns)
  [dbin dbinStd dbinN dbinInd] = Math_bin(dbt(jj,:),abt(jj,:)-dbt(jj,:),bt.qad(jj,:)); 
  [abin abinStd abinN abinInd] = Math_bin(abt(jj,:),abt(jj,:)-dbt(jj,:),bt.qad(jj,:));
  num_abin(jj,:) = abinN;
  num_dbin(jj,:) = dbinN;
  for i = 1:length(dbin)                                                       
    ubinInd(i,:) = {union(dbinInd{i},abinInd{i})};                  
  end

  % radiance space
  [dbin dbinStd dbinN dbinInd] = Math_bin(xrd(jj,:),xra(jj,:)-xrd(jj,:),rd.qad(jj,:)); 
  [abin abinStd abinN abinInd] = Math_bin(xra(jj,:),xra(jj,:)-xrd(jj,:),rd.qad(jj,:));
  num_abin(jj,:) = abinN;
  num_dbin(jj,:) = dbinN;
  for i = 1:length(dbin)                                                       
    rbinInd(i,:) = {union(dbinInd{i},abinInd{i})};                  
  end

  for i = 1:length(dbin)
    % in radiance space
    q.rbinsz(jj,i)  = length(rbinInd{i});
    rdelta          = xra(jj,rbinInd{i}) - xrd(jj,rbinInd{i}); 
    r_mn            = 0.5*( xra(jj,rbinInd{i}) + xrd(jj,rbinInd{i}) );
    bt_mn           = rad2bt(fd(jj), r_mn);
    bias_bt         = bt_mn - rad2bt(fd(jj), r_mn - rdelta);
    q.bias_mn(jj,i) = nanmean(real(bias_bt),2);
    
    % in BT space 
    q.bbinsz(jj,i)  = length(ubinInd{i});
    q.btbias(jj,i)  = nanmean( abt(jj,ubinInd{i}) - dbt(jj,ubinInd{i}) );
    q.radstd(jj,i)  = nanstd( xra(jj,ubinInd{i}) - xrd(jj,ubinInd{i}) );   % s.rc
    adbm(i)    = 0.5*( nanmean(abt(jj,ubinInd{i})) + nanmean(dbt(jj,ubinInd{i})) );
      mdr      = 1E0*( 1./drdbt(fd(jj),adbm(i)) );                % was 1E-3*()
    q.btstd(jj,i)   = mdr.*q.radstd(jj,i);  
    q.btser(jj,i)   = q.btstd(jj,i)./sqrt(q.bbinsz(jj,i));
    %%bias250(jj,i) = q.btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  jtot  = sum(q.bbinsz(jj,:));
  jmdr  = 1E-3*( 1./drdbt(fd(jj),abm(jj)) );
  jbtse = jmdr.* nanstd(xra(jj,:) - xrd(jj,:),1,2) / sqrt(jtot);
  fprintf(1,'.');
end
fprintf(1,'\n');
q.rd_qn = rd.qad;
q.bt_qn = bt.qad;
%  whos q     % binsz btbias btstd btser

% ------------------- Hot Bin Investigation ----------------------
%{
disp('working on hot scenes');
btbins  = [190.0: 0.2: 330]; btcens = [190.1: 0.2: 329.9];
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

clear pdf_*bin*
for i=1:size(cbt,1) pdf_cbt_bin290(i,:) = histcounts(cbt(i,uHot290), btbins); end
for i=1:size(abt,1) pdf_abt_bin290(i,:) = histcounts(abt(i,uHot290), btbins); end
for i=1:size(dbt,1) pdf_dbt_bin290(i,:) = histcounts(dbt(i,uHot290), btbins); end
%
for i=1:size(cbt,1) pdf_cbt_bin304(i,:) = histcounts(cbt(i,uHot304), btbins); end
for i=1:size(abt,1) pdf_abt_bin304(i,:) = histcounts(abt(i,uHot304), btbins); end
for i=1:size(dbt,1) pdf_dbt_bin304(i,:) = histcounts(dbt(i,uHot304), btbins); end
%
%for i=1:size(cbt,1) pdf_nbr_cbt_bin290(i,:) = histcounts(nbr_cbt(i,:,uHot290), btbins); end
%for i=1:size(abt,1) pdf_nbr_abt_bin290(i,:) = histcounts(nbr_abt(i,:,uHot290), btbins); end
%for i=1:size(dbt,1) pdf_nbr_dbt_bin290(i,:) = histcounts(nbr_dbt(i,:,uHot290), btbins); end
%for i=1:size(cbt,1) pdf_nbr_cbt_bin304(i,:) = histcounts(nbr_cbt(i,:,uHot304), btbins); end
%for i=1:size(abt,1) pdf_nbr_abt_bin304(i,:) = histcounts(nbr_abt(i,:,uHot304), btbins); end
%for i=1:size(dbt,1) pdf_nbr_dbt_bin304(i,:) = histcounts(nbr_dbt(i,:,uHot304), btbins); end
%}
%{
figure(5);clf;semilogy(btcens, pdf_dbt_bin290(ich,:),'.-', btcens, pdf_cbt_bin290(ich,:),'.-');
  xlim([275 315]);grid on;hold on;
  plot(btcens, pdf_dbt_bin304(ich,:),'.-', btcens, pdf_cbt_bin304(ich,:),'.-');
  legend('AIRS','CrIS','AIRS','CrIS')

figure(6);clf;semilogy(btcens, pdf_nbr_cbt_bin304(ich,:),'.-');xlim([275 315]);grid on;hold on;
  semilogy(btcens, pdf_nbr_cbt_bin290(ich,:),'.-')
%}
% ----------------------------------------------------------------------
%                   SUBSET by IASI FOV 
% ----------------------------------------------------------------------

disp('working on IASI FOV subsets');
clear xFOV nFOV fov;
for i = 1:4 xFOV{i} = find(s.ifov(s.iok) == i); end
for i = 1:4 nFOV(i) = numel(xFOV{i}); end
%ra_tmp = s.ra(:,s.iok);
%rd_tmp = s.rd(:,s.iok);
for i=1:4
  fov(i).r_mn    = 0.5*( xra(:,xFOV{i}) + xrd(:,xFOV{i}));
  fov(i).r_diff  = xra(:,xFOV{i}) - xrd(:,xFOV{i});
  fov(i).bias_mn = nanmean(rad2bt(fd, fov(i).r_mn) - ...
                       rad2bt(fd, fov(i).r_mn + fov(i).r_diff), 2);
  fov(i).abt     = real( rad2bt(s.fa(s.achns), xra(:, xFOV{i})) );
  fov(i).dbt     = real( rad2bt(s.fd(s.dchns), xrd(:, xFOV{i})) );
end
%

btbins  = [190.0: 1.0: 330]; btcens = [190.5: 1.0: 329.5];
for i=1:4
  for j=1:length(s.dchns) fov(i).pdf(j,:) = histcounts(fov(i).dbt(j,:), btbins); end
end

%for i=1:4
%  fov(i).mbias = nanmean(fov(i).bias_bt, 2);
%end

junk = s.ra(:,s.iok) - s.rd;
for i = 1:4
  radstd  = nanstd(junk(:,xFOV{i}),0,2 );
   adbm    = 0.5*( nanmean(dbt(:,xFOV{i}),2) + nanmean(abt(:,xFOV{i}),2) );
   mdr     = 1E-3*( 1./drdbt(fd,adbm) );
  btstd    = mdr.*radstd;
  fov(i).btser = btstd./sqrt(numel(xFOV{i}));
end
clear junk;
fov(1).btbins = btbins;
fov(1).btcens = btcens;

%  -------- Quantile Analysis ---------
for i=1:4
  qn_a = quantile(fov(i).abt,s.prf,2);
  qn_d = quantile(fov(i).dbt,s.prf,2);
  fov(i).qn = qn_a; % (qn_c + qn_d)/2.0;
end
disp('computing quantiles for each FOV - will take a while!')
for i=1:4 fov(i).binsz = []; fov(i).btbias = []; fov(i).btstd = []; fov(i).btser = []; end
for i=1:4
  xrd = s.rd(:,xFOV{i});
  xra = s.ra(:,xFOV{i});
for j = 1:size(abt,1)
  sbins = fov(i).qn(j,:);
  [dbin dbinStd dbinN dbinInd] = ...
      Math_bin(fov(i).dbt(j,:),fov(i).abt(j,:) - fov(i).dbt(j,:),sbins); 
  [abin abinStd abinN abinInd] = ...
      Math_bin(fov(i).abt(j,:),fov(i).abt(j,:) - fov(i).dbt(j,:),sbins);
  for k = 1:length(dbin)                                                       
    ubinInd(k,:) = {union(dbinInd{k},abinInd{k})};                  
    fov(i).binsz(j,k)  = length(ubinInd{k});
    fov(i).btbias(j,k) = nanmean(fov(i).abt(j,ubinInd{k}) - fov(i).dbt(j,ubinInd{k}) );

    radstd  = nanstd( xra(j,ubinInd{k}) - xrd(j,ubinInd{k}) );   % s.rc
    cdbm    = 0.5*( nanmean(fov(i).abt(j,ubinInd{k})) + ...
                    nanmean(fov(i).dbt(j,ubinInd{k})) );
      mdr      = 1E-3*( 1./drdbt(fd(j),cdbm) );
    fov(i).btstd(j,k)   = mdr.*radstd;  
    fov(i).btser(j,k)   = fov(i).btstd(j,k)./sqrt(fov(i).binsz(j,k));

  end
end
fprintf(1,'.')
end

% ----------------- choose which variables to return ------------
r.src   = s.src;
r.res   = '';  % s.res;
r.vers  = vers;
r.band  = band;
r.sdate = s.sdate;     r.edate = s.edate;
r.nsam  = size(s.iok,2);
r.fa    = fa;         r.fi = fi;         r.fd = fd;
r.achns = s.achns; r.ichns = s.ichns; r.dchns = s.dchns;
r.abt      = single(abt);
r.btbias   = single(btbias);
r.abm = abm;  r.ibm = ibm;  r.dbm = dbm;
r.bias_mn  = rbias_mn;  r.bias_sd = bias_sd; r.btser = btser;  r.btstd = btstd;
r.pdf      = pdf;
r.q        = q;
r.fov      = fov;

%


disp('completed and return structure r');

%{
% ----------------------------------------------------------------
%                     PLOTTING SECTION 
% ----------------------------------------------------------------
%
% plot options
% set(gcf,'Resize','off');

cc=fovcolors;       % Howard's 9 line colors uas as: plot(...,'-','color',cc(i,:))

figure(1);clf;plot(fa,abm,'-',fc,cbm,'-',fd,dbm,'-');
  grid on;legend('AIRS','CrIS','A2C','Location','southEast');
% ------------ Maps ----------------------
figure(1);clf;simplemap(s.cLat, s.cLon, s.tdiff*24*60);title('Delay AIRS-CrIS mins');
figure(1);clf;simplemap(s.cLat, s.cLon, s.dist); title('Separation deg');
ich = 402;
figure(1);clf;simplemap(s.cLat(s.ig), s.cLon(s.ig), cbt(ich,:)');title('CrIS BT (K)');
figure(1);clf;simplemap(s.cLat(s.ig), s.cLon(s.ig), (dbt(ich,:) - cbt(ich,:))');
hcb = colorbar; ylabel(hcb,'900 cm^{-1} dBT (K)');
  title('2018d005 A:C2 SNO Bias AIRS - CrIS-2 (K)');

% ------------ Histograms -----------------
ach = find(fa>900,1); cch = find(fc>900,1); dch = find(fd>900,1);     % LW band
ach = find(fa>1231,1); cch = find(fc>1231,1); dch = find(fd>1231,1);  % MW band
pc_diff_pdf = 100.*(pdf_cbt - pdf_abt)./pdf_cbt;
figure(2);clf;plot(btcens,pdf_cbt(cch,:),'.-', btcens,pdf_dbt(dch,:),'.-',...
      btcens,pdf_abt(ach,:),'-'); grid on;xlim([190 330]);
  xlabel('Scene BT bin (K)');ylabel('Number in bin');legend('CrIS.2','AIRStoCrIS','AIRS')
  title(['2018d021e048 AC2 SNO 900 cm^{-1} ' vers])
  
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

figure(3);clf;plot(biascens, pdf_bias(cch,:), '.-');grid on;
  xlabel('bin BT (K)');ylabel('population');legend('900 cm^{-1}');
  title('2018d005 AIRS:CrIS-2 SNO bias hist');
  
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

% ------------ Simple Bias Stats -------------------------
figure(2);cla;plot(fd,bias_mn,'-');
%
wnbnd = [floor(fc(1)-10) ceil(fc(end)+10)]
figure(3);clf;plot(fc,bias_mn,'-', fc,10*btser,'-');  
  axis([wnbnd(1) wnbnd(2) -0.8 0.8]);grid on;legend('CrIS.2 - AIRS','10*std.err.');
  xlabel('wavenumber cm^{-1}');ylabel('CrIS.2 minus AIRS (K)');
  title(['2018Mar AC2 SNO mean bias MW ' vers]);
  %saveas(gcf,[phome '2018mar_ac2_sno_mean_bias_stderr_mw_v20a.png'],'png');

figure(4);clf;hold on; for i=1:9 plot(fc,fov(i).mbias,'-','color',cc(i,:)); end
  axis([wnbnd(1) wnbnd(2) -0.8 0.8]);grid on;
  legend('1','2','3','4','5','6','7','8','9','Location','north',...
         'orientation','horizontal');
  xlabel('wavenumber cm^{-1}');ylabel('CrIS.2 minus AIRS (K)');
  title(['2018d021e048 AC2 SNO mean bias vs FOV MW ' vers]);
  
% Double difference (have loaded two sets: ac1_fov and ac2_fov)
figure(5);clf;hold on; 
  for i=1:9 plot(fc, r3.fov(i).mbias - r1.fov(i).mbias,'-','color',cc(i,:));end
  axis([wnbnd(1) wnbnd(2) -0.4 0.4]);grid on;
    legend('1','2','3','4','5','6','7','8','9',...
           'Location','north','orientation','horizontal');
  xlabel('wavenumber cm^{-1}');ylabel('A:CrIS.1 minus A:CrIS.2 (K)');
  title([{'2018Jan SNO mean bias of'} {'A:C1 minus A:C2.a2.test1 vs FOV MW'}]);
  %saveas(gcf, [phome '2018Jan_ac1_ac2_sno_mean_bias_dble_diff_lw.png'],'png')

% with FOV 5 as the reference
figure(8);clf;hold on;
  for i=[1:4 6:9] plot(fc,fov(i).mbias - fov(5).mbias,'-'); end
  grid on; axis([645 1100 -0.4 0.4]); legend('1','2','3','4','6','7','8','9');
  xlabel('wavenumber cm^{-1}');ylabel('dBT (K)');
nf7 = figure(7);  set(gcf,'Resize','off');
set(nf7,'Position',nf7.Position+[0,0,280 210]);
h1=subplot(3,1,1);hold on;for i=1:9 plot(fc, r3.fov(i).btser,'-','color',cc(i,:));end;
  grid on;legend('1','2','3','4','5','6','7','8','9',...
                 'Location','north','orientation','horizontal');
  annotation('textbox',[.2 .5 .3 .3],'String','A:C1','FitBoxToText','on');
h2=subplot(3,1,2);hold on;for i=1:9 plot(fc, r2.fov(i).btser,'-','color',cc(i,:));end;
  grid on;legend('A:C2.a2v4ref')
h3=subplot(3,1,3);hold on;for i=1:9 plot(fc, r1.fov(i).btser,'-','color',cc(i,:));end;
  grid on;legend('A:C2.a2test1')
  xlabel('wavenumber cm^{-1}');ylabel('AIRS - CrIS (K)');
  title(h1,'2018Jand021e048 A:C SNO std.err vs FOV 3 sets')
%saveas(gcf,[phome '2018d021e048_ac_sno_stderr_vs_fov_ac1_ac2t1_ac2ref_mw.png'],'png')
  
% ---------------- choose hot subset (or no subset here)
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
ich = find(fc > 900,1); % LW(713) = 402;   % or 16
figure(4);clf;h1=subplot(2,1,1);plot(qsBins(ich,1:end-1), btbias(ich,:), '.-');grid on;
  axis([190 330 -1 0.3]);title('2013 ASL AC SNO 902 wn Bias'); ylabel('AIRS - CrIS K')
  h2=subplot(2,1,2);semilogy(qsBins(ich,1:end-1), cbinN,'-');grid on;xlim([190 330]);
  xlabel('Scene BT at 902 wn (K)')

  % saveas(gcf, [phome  '2013_aslSNO_900wn_bias_vs_scene.png'], 'png');

% --- quantiles subset by FOV  -------
cch = find(fc > 900,1);
figure(5);clf; h1=subplot(2,1,1); hold on;
  for i=1:9 plot(fov(i).qn(cch,1:end-1),fov(i).btbias(cch,:),'.-','color',cc(i,:)); end
  grid on;axis([190 320 -1.0 1.5]);legend('1','2','3','4','5','6','7','8','9',...
     'Location','north','orientation','horizontal');
h2=subplot(2,1,2);
  semilogy(fov(1).qn(cch,1:end-1), fov(1).binsz(cch,:),'.-');grid on;xlim([190 320]);
   legend('bin population')
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

%}
