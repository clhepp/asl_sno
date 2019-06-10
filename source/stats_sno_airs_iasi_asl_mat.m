function r = stats_sno_airs_iasi_asl_mat(s,band)

% stats_sno_airs_iasi_asl_mat.m
%
% first run:  [s] = load_sno_airs_iasi_asl_mat(sdate1, sdate2, xchns, src);


% Hardwire spectral band
if(~ismember(band,{'LW','MW','SW'})) error('Invalid band'); return; end

% plot options
% set(gcf,'Resize','off');

s.rd = s.ri2a;
fa   = s.fa(s.achns);
fi   = s.fi(s.ichns);
fd   = s.fd(s.dchns);
vers = '';

phome = '/home/chepplew/projects/sno/airs_iasi/figs/';
hamm = 0;

% --------------- convert to BT ----------------------------
disp('converting to BT')
abt      = real(rad2bt(fa, s.ra(:,s.iok)));
if hamm
   junk  = hamm_app(double(s.ri2a(:,s.iok)));
   dbt   = real(rad2bt(fd, junk));
else
  ibt      = abs(rad2bt(fi, s.ri(:,s.iok)));         % when complex
  dbt      = real(rad2bt(fd, s.rd(:,s.iok)));
end
%nbr_cbt  = real(rad2bt(fc, s.nbr_rLW(:,:,ig)));     clear junk;
%nbr_abt  = real(rad2bt(fa, s.nbr_ra(:,:,ig)));
%nbr_dbt  = real(rad2bt(fd, s.nbr_rd(:,:,ig)));
% ---------------- Basic Stats ------------------------------
izx = find(ibt == 0);
ibt(izx) = NaN;

disp('calculating global basic stats');
btbias   = abt - dbt;
 % whos abt cbt dbt nbr_abt nbr_cbt nbr_dbt btbias
abm      = nanmean(abt,2);
ibm      = nanmean(ibt,2);
dbm      = nanmean(dbt,2);
bias_mn  = nanmean(abt - dbt,2);
bias_sd  = nanstd(abt - dbt, 0,2);
radstd   = nanstd( s.ra(:,s.iok) - s.ri2a(:,s.iok),0,2 );
 adbm    = 0.5*( nanmean(abt,2) + nanmean(dbt,2) );
 mdr     = 1E-3*( 1./drdbt(fd,adbm) );
btstd    = mdr.*radstd;
btser    = btstd./sqrt(size(abt,2));

% Work in radiance space (for IASI > 2400 cm-1 radiance becomes -ve due to noise)
r_mn     = 0.5*( s.ra(:,s.iok) + s.rd(:,s.iok));
bt_mn    = rad2bt(fd, r_mn);
r_diff   = s.ra(:,s.iok) - s.rd(:,s.iok);
bias_bt  = bt_mn - rad2bt(fd, r_mn + r_diff);
bias_mn  = nanmean(bias_bt,2);


%nbr_cbm  = nanmean(nbr_cbt,3);
%nbr_abm  = nanmean(nbr_abt,3);
%nbr_dbm  = nanmean(nbr_dbt,3);
%nbr_cbsd = nanstd(nbr_cbt,0,3);                     % global std.dev
%nbr_absd = nanstd(nbr_abt,0,3);
%nbr_dbsd = nanstd(nbr_dbt,0,3);
 whos *bt *bm bt* bias_* btser nbr_*

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

 whos nsd_abt nsd_cbt nsd_dbt nex_cbt pdf_nex_cbt pdf_nsd_cbt;
 
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
  whos sim_rc sim_ra sim_cbt sim_abt;

% --------------------- full set PDFs --------------------------
disp('calculating pdfs');
clear pdf_*
btbins  = [190.0: 1.0: 330]; btcens = [190.5: 1.0: 329.5];
for i=1:size(ibt,1) pdf_ibt(i,:) = histcounts(ibt(i,:), btbins); end
for i=1:size(abt,1) pdf_abt(i,:) = histcounts(abt(i,:), btbins); end
for i=1:size(dbt,1) pdf_dbt(i,:) = histcounts(dbt(i,:), btbins); end
%for i=1:size(nbr_cbt,1) pdf_nbr_cbt(i,:) = histcounts(nbr_cbt(i,:), btbins); end
%for i=1:size(nbr_abt,1) pdf_nbr_abt(i,:) = histcounts(nbr_abt(i,:), btbins); end
%for i=1:size(nbr_dbt,1) pdf_nbr_dbt(i,:) = histcounts(nbr_dbt(i,:), btbins); end
%for i=1:size(sim_rc,1)  pdf_sim_cbt(i,:) = histcounts(sim_cbt(i,:), btbins); end
%%for i=1:size(abt,1)     pdf_sim_abt(i,:) = histcounts(sim_abt(i,:), btbins); end
%
biasbins = [-10:0.05:10];  biascens = [-9.975:0.05:9.975];
for i=1:size(dbt,1) pdf_bias(i,:) = histcounts(abt(i,:)-dbt(i,:), biasbins); end
%
%for i=1:size(nsd_cbt,1) pdf_nsd_cbt(i,:) = histcounts(nsd_cbt(i,:), biasbins); end
 whos pdf_*
 
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
xbt = dbt;      xra = s.ri2a;
%xbt = sim_abt;  xra = sim_ra;

clear qaBins qdBins qsBins qa qd qad;
qaBins.B = quantile(abt,s.prf,2);
qdBins.B = quantile(xbt,s.prf,2);
qa       = cell2mat(struct2cell(qaBins));
qd       = cell2mat(struct2cell(qdBins));
qad      = (qa + qd)/2.0;

for jj = 1:numel(s.achns)
  sbins = qad(jj,:);
  [dbin dbinStd dbinN dbinInd] = Math_bin(xbt(jj,:),abt(jj,:)-xbt(jj,:),sbins); 
  [abin abinStd abinN abinInd] = Math_bin(abt(jj,:),abt(jj,:)-xbt(jj,:),sbins);

  % diagnostics: record the separate number samples in each bin:
  num_abin(jj,:) = abinN;
  num_dbin(jj,:) = dbinN;
  
  for i = 1:length(dbin)                                                       
    ubinInd(i,:) = {union(dbinInd{i},abinInd{i})};                  
  end

  for i = 1:length(dbin)
    q.binsz(jj,i)   = length(ubinInd{i});
    q.btbias(jj,i)  = nanmean( abt(jj,ubinInd{i}) - xbt(jj,ubinInd{i}) );
    q.radstd(jj,i)  = nanstd( s.ra(jj,s.iok(ubinInd{i})) - xra(jj,s.iok(ubinInd{i})) );   % s.rc
    cdbm(i)    = 0.5*( nanmean(abt(jj,ubinInd{i})) + nanmean(xbt(jj,ubinInd{i})) );
      mdr      = 1E-3*( 1./drdbt(fd(jj),adbm(i)) );
    q.btstd(jj,i)   = mdr.*q.radstd(jj,i);  
    q.btser(jj,i)   = q.btstd(jj,i)./sqrt(q.binsz(jj,i));
    %%bias250(jj,i) = q.btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  jtot  = sum(q.binsz(jj,:));
  jmdr  = 1E-3*( 1./drdbt(fd(jj),abm(jj)) );
  jbtse = jmdr.* nanstd(s.ra(jj,s.iok) - xra(jj,s.iok),1,2) / sqrt(jtot);   % s.rc
  fprintf(1,'.');
end
fprintf(1,'\n');
q.qn = qad;
  whos q     % binsz btbias btstd btser

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

disp('working on CrIS FOV subsets');
clear xFOV nFOV fov;
for i = 1:4 xFOV{i} = find(s.ifov(s.iok) == i); end
for i = 1:4 nFOV(i) = numel(xFOV{i}); end
ra_tmp = s.ra(:,s.iok);
rd_tmp = s.rd(:,s.iok);
for i=1:4
  fov(i).r_mn    = 0.5*( ra_tmp(:,xFOV{i}) + rd_tmp(:,xFOV{i}));
  fov(i).r_diff  = ra_tmp(:,xFOV{i}) - rd_tmp(:,xFOV{i});
  fov(i).bias_mn = nanmean(rad2bt(fd, fov(i).r_mn) - ...
                       rad2bt(fd, fov(i).r_mn + fov(i).r_diff), 2);
  fov(i).abt     = real( rad2bt(s.fa(s.achns), ra_tmp(:, xFOV{i})) );
  fov(i).dbt     = real( rad2bt(s.fd(s.dchns), rd_tmp(:, xFOV{i})) );
end
clear ra_tmp rd_tmp;
%
btbins  = [190.0: 1.0: 330]; btcens = [190.5: 1.0: 329.5];

for i=1:4
  for j=1:length(s.achns) fov(i).pdf(j,:) = histcounts(fov(i).ibt(j,:), btbins); end
end

%for i=1:4
%  fov(i).mbias = nanmean(fov(i).bias_bt, 2);
%end

junk = s.ra(:,s.iok) - s.rd(:,s.iok);
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
  xrd = s.ri2a(:,xFOV{i});
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
r.nsam  = size(dbt,2);
r.fa    = fa;         r.fi = fi;         r.fd = fd;
r.achns = s.achns; r.ichns = s.ichns; r.dchns = s.dchns;
r.abt      = abt;
r.btbias   = single(btbias);
r.abm = abm;  r.ibm = ibm;  r.dbm = dbm;
r.bias_mn  = bias_mn;  r.bias_sd = bias_sd; r.btser = btser;  r.btstd = btstd;
r.pdf_abt  = pdf_abt;
r.pdf_ibt  = pdf_ibt;
r.pdf_dbt  = pdf_dbt;
r.btbins   = btbins;      r.btcens  = btcens;
r.pdf_bias = pdf_bias;
r.biasbins = biasbins;  r.biascens = biascens;
r.q        = q;
r.fov      = fov;

%


disp('completed and return structure r');

%{
% ----------------------------------------------------------------
%                     PLOTTING SECTION 
% ----------------------------------------------------------------
%
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
