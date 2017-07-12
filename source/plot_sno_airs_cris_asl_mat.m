% plot_sno_airs_cris_asl_mat.m

% first run:  [s] = load_sno_airs_cris_asl_mat_v2(sdate1, sdate2, xchns, src);

% plot options

phome = '/home/chepplew/projects/sno/airs_cris/LR/';

%btm   = 0.5*(dbtm + cbtm(incd));
%mdr   = 1E-3*(1./drdbt(fd,btm') );
%drse  = sqrt((a.gddrs.^2 + a.gcdrs(incd).^2))/sqrt(sum(a.nSam));
%dbse  = mdr.*drse';

% --- remove 6-sigma outliers -----
addpath /home/chepplew/myLib/matlib/math          % remove_6sigma
disp(['Removing outliers']);
crbias = s.rc - s.rd;
clear g;
for i=1:length(s.cchns)
   n  = remove_6sigma(crbias(i,:));
   nn = remove_6sigma(crbias(i,n));
   g(i).n = n(nn);
end

% Now find unique set of bad SNO samples
x = [];
[~, psz] = size(crbias);
for i=1:length(s.cchns)
   x = [x setdiff(1:psz,g(i).n)];
end
x  = unique(x);
ig = setdiff(1:psz,x);
disp(['Removed ' num2str(length(x)) ' outliers']);

% --------------- convert to BT ----------------------------
abt      = real(rad2bt(s.fa(s.achns), s.ra(:,ig)));
cbt      = real(rad2bt(s.fc(s.cchns), s.rc(:,ig)));
dbt      = real(rad2bt(s.fd(s.dchns), s.rd(:,ig)));
nbr_cbt  = real(rad2bt(s.fc(s.cchns), s.nbr_rLW(:,:,ig)));     clear junk;
nbr_abt  = real(rad2bt(s.fa(s.achns), s.nbr_ra(:,:,ig)));
nbr_dbt  = real(rad2bt(s.fd(s.dchns), s.nbr_rd(:,:,ig)));
btbias   = dbt - cbt;
 whos abt cbt dbt nbr_abt nbr_cbt nbr_dbt btbias
abm      = nanmean(abt,2);
cbm      = nanmean(cbt,2);
dbm      = nanmean(dbt,2);
nbr_cbm  = nanmean(nbr_cbt,3);
nbr_abm  = nanmean(nbr_abt,3);
nbr_dbm  = nanmean(nbr_dbt,3);
nbr_cbsd = nanstd(nbr_cbt,0,3);                     % global std.dev
nbr_absd = nanstd(nbr_abt,0,3);
nbr_dbsd = nanstd(nbr_dbt,0,3);
 whos abm cbm dbm nbr_abm nbr_cbm nbr_dbm nbr_cbsd nbr_absd nbr_dbsd

% --------------- neighbour stats at SNOs ----------------------------

nsd_abt  = nanstd(nbr_abt,0,2);
nsd_cbt  = nanstd(nbr_cbt,0,2);
nsd_dbt  = nanstd(nbr_dbt,0,2);
%
nmx_cbt  = max(nbr_cbt,[],2);
nmn_cbt  = min(nbr_cbt,[],2);
nex_cbt  = nmx_cbt - nmn_cbt;
%
inLowSD = find(squeeze(nsd_cbt <= 6);
inLowEX = find(squeeze(nex_cbt(ich,1,:)) <= 6);
xbins = [0:1:60];   xcens = [0.5:1:59.5];
for i=1:length(s.cchns) pdf_nex_cbt(i,:) = histcounts(squeeze(nex_cbt(i,1,:)), xbins); end

 whos nsd_abt nsd_cbt nsd_dbt nex_cbt pdf_nex_cbt pdf_nsd_cbt;
 
% ---------------- generate mean vals from neighbours ---------------
sim_ra   = squeeze(nanmean(s.nbr_ra,2));
sim_rc   = squeeze(nanmean(s.nbr_rLW,2));
sim_cbt  = real(rad2bt(s.fc(s.cchn), sim_rc'));
sim_abt  = real(rad2bt(s.fa(s.achn), sim_ra'));
%

% --------- generate weighted vals using SNOs & neighbours ---------------
clear sim_rc sim_ra;
for j=1:length(s.cchns)
  for i=1:size(s.rc,2) sim_rc(j,i) = 0.8*s.rc(j,i) + 0.1*sum(s.nbr_rLW(j,1:2,i)); end
end
sim_cbt  = real(rad2bt(s.fc(s.cchns), sim_rc(:,ig)));
%sim_abt  = real(rad2bt(s.fa(s.achn(4)), sim_ra'));
  whos sim_rc sim_ra sim_cbt sim_abt;

% --------------------- full set PDFs --------------------------
btbins  = [190.0: 1.0: 330]; btcens = [190.5: 1.0: 329.5];
for i=1:size(cbt,1) pdf_cbt(i,:) = histcounts(cbt(i,:), btbins); end
for i=1:size(abt,1) pdf_abt(i,:) = histcounts(abt(i,:), btbins); end
for i=1:size(dbt,1) pdf_dbt(i,:) = histcounts(dbt(i,:), btbins); end
for i=1:size(nbr_cbt,1) pdf_nbr_cbt(i,:) = histcounts(nbr_cbt(i,:), btbins); end
for i=1:size(nbr_abt,1) pdf_nbr_abt(i,:) = histcounts(nbr_abt(i,:), btbins); end
for i=1:size(nbr_dbt,1) pdf_nbr_dbt(i,:) = histcounts(nbr_dbt(i,:), btbins); end
for i=1:size(sim_rc,1)  pdf_sim_cbt(i,:) = histcounts(sim_cbt(i,:), btbins); end
%for i=1:size(abt,1)     pdf_sim_abt(i,:) = histcounts(sim_abt(i,:), btbins); end
%
biasbins = [-20:0.2:20];  biascens = [-19.9:0.2:19.9];
for i=1:size(dbt,1) pdf_bias(i,:) = histcounts(dbt(i,:)-cbt(i,:), biasbins); end
%
for i=1:size(nsd_cbt,1) pdf_nsd_cbt(i,:) = histcounts(nsd_cbt(i,:), biasbins); end
 whos pdf_abt pdf_cbt pdf_dbt pdf_nbr_abt pdf_nbr_cbt pdf_nbr_dbt pdf_bias
 
% -------------------------- quantiles --------------------------- %
disp('working on quantiles...');

% Get quantile profiler and set to prf.
%load('/home/strow/Matlab/Sno/prob_vector.mat');  yp = p(1:200:end);
load('/home/chepplew/projects/sno/prob_vector61.mat');               % p [1x61]
junk = 0.0:0.1:10; yp = sigmf(junk,[2 4]); clear junk yp; 
junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
% Choose which profiler to use (goes in prf)
s.prf  = yp;

% create the scene bins for each channel
clear qcBins qdBins qsBins qc qd;
qcBins.B = quantile(cbt,s.prf,2);
qdBins.B = quantile(dbt,s.prf,2);
qc       = cell2mat(struct2cell(qcBins));
qd       = cell2mat(struct2cell(qdBins));
qsBins   = (qc + qd)/2.0;                        % x:AIRS & d:IASI-to-AIRS

clear binsz btbias radstd btstd btser
for jj = 1:numel(s.cchns)
  sbins = qsBins(jj,:);
  [dbin dbinStd dbinN dbinInd] = Math_bin(dbt(jj,:),sim_cbt(jj,:)-dbt(jj,:),sbins); 
  [cbin cbinStd cbinN cbinInd] = Math_bin(cbt(jj,:),sim_cbt(jj,:)-dbt(jj,:),sbins);

  % diagnostics: record the separate number samples in each bin:
  num_cbin(jj,:) = cbinN;
  num_dbin(jj,:) = dbinN;
  
  for i = 1:length(dbin)                                                       
    ubinInd(i,:) = {union(dbinInd{i},cbinInd{i})};                  
  end

  for i = 1:length(dbin)
    binsz(jj,i)   = length(ubinInd{i});
    btbias(jj,i)  = nanmean( dbt(jj,ubinInd{i}) - sim_cbt(jj,ubinInd{i}) );
    radstd(jj,i)  = nanstd( s.rd(jj,ig(ubinInd{i})) - sim_rc(jj,ig(ubinInd{i})) );   % s.rc
    cdbm(i)    = 0.5*( nanmean(dbt(jj,ubinInd{i})) + nanmean(cbt(jj,ubinInd{i})) );
      mdr      = 1E-3*( 1./drdbt(s.fd(s.dchns(jj)),cdbm(i)) );
    btstd(jj,i)   = mdr.*radstd(jj,i);  
    btser(jj,i)   = btstd(jj,i)./sqrt(binsz(jj,i));
    %%bias250(jj,i) = btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  jtot  = sum(binsz(jj,:));
  jmdr  = 1E-3*( 1./drdbt(s.fd(s.dchns(jj)),cbm(jj)) );
  jbtse = jmdr.* nanstd(s.rd(jj,ig) - sim_rc(jj,ig),1,2) / sqrt(jtot);   % s.rc
  fprintf(1,'.');
end
fprintf(1,'\n');
  whos binsz btbias btstd btser

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

% ----------------------------------------------------------------
%                     PLOTTING SECTION 
% ----------------------------------------------------------------
%
figure(1);clf;plot(s.fa(s.achns),abm,'.-',s.fc(s.cchns),cbm,'.-',s.fd(s.dchns),dbm,'.-');
  grid on;legend('AIRS','CrIS','A2C','Location','southEast');
% ------------ Maps ----------------------
figure(1);clf;simplemap(s.cLat, s.cLon, s.tdiff*24*60);title('Delay AIRS-CrIS mins');
figure(1);clf;simplemap(s.cLat, s.cLon, s.dist); title('Separation deg');
ich = 16;
figure(1);clf;simplemap(s.cLat(ig), s.cLon(ig), cbt(ich,:)');title('CrIS BT (K)');
figure(1);clf;simplemap(s.cLat(ig), s.cLon(ig), (dbt(ich,:) - cbt(ich,:))');
  title('Bias AIRS - CrIS (K)');

% ------------ Histograms -----------------
ich = 12;
pc_diff_pdf = 100.*(pdf_cbt - pdf_abt)./pdf_cbt;
figure(2);clf;plot(btcens,pdf_cbt(ich,:),'.-', btcens,pdf_dbt(ich,:),'.-'); grid on;
figure(2);clf;
 h1=subplot(2,1,1);plot(btcens,pdf_cbt(ich,:),'.-', btcens,pdf_dbt(ich,:),'.-'); grid on;
 title('Obs count by 900cm-1 temp. bin');axis([200 330 0 Inf]);
 h2=subplot(2,1,2);plot(btcens,pdf_cbt(ich,:),'.-', btcens,pdf_dbt(ich,:),'.-'); grid on;
   xlim([300 320]);legend('CrIS','AIRS');xlabel('Tb (K)');
 h2=subplot(2,1,2);plot(btcens, pc_diff_pdf(4,:),'.-');grid on;ylabel('% diff CrIS-AIRS');

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
  axis([190 330 -0.8 0.4]);
  h2=subplot(2,1,2);semilogy(qsBins(ich,1:end-1), cbinN,'-');grid on;xlim([190 330]);
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


