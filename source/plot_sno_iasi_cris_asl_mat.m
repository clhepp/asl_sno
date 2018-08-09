function plot_sno_iasi_cris_asl_mat(s, band)

% first run:  s = load_sno_iasi_cris_asl_mat(sdate, edate, xchns, src);
% then run:   r = stats_sno_iasi_cris_asl_mat(2,'band');
%
% currently no neighbour FOVs
% Manual edits required for which missions (IASI-1,2 CrIS-1,2) and resolution.
%
% 

addpath /asl/matlib/aslutil               % drdbt.m

% Hardwire spectral band
if(~ismember(band,{'LW','MW','SW'})) error('Invalid band'); return; end

% plot options
% set(gcf,'Resize','off');

fi   = s.fi(s.ichns);
fc   = s.fc(s.cchns);
fd   = s.fd(s.dchns);
%vers = strrep(s.vers,'_','.');
res  = r.res;
src  = r.src;    %  (s/craft npp = 1, j01=2);
band = r.band;   %  LW, MW or SW

if(strcmp(res,'HIGH')) CR='HR'; end
if(strcmp(res,'LOW'))  CR='LR'; end
phome = ['/home/chepplew/projects/sno/airs_cris/' CR '/figs/'];

xyr  = year(datetime(datenum(r.sdate),'convertfrom','datenum'));
cyr  = num2str(xyr);
xmn  = month(datetime(datenum(r.sdate),'convertfrom','datenum'),'shortname');
xmn  = lower(cell2mat(xmn));
part = 'a';
% Initialization
phome = '/home/chepplew/projects/sno/iasi_cris/figs/';
hamm  = 0;

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
figure(2);clf;simplemap(s.clat(s.iok), s.clon(s.iok), cbt(ich,:)');title('CrIS BT (K)');
figure(2);clf;simplemap(s.clat(s.iok), s.clon(s.iok), (dbt(ich,:) - cbt(ich,:))');
  title('900wn Bias IASI - CrIS (K)');

% ------------ Histograms -----------------
ich = find(s.fi(s.ichns)>900,1); cch = find(s.fc(s.cchns)>900,1); 
dch = find(s.fd(s.dchns)>900,1);
pc_diff_pdf = 100.*(pdf_cbt - pdf_dbt)./pdf_cbt;
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
figure(4);clf;h1=subplot(2,1,1);plot(q.qn(ich,1:end-1), q.btbias(ich,:), '.-');grid on;
  axis([190 330 -1 0.3]);title('2013 ASL AC SNO 902 wn Bias'); ylabel('AIRS - CrIS K')
  h2=subplot(2,1,2);semilogy(q.qn(ich,1:end-1), cbinN,'-');grid on;xlim([190 330]);
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

