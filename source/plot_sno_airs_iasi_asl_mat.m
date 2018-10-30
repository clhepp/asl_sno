% plot_sno_airs_iasi_asl_mat.m

% first run:  [s] = load_sno_airs_iasi_asl_mat(sdate, edate, xchns, src);
% and          r  = stats_sno_airs_iasi_asl_mat(s,<band>);
% plot options
%

addpath /home/chepplew/gitLib/asl_sno/source    

fi   = s.fi(s.ichns);
fa   = s.fa(s.achns);
fd   = s.fd(s.dchns);
%vers = strrep(s.vers,'_','.');
vers = '';
res  = '';
src  = r.src;    %  (s/craft npp = 1, j01=2);
band = r.band;   %  LW, MW or SW

if(src == 1) SRC='1'; end
if(src == 2) SRC='2'; end
phome = ['/home/chepplew/projects/sno/airs_iasi/figs/'];

xyr  = year(datetime(datenum(r.sdate),'convertfrom','datenum'));
cyr  = num2str(xyr);
smn  = month(datetime(datenum(r.sdate),'convertfrom','datenum'),'shortname');
smn  = lower(cell2mat(smn));
emn  = month(datetime(datenum(r.edate),'convertfrom','datenum'),'shortname');
emn  = lower(cell2mat(emn));
part = '';

% Initialization
cc=fovcolors;       % Howard's 9 line colors uas as: plot(...,'-','color',cc(i,:))

pfnam_pref = ['sno_ai' num2str(src) '_' lower(band) '_' cyr smn '-' emn '_'];

wnbnd = [floor(fa(1)-10) ceil(fa(end)+10)];

figure(1);clf;plot(fi, r.ibm,'-', fa, r.abm,'-', fd, r.dbm,'-');
  grid on;legend('IASI','AIRS','I2A','Location','southEast');
% ------------ Maps ----------------------
if(strcmp(band,'LW'))
  ich = find(fi>900,1);
  ach = find(fa>900,1);
  dch = find(fd>900,1);
  wvn='900';
end
if(strcmp(band,'MW'))
  ich = find(fi>1231,1);
  ach = find(fa>1231,1);
  dch = find(fd>1231,1);
  wvn='1231';
end
if(strcmp(band,'SW'))
  ich = find(r.fi>2410,1);
  ach = find(r.fa>2410,1);
  dch = find(r.fd>2410,1);
  wvn='2410';
end


% ----------------------------------------------------------------
%                     PLOTTING SECTION 
% ----------------------------------------------------------------
%
figure(1);clf;plot(fa,abm,'-',fc,cbm,'-',fd,dbm,'-');
  grid on;legend('AIRS','CrIS','A2C','Location','southEast');
% ------------ Maps ----------------------
title2=['ASL AI.' num2str(src) ' SNO ' r.sdate ' to ' r.edate ' ' ...
         band ' overview'];
fh2=figure(2);clf;set(gcf,'Resize','Off'); set(gcf,'Position',fh2.Position+[0 0 420 240]);
  subplot(2,2,1); simplemap(s.alat, s.alon, s.tdiff*24*60);title('Delay AIRS-IASI mins');
  subplot(2,2,2); simplemap(s.alat, s.alon, s.dist); title('Separation deg');

  junk = r.abt(ach,:)';
  subplot(2,2,3);simplemap(s.alat(s.iok), s.alon(s.iok), junk);title('AIRS BT (K)');
  subplot(2,2,4);simplemap(s.alat(s.iok), s.alon(s.iok), r.btbias(ach,:)');
    hcb = colorbar;ylabel(hcb,[wvn ' cm^{-1} dBT (K)']);
  title('SNO Bias AIRS minus IASI (K)');
  annotation('textbox', [0 0.9 1 0.1], 'String', title2, 'fontSize',16,...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
  %saveas(gcf,[phome pfnam_pref 'maps.fig'],'fig');


% ------------ Histograms -----------------
lat_pdf = histcounts(s.alat,[-90:2:90]);
figure(3);clf;plot([-89:2:89], lat_pdf,'.-')
   xlabel('latitude');ylabel('population');title('AIRS:IASI SNO density vs latitude')
   grid on;
   %aslprint([phome pfnam_pref '_pop_vs_lat.png'])
   
pc_diff_pdf = 100.*(r.pdf_abt - r.pdf_dbt)./r.pdf_abt;
title3=['ASL AI.' num2str(src) ' SNO ' r.sdate ' to ' r.edate ' ' wvn ...
        ' cm^{-1} pdfs ' vers];
fh3=figure(3);clf;set(gcf,'Resize','off');set(fh3,'Position',fh3.Position+[0 0 420 240]);
  h1 = subplot(221);plot(r.btcens,r.pdf_abt(ach,:),'.-', r.btcens,r.pdf_dbt(dch,:),'.-',...
    r.btcens,r.pdf_ibt(ich,:),'-'); grid on;xlim([190 330]);
    xlabel('Scene BT bin (K)');ylabel('Number in bin');legend('AIRS','IASItoAIRS','IASI')
    title('')
  h2=subplot(223);plot(r.biascens, r.pdf_bias(ach,:), '.-');grid on;
    xlabel('bin BT (K)');ylabel('population');legend([wvn ' cm^{-1}']);
    title('AIRS:IASI SNO bias');
  h3=subplot(222);plot(r.btcens,r.pdf_abt(ach,:),'.-', r.btcens,r.pdf_dbt(dch,:),'.-');
    grid on;xlim([200 300]);legend('AIRS','I2A');xlabel('Tb (K)');
    title('');
  h4=subplot(224);plot(r.btcens, pc_diff_pdf(dch,:),'.-');
    grid on;xlim([200 300]);ylabel('% diff AIRS-IASI');
    title('PDF difference')
  annotation('textbox', [0 0.9 1 0.1], 'String', title3, 'fontSize',16,...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
  %saveas(gcf,[phome  pfnam_pref 'pdfs.fig'],'fig');

  
%{
figure(2);clf;plot(btcens, (pdf_sim_cbt(4,:) - pdf_cbt(4,:))./pdf_cbt(4,:),'.-');hold on;
   plot(btcens, (pdf_sim_abt(4,:) - pdf_abt(4,:))./pdf_abt(4,:),'.-');                      
   xlim([200 330]); grid on; legend('CrIS','AIRS');
   title('Fraction difference AIRS and CrIS neighbors vs SNO'); 
%}
% ----------------------- Bias Stats -------------------------
figure(4);clf;plot(fd,r.bias_mn,'-');
%
title4=['ASL A:I.' num2str(src) ' SNO ' r.sdate ' to ' r.edate ...
        ' Mean Bias ' band ' ' vers ''];
fh4=figure(4);clf;set(fh4,'Resize','Off');set(fh4,'Position',fh4.Position+[0 0 280 210]);
  h1=subplot(211);plot(fd,r.bias_mn,'-', fd,10*r.btser,'-');
  axis([wnbnd(1) wnbnd(2) -0.8 0.8]);grid on;
  legend('AIRS - IASI','10*std.err.');
  xlabel('wavenumber cm^{-1}');ylabel('AIRS minus IASI (K)');

  h2=subplot(212);hold on; for i=1:4 plot(fa,r.fov(i).mbias,'-','color',cc(i,:)); end
  axis([wnbnd(1) wnbnd(2) -0.8 0.8]);grid on;
  legend('1','2','3','4','Location','best',...
         'orientation','vertical');
  xlabel('wavenumber cm^{-1}');ylabel('AIRS minus IASI (K)');

  annotation('textbox', [0 0.9 1 0.1], 'String', title4, 'fontSize',16,...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
  %saveas(gcf,[phome pfnam_pref 'mean_bias.fig'],'fig');

%{  
% Double difference (must have loaded two sets: r1: A:I.1, r2: A:I.2. )
title6=['A:I SNO dble diff ' r1.sdate ' to ' r1.edate ' Mean Bias ' band ' ' vers ''];
figure(6);clf;hold on; 
  for i=1:4 plot(fa, r1.fov(i).mbias - r2.fov(i).mbias,'-','color',cc(i,:));end
  plot(fa, r1.bias_mn -  r2.bias_mn,'k-')
  axis([wnbnd(1) wnbnd(2) -0.4 0.4]);grid on;
    legend('1','2','3','4','all',...
           'Location','best');
  xlabel('wavenumber cm^{-1}');ylabel('A:I.1 minus A:I.2 (K)');grid on;
  title(title6);

  annotation('textbox', [0 0.9 1 0.1], 'String', title4, 'fontSize',16,...
    'EdgeColor', 'none','HorizontalAlignment', 'center')

  %saveas(gcf, [phome 'sno_ai1_ai2_dble_diff_lw_2018jan-may.fig'],'fig')

nf4 = figure(4);clf;  set(gcf,'Resize','off');
set(nf4,'Position',nf4.Position+[0,0,280 210]);
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
%saveas(gcf,[phome '_dble_diff.fig'],'fig')
%}  



  
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
%{
% remove fake, dead and popping channels
ig = a2c_good_chans(fairs);          % fairs = L1c set: 2645

fa_good = fairs(ig.all);
[ixx ixy] = seq_match(fd, fa_good);


f_fill   = fairs(ig.fill);
fg       = setdiff(fairs, f_fill);
ig_lw    = seq_match(fd,fg);
fg_lw    = fd(ig_lw);

figure(1);clf;plot(fg_lw, r.bias_mn(ig_lw),'-')

ig_lw  = seq_match(fd, sort(unique([fairs(ig.fill); fairs(ig.edgechans); fairs(ig.dead); ...
         fairs(ig.pop); fairs(ig.noise)])) );
rbg_lw = r.bias_mn;
rbg_lw(ig_lw) = NaN;
figure(1);clf;plot(fd, rbg_lw,'-');xlim([640 1110])



%}
