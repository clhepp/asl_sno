function plot_sno_iasi_cris_asl_mat(r, s)

% first run:  s = load_sno_iasi_cris_asl_mat(sdate, edate, xchns, src);
% then run:   r = stats_sno_iasi_cris_asl_mat(2,'band');
%
% currently no neighbour FOVs
% Manual edits required for which missions (IASI-1,2 CrIS-1,2) and resolution.
%
% 
% plot options
% set(gcf,'Resize','off');

addpath /home/chepplew/gitLib/asl_sno/source    

fi   = s.fi(s.ichns);
fc   = s.fc(s.cchns);
fd   = s.fd(s.dchns);
%vers = strrep(s.vers,'_','.');
vers = '';
res  = r.res;
src  = r.src;    %  (s/craft npp = 1, j01=2);
band = r.band;   %  LW, MW or SW

if(strcmp(upper(res),'HIGH')) CR='HR'; end
if(strcmp(upper(res),'LOW'))  CR='LR'; end
phome = ['/home/chepplew/projects/sno/iasi_cris/' CR '/figs/'];

xyr  = year(datetime(datenum(r.sdate),'convertfrom','datenum'));
cyr  = num2str(xyr);
smn  = month(datetime(datenum(r.sdate),'convertfrom','datenum'),'shortname');
smn  = lower(cell2mat(smn));
emn  = month(datetime(datenum(r.edate),'convertfrom','datenum'),'shortname');
emn  = lower(cell2mat(emn));
part = '';

% Initialization
cc=fovcolors;       % Howard's 9 line colors uas as: plot(...,'-','color',cc(i,:))

pfnam_pref = ['sno_i' num2str(src(1)) 'c' num2str(src(2)) '_' lower(CR) '_' lower(band) ...
              '_' cyr smn '-' emn '_'];

wnbnd = [floor(fc(1)-10) ceil(fc(end)+10)];

figure(1);clf;plot(fi, r.ibm,'-', fc, r.cbm,'-', fd, r.dbm,'-');
  grid on;legend('IASI','CrIS','I2C','Location','southEast');
% ------------ Maps ----------------------
if(strcmp(band,'LW'))
  ich = find(fi>900,1);
  cch = find(fc>900,1);
  dch = find(fd>900,1);
  wvn='900';
end
if(strcmp(band,'MW'))
  ich = find(fi>1231,1);
  cch = find(fc>1231,1);
  dch = find(fd>1231,1);
  wvn='1231';
end
if(strcmp(band,'SW'))
  ich = find(r.fi>2410,1);
  cch = find(r.fc>2410,1);
  dch = find(r.fd>2410,1);
  wvn='2410';
end


% ----------------------------------------------------------------
%                     PLOTTING SECTION 
% ----------------------------------------------------------------
%
figure(1);clf;plot(s.fi(s.ichns),r.ibm,'-',s.fc(s.cchns),r.cbm,'-',s.fd(s.dchns),r.dbm,'-');
  grid on;legend('IASI','CrIS','I2C','Location','southEast');
  xlabel('wavenumber cm^{-1}');ylabel('B.T. (K)');title('');
  %saveas(gcf,[phome '2017_iasi2_cris1_mean_bt_lw_spectrum.png'],'png')

% ------------ Maps ----------------------
title2=['ASL I.' num2str(src(1)) 'C.' num2str(src(2)) ' SNO ' r.sdate ' to ' r.edate ' ' ...
         band ' overview'];
fh2=figure(2);clf;set(gcf,'Resize','Off'); set(gcf,'Position',fh2.Position+[0 0 420 240]);
  subplot(2,2,1); simplemap(s.clat, s.clon, s.tdiff*24*60);title('Delay IASI-CrIS mins');
  subplot(2,2,2); simplemap(s.clat, s.clon, s.dist); title('Separation deg');

  junk = r.cbt(cch,:)';
  subplot(2,2,3);simplemap(s.clat(s.iok), s.clon(s.iok), junk);title('CrIS BT (K)');
  subplot(2,2,4);simplemap(s.clat(s.iok), s.clon(s.iok), r.btbias(cch,:)');
    hcb = colorbar;ylabel(hcb,[wvn ' cm^{-1} dBT (K)']);
  title('SNO Bias IASI minus CrIS (K)');
  annotation('textbox', [0 0.9 1 0.1], 'String', title2, 'fontSize',16,...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
  %saveas(gcf,[phome pfnam_pref 'maps.fig'],'fig');


% ------------ Histograms -----------------
lat_pdf = histcounts(s.clat,[-90:2:90]);
figure(3);clf;plot([-89:2:89], lat_pdf,'.-')
   xlabel('latitude');ylabel('population');title('AIRS:CrIS SNO density vs latitude')
   grid on;
   %aslprint([phome pfnam_pref '_pop_vs_lat.png'])
   
pc_diff_pdf = 100.*(r.pdf_cbt - r.pdf_dbt)./r.pdf_cbt;
title3=['ASL I.' num2str(src(1)) ':C.' num2str(src(2)) ' SNO ' r.sdate ' to ' r.edate ' ' wvn ...
        ' cm^{-1} pdfs ' vers];
fh3=figure(3);clf;set(gcf,'Resize','off');set(fh3,'Position',fh3.Position+[0 0 420 240]);
  h1 = subplot(221);plot(r.btcens,r.pdf_cbt(cch,:),'.-', r.btcens,r.pdf_dbt(dch,:),'.-',...
    r.btcens,r.pdf_ibt(ich,:),'-'); grid on;xlim([190 330]);
    xlabel('Scene BT bin (K)');ylabel('Number in bin');legend('CrIS','IASItoCrIS','IASI')
    title('')
  h2=subplot(223);plot(r.biascens, r.pdf_bias(cch,:), '.-');grid on;
    xlabel('bin BT (K)');ylabel('population');legend([wvn ' cm^{-1}']);
    title('IASI:CrIS SNO bias');
  h3=subplot(222);plot(r.btcens,r.pdf_cbt(cch,:),'.-', r.btcens,r.pdf_dbt(dch,:),'.-');
    grid on;xlim([200 320]);legend('CrIS','I2C');xlabel('Tb (K)');
    title('');
  h4=subplot(224);plot(r.btcens, pc_diff_pdf(dch,:),'.-');
    grid on;xlim([200 320]);ylabel('% diff CrIS-IASI');
    title('PDF difference')
  annotation('textbox', [0 0.9 1 0.1], 'String', title3, 'fontSize',16,...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
  %saveas(gcf,[phome  pfnam_pref 'pdfs.fig'],'fig');

   
%figure(2);clf;plot(btcens,pdf_sim_cbt(ich,:),'.-', btcens,pdf_dbt(ich,:),'.-'); grid on;

%figure(2);clf;plot(btcens, 0.25*pdf_nbr_cbt(ich,:),'.-', btcens, pdf_cbt(ich,:),'.-');
%figure(2);clf;plot(btcens, 0.50*pdf_nbr_abt(ich,:),'.-', btcens, pdf_abt(ich,:),'.-');
%figure(2);clf;plot(btcens, 0.50*pdf_nbr_dbt(ich,:),'.-', btcens, pdf_dbt(ich,:),'.-');
%figure(3);clf;plot(biascens, pdf_bias(ich,:), '.-');grid on;

% ------------ fractional PDF differences --------------
%{
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
%}
% ------------ Basic Stats -------------------------
figure(4);clf;plot(fd,r.bias_mn,'-');
%
wnbnd = [floor(fc(1)-10) ceil(fc(end)+10)];
title4=['ASL I.' num2str(src(1)) ':C.' num2str(src(2)) ' SNO ' r.sdate ' to ' r.edate ...
        ' Mean Bias ' band ' ' vers ''];
fh4=figure(4);clf;set(fh4,'Resize','Off');set(fh4,'Position',fh4.Position+[0 0 280 210]);
  h1=subplot(221);plot(fd,r.bias_mn,'-', fd,10*r.btser,'-');
  axis([wnbnd(1) wnbnd(2) -0.8 0.8]);grid on;
  legend('CrIS - IASI','10*std.err.');
  xlabel('wavenumber cm^{-1}');ylabel('CrIS minus IASI (K)');

  h2=subplot(222);hold on; for i=1:9 plot(fc,r.fov(i).mbias,'-','color',cc(i,:)); end
  axis([wnbnd(1) wnbnd(2) -0.8 0.8]);grid on;
  legend('1','2','3','4','5','6','7','8','9','Location','eastOutside',...
         'orientation','vertical');
  xlabel('wavenumber cm^{-1}');ylabel('CrIS minus IASI (K)');
% with FOV 5 as the reference
  h3=subplot(223);hold on;
  for i=[1:4 6:9] plot(fc,r.fov(i).mbias - r.fov(5).mbias,'-','color',cc(i,:)); end
  grid on; axis([wnbnd(1) wnbnd(2) -0.4 0.4]);
  legend('1','2','3','4','6','7','8','9','Location','eastOutside');
  xlabel('wavenumber cm^{-1}');ylabel('dBT (K)');
  title('C minus I rel. FOV 5')
  annotation('textbox', [0 0.9 1 0.1], 'String', title4, 'fontSize',16,...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
  %saveas(gcf,[phome pfnam_pref 'mean_bias.fig'],'fig');


%{  
% Double difference (must have loaded two sets: r1: I:C.1, r2: I:C.2. )
title6=[' SNO dble diff ' r1.sdate ' to ' r1.edate ' Mean Bias ' band ' ' vers ''];
figure(6);clf;hold on; 
  for i=1:9 plot(fc, r1.fov(i).mbias - r2.fov(i).mbias,'-','color',cc(i,:));end
  axis([wnbnd(1) wnbnd(2) -0.4 0.4]);grid on;
    legend('1','2','3','4','5','6','7','8','9',...
           'Location','eastOutside','orientation','vertical');
  xlabel('wavenumber cm^{-1}');ylabel('I.1:C.1 minus I.1:C.2 (K)');grid on;
  title(title6);
  %saveas(gcf, [phome 'sno_i1c1_i1c2_dble_diff_lr_mw_2018feb-jun.fig'],'fig')

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

