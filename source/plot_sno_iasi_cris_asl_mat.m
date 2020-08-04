function plot_sno_iasi_cris_asl_mat(r, s, vis)

% first run:  s = load_sno_iasi_cris_asl_mat(sdate, edate, xchns, src);
% then run:   r = stats_sno_iasi_cris_asl_mat(s,'band');
%
% currently no neighbour FOVs
% Manual edits required for which missions (IASI-1,2 CrIS-1,2) and resolution.
%
% 
% plot options
% set(gcf,'Resize','off');

addpath /home/chepplew/gitLib/asl_sno/source    
addpath /asl/matlib/aslutil                      % simplemap

% Check number of input arguments
if(nargin ~= 3) error('Please enter all 3 input arguments'); return; end
if(~ismember(vis,[0 1])) error('vis must be 0 or 1'); return; end
if(vis == 0) VIS=false; end
if(vis == 1) VIS=true;  end

% close any open figures from previous runs
close all;

% copy over useful variables
fi   = r.fi(r.ichns);
fc   = r.fc(r.cchns);
fd   = r.fd(r.dchns);
vers = strrep(r.vers,'_','.');
%vers = '';
res  = r.res;
src  = r.src;    %  (1st: s/craft npp = 1, j01=2, 2nd: MetOp-A, MetOp-B);
band = r.band;   %  LW, MW or SW

if(strcmp(upper(res),'HIGH'))   CR='HR'; end
if(strcmp(upper(res),'MEDIUM')) CR='MR'; end
if(strcmp(upper(res),'LOW'))    CR='LR'; end
if(src(1) == 1) IR = '';  end
if(src(1) == 2) IR = '2'; end
if(src(2) == 1) SR = '';  end
if(src(2) == 2) SR = '2'; end 

phome = ['/home/chepplew/projects/sno/iasi' IR '_cris' SR '/' CR '/figs/'];

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
              '_' cyr smn '-' emn '_noaa_poff_'];

wnbnd = [floor(fc(1)-10) ceil(fc(end)+10)];

%Get a suitable window channel for each band
if(strcmp(band,'LW'))
  wvn = 900;
  ich = find(fi>wvn,1);
  cch = find(fc>wvn,1);
  dch = find(fd>wvn,1);
  wvn = num2str(wvn);
end
if(strcmp(band,'MW'))
  wvn=1231;
  ich = find(fi>wvn,1);
  cch = find(fc>wvn,1);
  dch = find(fd>wvn,1);
  wvn = num2str(wvn);
end
if(strcmp(band,'SW'))
  wvn=2360;
  ich = find(fi>wvn,1);
  cch = find(fc>wvn,1);
  dch = find(fd>wvn,1);
  wvn = num2str(wvn);
end


% ----------------------------------------------------------------
%                     PLOTTING SECTION 
% ----------------------------------------------------------------
%
figure(1);clf;plot(fi,r.ibm,'-', fc,r.cbm,'-', fd,r.dbm,'-');
  grid on;legend('IASI','CrIS','I2C','Location','southEast');
  xlabel('wavenumber cm^{-1}');ylabel('B.T. (K)');title('');
  %saveas(gcf,[phome pfnam_pref 'spectrum.fig'],'fig')
  %saveas(gcf,[phome pfnam_pref 'spectrum.png'],'png')

% ------------ Maps ----------------------
title2=['ASL I.' num2str(src(1)) ':C.' num2str(src(2)) ' SNO ' r.sdate ' to ' r.edate ' ' ...
         band ' overview'];
if(~VIS) fh2 = figure('visible','off'); end
if(VIS)  fh2 = figure(2);clf; end
  set(fh2,'Resize','Off'); set(fh2,'Position',fh2.Position+[0 0 380 210]);
  subplot(2,2,1); simplemap(s.clat, s.clon, s.tdiff*24*60);title('Delay IASI-CrIS mins');
  subplot(2,2,2); simplemap(s.clat, s.clon, s.dist); title('Separation deg');

  junk = r.cbt(cch,:)';
  subplot(2,2,3);simplemap(s.clat(s.iok), s.clon(s.iok), junk);title('CrIS BT (K)');
  subplot(2,2,4);simplemap(s.clat(s.iok), s.clon(s.iok), r.btbias(cch,:)');
    hcb = colorbar;ylabel(hcb,[wvn ' cm^{-1} dBT (K)']);
  title('SNO Bias IASI minus CrIS (K)');
  annotation('textbox', [0 0.9 1 0.1], 'String', title2, 'fontSize',16,...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
  %saveas(fh2,[phome pfnam_pref 'maps.fig'],'fig');
  %saveas(fh2,[phome pfnam_pref 'maps.png'],'png');


% ------------ Histograms -----------------
lat_pdf = histcounts(s.clat,[-90:2:90]);
figure(3);clf;plot([-89:2:89], lat_pdf,'.-')
   xlabel('latitude');ylabel('population');title('AIRS:CrIS SNO density vs latitude')
   grid on;
   %aslprint([phome pfnam_pref '_pop_vs_lat.png'])
   
pc_diff_pdf = 100.*(r.pdf_crad - r.pdf_drad)./r.pdf_crad;
title3=['ASL I.' num2str(src(1)) ':C.' num2str(src(2)) ' SNO ' r.sdate ' to ' ...
    r.edate ' ' num2str(wvn) ' cm^{-1} pdfs ' vers];
if(~VIS) fh3 = figure('visible','off'); end
if(VIS)  fh3 = figure(3);clf; end
  set(fh3,'Resize','off');set(fh3,'Position',fh3.Position+[0 0 360 180]);
  h1 = subplot(221);plot(r.btcens,r.pdf_crad(cch,:),'.-', r.btcens,r.pdf_drad(dch,:),'.-',...
    r.btcens,r.pdf_irad(ich,:),'-'); grid on;xlim([190 330]);
    xlabel('Scene BT bin (K)');ylabel('Number in bin');legend('CrIS','IASItoCrIS','IASI')
    title('')
  h2=subplot(223);plot(r.biascens, r.pdf_bias(cch,:), '.-');grid on;
    xlabel('bin BT (K)');ylabel('population');legend([wvn ' cm^{-1}']);
    title('IASI:CrIS SNO bias');
  h3=subplot(222);plot(r.btcens,r.pdf_crad(cch,:),'.-', r.btcens,r.pdf_drad(dch,:),'.-');
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
title4=['ASL I.' num2str(src(1)) ':C.' num2str(src(2)) ' SNO ' r.sdate ' to ' r.edate ...
        ' Mean Bias ' band ' ' vers ''];
if(~VIS) fh4 = figure('visible','off'); end
if(VIS)  fh4 = figure(4);clf; end
  set(fh4,'Resize','Off');set(fh4,'Position',fh4.Position+[0 0 280 210]);
  h1=subplot(221);plot(fd,r.bias_mn,'-', fd,10*r.btser,'-');
  axis([wnbnd(1) wnbnd(2) -0.8 0.8]);grid on;
  legend('CrIS - IASI','10*std.err.');
  xlabel('wavenumber cm^{-1}');ylabel('CrIS minus IASI (K)');

  h2=subplot(222);hold on; for i=1:9 plot(fc,r.fov(i).bias_mn,'-','color',cc(i,:)); end
  axis([wnbnd(1) wnbnd(2) -0.8 0.8]);grid on;
  legend('1','2','3','4','5','6','7','8','9','Location','eastOutside',...
         'orientation','vertical');
  xlabel('wavenumber cm^{-1}');ylabel('CrIS minus IASI (K)');
% with FOV 5 as the reference
  h3=subplot(223);hold on;
  for i=[1:4 6:9] plot(fc,r.fov(i).bias_mn - r.fov(5).bias_mn,'-','color',cc(i,:)); end
  grid on; axis([wnbnd(1) wnbnd(2) -0.4 0.4]);
  legend('1','2','3','4','6','7','8','9','Location','eastOutside');
  xlabel('wavenumber cm^{-1}');ylabel('dBT (K)');
  title('C minus I rel. FOV 5')
  annotation('textbox', [0 0.9 1 0.1], 'String', title4, 'fontSize',16,...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
  %saveas(fh4,[phome pfnam_pref 'mean_bias.fig'],'fig');
  %saveas(fh4,[phome pfnam_pref 'mean_bias.png'],'png');


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
title5 = ['ASL I.' num2str(src(1)) ':C.' num2str(src(2)) ' SNO ' r.sdate ...
          ' to ' r.edate ' Quantiles ' wvn ' wn Bias'];
if(~VIS) fh5 = figure('visible','off'); end
if(VIS)  fh5 = figure(5);clf; end
  set(fh5,'Resize','Off');set(fh5,'Position',fh5.Position+[0 0 240 0]);
  h1=subplot(221);plot(r.q.qn(cch,1:end-1), r.q.btbias(cch,:), '.-');grid on;
  axis([190 300 -1 1]); ylabel('IASI - CrIS K')
  h2=subplot(223);semilogy(r.q.qn(cch,1:end-1), r.q.binsz(cch,:),'-');grid on;xlim([190 300]);
  xlabel(['Scene BT at ' sprintf('%5.1f',fc(cch)) ' wn (K)']);ylabel('population');
    linkaxes([h1 h2],'x');set(h1,'xticklabel','');
    pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
    pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])

% --- quantiles subset by FOV  -------
  h1=subplot(222); hold on;
  for i=1:9
    plot(r.fov(i).qn(cch,1:end-1),r.fov(i).qbtbias(cch,:),'.-','color',cc(i,:));
  end
  grid on;axis([190 300 -1.9 1.9]);
  hleg = legend({'1','2','3','4','5','6','7','8','9'},'orientation','vertical',... %'Location','north',,...
     'Position',[0.9176 0.3373 0.0625 0.3833]);
  h2=subplot(224);hold on;
  for i=1:9 semilogy(r.fov(1).qn(cch,1:end-1), r.fov(i).qbinsz(cch,:),'.-');end
  grid on;xlim([190 300]); legend('bin population')
    linkaxes([h1 h2],'x');set(h1,'xticklabel','');
    pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
    pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])
  annotation('textbox', [0 0.9 1 0.1], 'String', title5, 'fontSize',14,...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
  %saveas(fh5,[phome pfnam_pref 'quantiles.fig'],'fig');
  %saveas(fh5,[phome pfnam_pref 'quantiles.png'],'png');

  
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

