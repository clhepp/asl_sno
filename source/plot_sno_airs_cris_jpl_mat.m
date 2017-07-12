% plot_sno_airs_cris_jpl_mat.m
%
% Working with data supplied by: allchns_sno_airs_cris_jpl_mat.m
%  [ sa qn ] = allchns_sno_airs_cris_jpl_mat('2013/01/01', 1);

addpath /asl/matlib/plotutils              % aslprint.m

% Home for plots (set to where the paper is being written)
phome = '/home/chepplew/projects/sno/airs_cris/figs/';
phome = '/home/chepplew/projects/sno/sno_paper_2016/figs/';

%
% Prep some data:
igrp  = 1;
bands = [650, 1095; 1210, 1615; 2182 2550];

% load the frequency grids:
[xd, xc] = seq_match(sa.fd, sa.fc);                             % track indexes for both grids

nLW = 713; nMW = 324; nSW = 148;
if (igrp < 1 || igrp > 3) fprintf(1,'igrp out of range (1 to 3)\n'); exit; end
if(igrp == 1) sa.ichns = [1:nLW];          cband = 'LW'; nchns = nLW;  end
if(igrp == 2) sa.ichns = [nLW+1: nLW+nMW]; cband = 'MW'; nchns = nMW;  end
if(igrp == 3) sa.ichns = [nLW+nMW+1:1185]; cband = 'SW'; nchns = nSW;  end

% Convert rad to BT, do global stats (also used in next section)  
% remove bad data first:
idx = sa.ing;
cbt = single(rad2bt(sa.fc(xc(sa.ichns)), sa.crad(:,idx)));
dbt = single(rad2bt(sa.fd(xd(sa.ichns)), sa.drad(:,idx)));
% qnstats.cbm 
% qnstats.dbm
% qnstats.stderr

figure(1);clf;plot(qn.wn, qn.cbm,'b-',qn.wn,qn.dbm,'g-');grid on;
figure(1);clf;[hax,hl1,hl2]=plotyy(qn.wn,qn.cbm - qn.dbm, ... 
  qn.wn,qn.stderr);grid on;
  xlim(hax(:),[bands(igrp,:)]);ylim(hax(1),[-0.3 0.6]);ylim(hax(2),[0 0.02]);
  hax(1).YTick = [-0.4:0.1:0.6]; hax(2).YTick = [0:0.002:0.008];
  xlabel('wavenumber cm^{-1}');ylabel(hax(1),'Mean Bias K');ya2=ylabel(hax(2),'Std. error K');
  set(ya2, 'Units', 'Normalized', 'Position', [1.05, 0.7, 0]);
  title('2013 AIRS CrIS SNO LW Bias'); legend('bias','std. error','Location','north');
  % aslprint('../figs/AC_jplSNO_Bias_stdErr_LW.png',1);


%% --------------------------------------------------------------------- %%
% Choose which subset to use
clear day nit nor sou;
if(LDAY) clear idx; idx = intersect(sa.ing, sa.idd); disp([num2str(numel(idx))]); end
if(LNIT) clear idx; idx = intersect(sa.ing, sa.idn); disp([num2str(numel(idx))]); end
if(LNOR) clear idx; idx = intersect(sa.ing, sa.inh); disp([num2str(numel(idx))]); end
if(LSOU) clear idx; idx = intersect(sa.ing, sa.ish); disp([num2str(numel(idx))]); end

clear cbt dbt cbm dbm radstd;
cbt     = single(rad2bt(sa.fc(xc(sa.ichns)), sa.crad(:,idx)));
dbt     = single(rad2bt(sa.fd(xd(sa.ichns)), sa.drad(:,idx)));
cbm     = nanmean(cbt,2);
dbm     = nanmean(dbt,2);
radstd  = nanstd( sa.drad(:,idx) - sa.crad(:,idx), 0, 2 );
cdbm    = 0.5*( nanmean(dbt,2) + nanmean(cbt,2) );
  mdr   = 1E-3*( 1./drdbt(qn.wn,cdbm') );       % s.fc(xc(s.ichns))
btstd   = mdr.*radstd';  
stderr  = btstd'./sqrt(numel(idx));

if(LDAY) day.cbm = cbm; day.dbm = dbm; day.stderr = stderr; day.btbias = cbm-dbm; 
  day.cbt = cbt; day.dbt = dbt;  end
if(LNIT) nit.cbm = cbm; nit.dbm = dbm; nit.stderr = stderr; nit.btbias = cbm-dbm; 
  nit.cbt = cbt; nit.dbt = dbt;  end
if(LNOR) nor.cbm = cbm; nor.dbm = dbm; nor.btser = btser; nor.btbias = btbias; 
  nor.cbt = cbt; nor.dbt = dbt;  end
if(LSOU) sou.cbm = cbm; sou.dbm = dbm; sou.btser = btser; sou.btbias = btbias; 
  sou.cbt = cbt; sou.dbt = dbt;  end

%{ 
% ----------------- Geolocation -------------------------
 ich=402;
 figure(1);clf;simplemap(s.clat(idx),s.clon(idx),cbt(ich,:)); 
 hc = colorbar; ylabel(hc,'Kelvin')
 % aslprint([phome 'AC_jplSNO_900wn_map.png'])
 
% ------------------ Global Stats Spectra ----------------------
j = igrp; wn = qn.wn;
figure(3);clf;[hax,hl1,hl2]=plotyy(s.fc(xc(s.ichns)),cbm - dbm,s.fc(xc(s.ichns)), stderr);
  grid on;xlim(hax(:),[bands(igrp,:)]);ylim(hax(1),[-0.3 0.6]);ylim(hax(2),[0 0.02]);
  hax(1).YTick = [-0.4:0.1:0.6]; hax(2).YTick = [0:0.002:0.006];
  xlabel('wavenumber cm^{-1}');ylabel(hax(1),'Mean Bias K');ya2=ylabel(hax(2),'Std. error K');
  set(ya2, 'Units', 'Normalized', 'Position', [1.05, 0.5, 0]);
  title('2013 IASI CrIS SNO LW Bias'); legend('bias','std. error','Location','north');
  % aslprint([phome 'AC_jplSNO_Bias_stdErr_MW.png']);

% ------------------ Subset Plots ----------------------- %
j = igrp; wn = qn.wn;
% Choose which subset to plot, mean spectra, bias and std-error.
a = nor; b = sou; %% a = day; b = nit;
figure(1);clf;h1=subplot(2,1,1);plot(wn,a.cbm,'-',wn,a.dbm,'-',wn,b.cbm,'-',wn,b.dbm,'-');
    ylabel('BT (K)');grid on;axis([bands(igrp,:) 205 270]);
    %title('AC SNO Day,Night Mean, Bias, std error');
    legend('CrIS day','AIRS day','CrIS night','AIRS night','Location','south');
  h2=subplot(2,1,2);plot(wn,a.btbias,'-',wn, a.stderr,'-',wn,b.btbias,'-',wn, b.stderr,'-');
    xlabel('wavenumber cm^{-1}'); grid on;axis([bands(igrp,:),-0.4 0.4]);ylabel('dBT (K)');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp1=get(h1,'position');
  set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position'); set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1]);

figure(2);clf;[hax,hl1,hl2]=plotyy(wn,qn.cbm - qn.dbm, wn, qn.stderr);
  grid on;xlim(hax(:),[bands(igrp,:)]);ylim(hax(1),[-0.3 0.6]);ylim(hax(2),[0 0.02]);
  hax(1).YTick = [-0.4:0.1:0.6]; hax(2).YTick = [0:0.002:0.008];
  xlabel('wavenumber cm^{-1}');ylabel(hax(1),'Mean Bias K');ya2=ylabel(hax(2),'Std. error K');
  set(ya2, 'Units', 'Normalized', 'Position', [1.05, 0.7, 0]);
  title('2013 AIRS CrIS SNO LW Bias'); legend('bias','std. error','Location','north');
  % aslprint('../figs/AC_jplSNO_Bias_stdErr_MW.png',1);

% ---------------- Add module boundaries: ------------------- %
x = load('/home/chepplew/gitLib/asl_sno/data/L2.chan_prop.2002.10.22.v9.5.3.mat');
cmods = unique(x.cmod);                   % <- 17 cell
for i = 1:length(cmods)
  j = ismember(x.cmod,cmods{i});            % provides index to use for valid cmods{i}
  ModEds(i,:) = [max(x.cfreq(j)) min(x.cfreq(j))];
end
junk = sort(reshape(ModEds,34,1));
clear ModEds; ModEds = junk;
%ModEd = [652.16,  679.27,  712.06,  727.13,  731.91,  751.02,  790.58,  841.40, ...
%         869.87,  890.29,  921.15,  952.30,  1029.47, 1113.44, 1118.19, 1122.45, ...
%         1253.88, 1611.43, 1384.30, 1397.51, 2218.02, 2560.46, 2342.81, 2642.94];
% Get the AIRS detector module list:
%FH = fopen('/home/chepplew/gitLib/asl_sno/data/AIRS_DetectorModule_Table.csv','r');
%  for i = 1:2 hdr = fgetl(FH); end
%  AMods = textscan(FH,'%f %s %f %f %f %f %f %f %f %f','Delimiter',',','EmptyValue',NaN);
%fclose(FH);
%ModEds = [AMods{9} AMods{10}];           % <- [34 x 2]
%
% use [34:-1:17] for LW. [16:-1:5] for MW. [4:-1:1] for SW
figure(2); hold on; 
  %for i=1:length(ModEds) line([ModEds(i,1) ModEds(i,1)],[-1 1],'Color','green','Linestyle','--'); end
  if(igrp == 1)
   for i=1:1:17 line([ModEds(i,1) ModEds(i,1)],[-1 1],'Color','green','Linestyle','-'); end
  elseif(igrp == 2)
   for i=17:1:30 line([ModEds(i,1) ModEds(i,1)],[-1 1],'Color','green','Linestyle','-'); end
  elseif(igrp == 3)
   for i=4:-1:1 line([ModEds(i,1) ModEds(i,1)],[-1 1],'Color','green','Linestyle','--'); end
  end   
 % aslprint([phome 'AC_jplSNO_Bias_stdErr_wModEdges_LW.png']);

% -------------- HISTOGRAMS: for 900 wn channel:  ----------------- %
  ich = 403; ich = [360:480];
  btbins   = [180:1:320];
  c403_pdf = histcounts(cbt(ich,:),btbins);    d403_pdf = histcounts(dbt(ich,:),btbins);
  btbinx   = (btbins(1:end-1)+btbins(2:end))/2;
  figure(4);clf;plot(btbinx,c403_pdf,'+-', btbinx,d403_pdf,'o-');grid on;
    xlabel('Scene BT (K)');ylabel('population');    legend('CrIS','AIRS');
    %title('AIRS CrIS SNO 900wn channel population');
  % aslprint([phome 'AC_jplSNO_900wn_population_vs_scene.png'])
% day/night
  c404d_pdf = histcounts(day.cbt(ich,:),btbins); d404d_pdf = histcounts(day.dbt(ich,:),btbins);
  c404n_pdf = histcounts(nit.cbt(ich,:),btbins); d404n_pdf = histcounts(nit.dbt(ich,:),btbins);
  figure(2);clf;h1=subplot(2,1,1);hold on;plot(btbinx,c404d_pdf,'+-');plot(btbinx,d404d_pdf,'o-');
    grid on;legend('CrIS day','AIRS day');title('AC SNO 900wn sample pop vs scene (day night)');
    h2=subplot(2,1,2);hold on;plot(btbinx,c404n_pdf,'+-');plot(btbinx,d404n_pdf,'o-');grid on;
    xlabel('Scene BT (K)');ylabel('population');legend('CrIS night','AIRS night');
    % aslprint('../figs/AC_jplSNO_900wn_population_vs_scene_day_night.png');
% north (boreal)/south (austral) hemisphere
  c404b_pdf = histcounts(nor.cbt(ich,:),btbins); d404b_pdf = histcounts(nor.dbt(ich,:),btbins);
  c404a_pdf = histcounts(sou.cbt(ich,:),btbins); d404a_pdf = histcounts(sou.dbt(ich,:),btbins);
  figure(4);clf;h1=subplot(2,1,1);hold on;plot(btbinx,c404b_pdf,'+-');plot(btbinx,d404b_pdf,'o-');
    grid on;legend('CrIS NH','AIRS NH');title('AC SNO 900wn sample pop vs scene (north south)');
    h2=subplot(2,1,2);hold on;plot(btbinx,c404a_pdf,'+-');plot(btbinx,d404a_pdf,'o-');grid on;
    xlabel('Scene BT (K)');ylabel('population');legend('CrIS SH','AIRS SH');
    % aslprint('../figs/AC_jplSNO_900wn_population_vs_scene_NH_SH.png');

%  Assymetric binning - change test range and bins as desired.
ibinc  = find(cbt(ich,:) > 300 & cbt(ich,:) <= 302);  
ibina  = find(dbt(ich,:) > 300 & dbt(ich,:) <= 302);
cbinc_pdf = histcounts(cbt(ich,ibinc),[287:0.2:327]); bincen = [287.1:0.2:326.9];
dbinc_pdf = histcounts(dbt(ich,ibinc),[287:0.2:327]);
cbind_pdf = histcounts(cbt(ich,ibina),[287:0.2:327]);
dbind_pdf = histcounts(dbt(ich,ibina),[287:0.2:327]);
figure(2);clf;plot(bincen,cbinc_pdf,'+-',bincen,dbinc_pdf,'o-');grid on;xlim([296 316])
  legend('CrIS','AIRS');title('CrIS bins, AIRS matches');xlabel('Scene BT K');
figure(3);clf;plot(bincen,cbind_pdf,'+-',bincen,dbind_pdf,'o-');grid on;xlim([296 316])
  legend('CrIS','AIRS');title('AIRS bins, CrIS matches');xlabel('Scene BT K');
%}


% --------------- Quantiles Plotting --------------------- %
ich = find(qn.wn > 900,1);     % ich = find(qn.wn > 2.4975e+03,1)  % 900wn:401 2498:128
bincen = qn.binqa{ich};
figure(1);clf;
  h1=subplot(2,1,1);plot(bincen,qn.bias{ich},'k.-',bincen,1.42*qn.btser{ich},'r-',...
   bincen,-1.42*qn.btser{ich},'r-');grid on;
    axis([190 310 -1.2 0.3]);ylabel('d(BT) K');
    %title('2013 AC SNO 900wn channel bias');
  h2=subplot(2,1,2);semilogy(bincen,qn.binsz{ich});axis([190 310 100 50000]);grid on;
    ylabel('population');xlabel('Scene BT (K)');
figure(2);clf;subplot(2,1,1)
  plot(bincen(blo(jj,1):bhi(jj,1)),btbias(jj,blo(jj,1):bhi(jj,1)));grid on;
 subplot(2,1,2);plot(bincen(blo(jj,2):bhi(jj,2)),btbias(jj,blo(jj,2):bhi(jj,2)));grid on;

figure(3);clf;h1=subplot(2,1,1);plot(wmstats.wn,wmstats.cbm,'b-',wmstats.wn,wmstats.dbm,'g-');
  grid on;title('Airs (g) CrIS (b) 2013 6mos SNO mean BT');ylabel('BT K');
  h2=subplot(2,1,2);plot(wmstats.wn,wmstats.dbm - wmstats.cbm,'m.-');ylim([-1 1]);
  grid on;xlabel('wavenumber');ylabel('BT Bias K');legend('A2C - CrIS','Location','north');
  %%ha=findobj(gcf,'type','axes');set(ha(1),'ylim',[-1 1]);
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp1=get(h1,'position');
  set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position'); set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1]);
  %aslprint([phome 'AC_jplSNO_900wn_bias_vs_scene.png']);


%    -------------- HOT scene Investigation -------------------- %
%
ich = 402;
inhot_c = find(cbt(ich,:) > 290 & cbt(ich,:) <= 295); 
inhot_d = find(dbt(ich,:) > 290 & dbt(ich,:) <= 295);
inhot_u = union(inhot_c, inhot_d);
indhot  = sa.ing(inhot_u);
btbins  = [-9:0.1:9];
btbinx  = (btbins(1:end-1)+btbins(2:end))/2;
bias404_305t320_pdf = histcounts(cbt(ich,inhot_u) - dbt(ich,inhot_u),btbins);
addpath /asl/matlib/time
dtims_c = datetime(airs2dnum(sa.ctim),'convertfrom','datenum');

figure(1);clf;plot(btbinx,bias404_305t320_pdf,'o-');grid on;
nanmean(cbt(ich,inhot_u) - dbt(ich,inhot_u))
nanstd((cbt(ich,inhot_u) - dbt(ich,inhot_u)),0,2)
figure(2);clf;plot(dtims_c(indhot), cbt(ich,inhot_u) - dbt(ich,inhot_u),'.');grid on;
figure(3);clf;simplemap(sa.clat(indhot),sa.clon(indhot),cbt(ich,inhot_u)-dbt(ich,inhot_u));

inLat = find(sa.clat < -32 & sa.clat > -48 & sa.clon < -65 & sa.clon > -75);
inhot_Lat = intersect(inLat, inhot_u);

figure(4);simplemap(sa.clat(sa.ing(inhot_Lat)), sa.clon(sa.ing(inhot_Lat)), cbt(ich,inhot_Lat)); 
figure(5);simplemap(sa.clat(sa.ing(inhot_Lat)), sa.clon(sa.ing(inhot_Lat)), dbt(ich,inhot_Lat));
figure(6);simplemap(sa.clat(sa.ing(inhot_Lat)), sa.clon(sa.ing(inhot_Lat)), cbt(ich,inhot_Lat) - dbt(ich,inhot_Lat));
axis([-95 -35 -60 -5])

% ----------------- Bias for each CrIS IFOV -----------------
for i=1:9 cfov_ind{i} = find(sa.cifov == i); end

% Convert rad to BT, do global stats (also used in next section)  
% remove bad data first:
clear idx cbt dbt cbm dbm cbsd dbsd;
for i = 1:9
  %idx = cfov_ind{i};
  idx = intersect(inhot_u,cfov_ind{i});
  disp(['number in cifov: ' num2str(i) ' = ' num2str(numel(idx))])
  cbt = single( real(rad2bt(sa.fc(xc(sa.ichns)), sa.crad(:,idx))) );
  dbt = single( real(rad2bt(sa.fd(xd(sa.ichns)), sa.drad(:,idx))) );
  cbm(i,:) = nanmean(cbt,2);   cbsd(i,:) = nanstd(cbt,0,2);
  dbm(i,:) = nanmean(dbt,2);   dbsd(i,:) = nanstd(dbt,0,2);
end
figure(21);clf;hold on; for i=1:9 plot(sa.fc(sa.ichns) ,cbm(i,:)-dbm(i,:),'-'); end
  grid on; axis([640 1100 -1 1.5]);
  xlabel('wavenumber cm^{-1}');ylabel('d(BT) K');legend('1','2','3','4','5','6','7','8','9');
  title('2013 AC allsky SNO Mean Bias CrIS-AIRS x CrIS FOV num. 300Kbin');
  % saveas(gcf,[phome '2013_AC_SNO_Mean_Bias_by_CrIS_FOV_LW.png'],'png')
