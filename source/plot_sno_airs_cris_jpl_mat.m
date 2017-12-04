function qnstats = plot_sno_airs_cris_jpl_mat(s)
%
% plot_sno_airs_cris_jpl_mat.m
%
% Produces stats with data supplied by: load_sno_airs_cris_jpl_mat.m
% and then selection of plts
%
%  first run: s = load_sno_airs_cris_jpl_mat(sdate1,sdate2, xchns);
%
% Results: 
%        qnstats: structure with fields of data derived from quantiles:
%        qnstats.cbm            CrIS mean spectral BT of the global sample.
%        qnstats.dbm            A2C  mean spectral BT of the global sample.
%        qnstats.wn             wavenumbers in group. [nv size]
%        qnstats.stderr         standard error of the mean bias of the global sample.
%        qnstats.bias{jj}       bias in quantile bin [nq x nv cells]
%        qnstats.btser{jj}      std error of mean bias in quantile bin [nq x nv cells]
%        qnstats.binqa{jj}      Quantile BT bins used [nq x nv cells] 
%        qnstats.binsz{jj}      Sample population in each quantile bin [nq x nv cells]
%        qnstats.b:             weighted mean bias from all quantile bins [713 x 3 single]
%        qnstats.mx:            Maximum BT bias found in each quantile bin [713 x 3 single]
%        qnstats.mn:            Minimum BT bias found in each quantile bin [713 x 3 single]
%        qnstats.ph:            the scene temp at which the max BT bias occurs [713 x 3 single]
%        qnstats.pl:            the scene temp at which the min BT bias occurs[713 x 3 single]
%        qnstats.blo:           low-end BT of quantile bin corresp. to qlo [713 x 3 single]
%        qnstats.bhi:           high-end BT of quantile bin corresp. to qhi.[713 x 3 single]
%        qnstats.qlo:           low-end quantile bin number [713 x 3 double]
%        qnstats.qhi:           high-end quantile bin number [713 x 3 double]
%


addpath /asl/matlib/plotutils                     % aslprint.m
addpath /home/strow/Git/breno_matlab/Math         % Math_bin.m
addpath /asl/matlab2012/aslutil/                  % drdbt.m

% Home for plots (set to where the paper is being written)
phome = '/home/chepplew/projects/sno/airs_cris/figs/';
phome = '/home/chepplew/projects/sno/sno_paper_2016/figs/';

% Initialize
idx = struct; 
abm = struct;
cbm = struct;
dbm = struct;

% Prep bands for convenience:
igrp  = 1;
bands = [650, 1095; 1210, 1615; 2182 2550];

% -----------------------------------------------------------------------------------
% Convert rad to BT, General stats used in next sections. Remove bad data first:
idx.all = s.ing;
idx.day = intersect(s.ing, s.idd); disp(['day: ' num2str(numel(idx.day))]);
idx.nit = intersect(s.ing, s.idn); disp(['night: ' num2str(numel(idx.nit))]);
idx.nor = intersect(s.ing, s.inh); disp(['north: ' num2str(numel(idx.nor))]);
idx.sou = intersect(s.ing, s.ish); disp(['south: ' num2str(numel(idx.sou))]);

% **** Choose which subset to use  ****
iidx = idx.all;   csub='all';

nz  = length(iidx);
abt = single(rad2bt(s.fa, s.ra(:,iidx)));
cbt = single(rad2bt(s.fc, s.rc(:,iidx)));
dbt = single(rad2bt(s.fa2c, s.ra2c(:,iidx)));
abm = nanmean(abt,2);
cbm = nanmean(cbt,2);
dbm = nanmean(dbt,2);

radstd   = nanstd( s.ra2c(:,iidx) - s.rc(:,iidx), 0, 2 );
cdbm     = 0.5*( nanmean(dbt,2) + nanmean(cbt,2) );
  mdr    = 1E-3*( 1./drdbt(s.fc,cdbm) );
btstd    = mdr.*radstd;  
stderr   = btstd./sqrt(nz);
  whos iidx cbt dbt crad drad cbm dbm radstd btstd stderr

%% ------------------------------------------------------------ %%

%cbm.day  = nanmean(single(rad2bt(s.fc,   s.rc(:,idx.day))),2);
%cbm.nit  = nanmean(single(rad2bt(s.fc,   s.rc(:,idx.day))),2);
%cbm.nor  = nanmean(single(rad2bt(s.fc,   s.rc(:,idx.day))),2);
%cbm.sou  = nanmean(single(rad2bt(s.fc,   s.rc(:,idx.day))),2);
%dbm.day  = nanmean(single(rad2bt(s.fa2c, s.ra2c(:,idx.day))),2);


% -------------------- histograms --------------------------- %
btbins = [180:1:340];
btcens = [btbins(1)+0.5:1:btbins(end)-0.5];
nach = length(s.fa); ncch = length(s.fa2c); cich = 402; 
for i=1:nach pdf_abt(i,:) = histcounts(abt(i,:),btbins);  end 
for i=1:ncch pdf_cbt(i,:) = histcounts(cbt(i,:),btbins);  end
for i=1:ncch pdf_dbt(i,:) = histcounts(dbt(i,:),btbins);  end

% ----------------------- quantiles ------------------------- %
disp('working on quantiles...');

% Get quantile profiler and set to prf.
%load('/home/strow/Matlab/Sno/prob_vector.mat');  yp = p(1:200:end);
load('/home/chepplew/projects/sno/prob_vector61.mat');               % p [1x61]
% Alternate profiler - softer tails than p.
junk = 0.0:0.1:10; yp = sigmf(junk,[2 4]); clear junk; 
junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
% Choose which profiler to use (goes in prf)
prf  = yp;

% create the scene bins for each channel
clear qsBins qxBins qdBins qx qd;
qxBins.B = quantile(cbt,prf,2);
qdBins.B = quantile(dbt,prf,2);
qx       = cell2mat(struct2cell(qxBins));
qd       = cell2mat(struct2cell(qdBins));
qsBins   = (qx + qd)/2.0;                        % x:AIRS & d:IASI-to-AIRS

% populate the scene bins for each channel (jj)
% NB the s.crad, s.drad have not been subset.
clear binsz btbias radstd cdbm btstd btser bias250 num_cbin num_dbin;
for jj = 1:length(s.fc)
  sbins = qsBins(jj,:);
  clear dbin dbinStd dbinN dbinInd xbin xbinStd xbinN xbinInd ubinInd;
  [dbin dbinStd dbinN dbinInd] = Math_bin(dbt(jj,:),cbt(jj,:)-dbt(jj,:),sbins); 
  [cbin cbinStd cbinN cbinInd] = Math_bin(cbt(jj,:),cbt(jj,:)-dbt(jj,:),sbins);

  % diagnostics: record the separate number samples in each bin:
  num_cbin(jj,:) = cbinN;
  num_dbin(jj,:) = dbinN;
  
  for i = 1:length(dbin)                                                       
    ubinInd{i} = {union(dbinInd{i},cbinInd{i})};                  
  end

  for i = 1:length(dbin)
    binsz(jj,i)   = cellfun('length',ubinInd{i});
    btbias(jj,i)  = nanmean( dbt(jj,cell2mat(ubinInd{i})) - cbt(jj,cell2mat(ubinInd{i})) );
    radstd(jj,i)  = nanstd( s.ra2c(jj,iidx(cell2mat(ubinInd{i}))) - ...
                            s.rc(jj,iidx(cell2mat(ubinInd{i}))) );
    cdbm(i)       = 0.5*( nanmean(dbt(jj,cell2mat(ubinInd{i}))) +  ...
                          nanmean(cbt(jj,cell2mat(ubinInd{i}))) );
         mdr      = 1E-3*( 1./drdbt(s.fa2c(jj),cdbm(i)) );
    btstd(jj,i)   = mdr.*radstd(jj,i);  
    btser(jj,i)   = btstd(jj,i)./sqrt(binsz(jj,i));
    %%bias250(jj,i) = btbias(jj,i)./drd250(jj);                 % option hard wired
  end
  jtot  = sum(binsz(jj,:));
  jmdr  = 1E-3*( 1./drdbt(s.fa2c(jj),cbm(jj)) );
  jbtse = jmdr.* nanstd(s.ra2c(jj,iidx) - s.rc(jj,iidx),1,2) / sqrt(jtot);
  fprintf(1,'.');
end
fprintf(1,'\n');

% parameter fitting section
% -------------------------
qnstats        = struct;
qnstats.cbm    = cbm;
qnstats.dbm    = dbm;
qnstats.wn     = s.fa2c;
qnstats.stderr = stderr;

for jj = 1:length(s.fa2c)
  clear junk;
  qnstats.bias{jj}  = btbias(jj,:);
  qnstats.btser{jj} = btser(jj,:);
  qnstats.binqa{jj} = qsBins(jj,1:end-1);
  qnstats.binsz{jj} = single(binsz(jj,:));

  % weighted mean & stats section
  % -----------------------------
  % 1. full range
  % ----------
  qlo(jj,1)  = 1;
  qhi(jj,1)  = 200;
  jtot = sum(binsz(jj,:));
  for i = 1:length(dbin)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  qnstats.b(jj,1)   = nansum(junk);  
  qnstats.mx(jj,1)  = nanmax(btbias(jj,:));  
  qnstats.mn(jj,1)  = nanmin(btbias(jj,:));
  qnstats.blo(jj,1) = qsBins(jj,1);
  qnstats.bhi(jj,1) = qsBins(jj,end);
  qnstats.ph(jj,1)  = qsBins(jj, find(btbias(jj,:) == nanmax(btbias(jj,:)),1) );
  qnstats.pl(jj,1)  = qsBins(jj, find(btbias(jj,:) == nanmin(btbias(jj,:)),1) );
  clear junk;
  % -------------------------------
  % 2. range set by min sample size
  % -------------------------------
  % range select by bin size > 500 samples
  inband    = find(binsz(jj,:) > 500);
  qlo(jj,2) = min(inband); 
  qhi(jj,2) = max(inband);                                      % was min(max(inband), find(qsBins(jj,:) > 297,1));
  jtot = sum(binsz(jj,qlo(jj,2):qhi(jj,2)));
  for i = qlo(jj,2):qhi(jj,2)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  qnstats.b(jj,2)   = nansum(junk);  clear junk;
  qnstats.mx(jj,2)  = nanmax(btbias(jj,qlo(jj,2):qhi(jj,2)));  
  qnstats.mn(jj,2)  = nanmin(btbias(jj,qlo(jj,2):qhi(jj,2)));
  qnstats.blo(jj,2) = qsBins(jj,qlo(jj,2));
  qnstats.bhi(jj,2) = qsBins(jj,qhi(jj,2));
  qnstats.ph(jj,2)  = qsBins(jj,find(btbias(jj,:) == nanmax(btbias(jj,qlo(jj,2):qhi(jj,2))),1) ); 
  qnstats.pl(jj,2)  = qsBins(jj,find(btbias(jj,:) == nanmin(btbias(jj,qlo(jj,2):qhi(jj,2))),1) );
  % ---------------------------------------
  % 3. range set by max std.err & bin size.
  % ---------------------------------------
  inband = find(btser(jj,:) < 0.04 & binsz(jj,:) > 500);        % highly tuned by trial n error
  if(numel(inband) < 2) fprintf(1,'ichn: %d\t Std Err too large\n',jj); continue; end
  qlo(jj,3) = min(inband);                                      % prob OK for all wns. 
  qhi(jj,3) = max(inband);
  jtot = sum(binsz(jj,qlo(jj,3):qhi(jj,3)));
  for i = qlo(jj,3):qhi(jj,3)
    junk(i) = binsz(jj,i).*btbias(jj,i)/jtot;
  end
  qnstats.b(jj,3)   = nansum(junk);  clear junk;
  qnstats.mx(jj,3)  = nanmax(btbias(jj,qlo(jj,3):qhi(jj,3)));  
  qnstats.mn(jj,3)  = nanmin(btbias(jj,qlo(jj,3):qhi(jj,3)));
  qnstats.blo(jj,3) = qsBins(jj,qlo(jj,3));
  qnstats.bhi(jj,3) = qsBins(jj,qhi(jj,3));
  qnstats.ph(jj,3)  = qsBins(jj,find(btbias(jj,:) == nanmax(btbias(jj,qlo(jj,3):qhi(jj,3))),1) );
  qnstats.pl(jj,3)  = qsBins(jj,find(btbias(jj,:) == nanmin(btbias(jj,qlo(jj,3):qhi(jj,3))),1) );
    
end
qnstats.qlo = qlo;
qnstats.qhi = qhi;

% save file
% ---------
savfn = ['AC_jplSNO_2013_qnstats_' csub '.mat'];

%{
fprintf(1,'Saving: %s\n',savfn);
save(['/home/chepplew/data/sno/airs_cris/' savfn],'s','qnstats','-v7.3');
%}
%{
% display summary of results for selected channel
% -----------------------------------------------
find(s.fc(s.ichns) > 900,1);
jj = 402;
sfnam = fieldnames(wmstats);
disp([wmstats.wn(jj)]);
% display stats for selected wavenumber. (the three methods should be similar)
for i = 9:length(sfnam)
  %disp([sfnam{i}  wmstats.(sfnam{i})(jj,:)]);
  fprintf(1,'%s    \t%8.4f\t%8.4f\t%8.4f\n', sfnam{i}, wmstats.(sfnam{i})(jj,:) );
end
%}
%{
% -------------------------------------------------------------------------------
%                     Plotting Section
% -------------------------------------------------------------------------------

figure(1);plot(s.fa,abm,'-', s.fc, cbm,'-',s.fa2c,dbm,'-');grid on;legend('AIRS','CrIS','A2C');
figure(1);plot(s.fc, cbm - dbm, '-', s.fc,100*stderr,'-');grid on;
  legend('Cris - Airs','100*std.error');axis([640 1100 -0.7 0.7]);
  xlabel('wavenumber cm^{-1}');ylabel('CrIS - AIRS (K)');
  title('2013 JPL AC SNO Bias LW');
  %aslprint([phome '2013_ac_jpl_sno_bias_stderr_lw_08112017.png'])

figure(1);clf;plot(qn.wn, qn.cbm,'b-',qn.wn,qn.dbm,'g-');grid on;
figure(1);clf;[hax,hl1,hl2]=plotyy(qn.wn,qn.cbm - qn.dbm, ... 
  qn.wn,qn.stderr);grid on;
  xlim(hax(:),[bands(igrp,:)]);ylim(hax(1),[-0.3 0.6]);ylim(hax(2),[0 0.02]);
  hax(1).YTick = [-0.4:0.1:0.6]; hax(2).YTick = [0:0.002:0.008];
  xlabel('wavenumber cm^{-1}');ylabel(hax(1),'Mean Bias K');ya2=ylabel(hax(2),'Std. error K');
  set(ya2, 'Units', 'Normalized', 'Position', [1.05, 0.7, 0]);
  title('2013 AIRS CrIS SNO LW Bias'); legend('bias','std. error','Location','north');
  % aslprint('../figs/AC_jplSNO_Bias_stdErr_LW.png',1);


 
% ----------------- Geolocation -------------------------
 ich=402;
 figure(1);clf;simplemap(s.clat(iidx),s.clon(iidx),cbt(ich,:) - dbt(ich,:)); 
   hc = colorbar; ylabel(hc,'Kelvin');title('2013 JPL AC SNO 900wn bias C-A');
 % aslprint([phome '2013_ac_jpl_sno_900wn_bias_map.png'])
 tdiff = s.atim - s.ctim; 
 figure(1);clf;simplemap(s.clat(iidx),s.clon(iidx),tdiff(s.ing));hc=colorbar;
   ylabel(hc,'Delay (secs)');title('2013 JPL AC SNO delay AIRS - CrIS');
 % aslprint([phome '2013_ac_jpl_sno_delay_map.png'])
   
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


% ----------------------------------------------------------------------- %
%                 HISTOGRAMS: for 900 wn channel:
% ----------------------------------------------------------------------- %
figure(2);clf;plot(btcens,pdf_abt(793,:),'-', btcens,pdf_cbt(402,:),'-',btcens, ...
           pdf_dbt(402,:),'-'); grid on; xlim([180 340]);
  xlabel('Scene BT (K)');ylabel('population');    legend('AIRS','CrIS','A2C');
  title('2013 JPL AC SNO 900wn channel population');
  % aslprint([phome '2013_ac_jpl_sno_900wn_population_vs_scene.png'])
% day/night
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

% ------------------------------------------------------------------- %
%                     Quantiles Plotting 
% ------------------------------------------------------------------- %
ich    = find(qnstats.wn > 900,1);     % ich = find(qn.wn > 2.4975e+03,1)  % 900wn:401 2498:128
bincen = qnstats.binqa{ich};
figure(3);clf;
  h1=subplot(2,1,1);plot(bincen,qnstats.bias{ich},'.-',bincen,1.42*qnstats.btser{ich},'m-',...
   bincen,-1.42*qnstats.btser{ich},'m-');grid on; axis([190 310 -0.8 0.3]);
   ylabel('AIRS - CrIS dBT (K)');title('2013 JPL AC SNO 900wn bias vs scene');
   legend('AIRS - CrIS','std. error');
  h2=subplot(2,1,2);semilogy(bincen,qnstats.binsz{ich});axis([190 310 100 50000]);grid on;
    ylabel('population');xlabel('Scene BT (K)');
    %aslprint([phome '2013_ac_jpl_sno_bias_900wn_vs_scene.png'])
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
%}
