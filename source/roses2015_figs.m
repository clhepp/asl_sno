% Summary plots ROSES December 2015

cd /home/chepplew/gitLib/asl_sno/run
addpath /home/chepplew/gitLib/asl_sno/source
addpath /asl/matlib/plotutils                                               % aslprint.m
addpath /asl/packages/ccast/source                                          % seq_match.m
addpath /asl/matlab2012/aslutil                                             % wmean

xx = load('/home/chepplew/projects/cris/cris_freq_2grd.mat');  fcris = xx.vchan; clear xx;    % 1317 chns (12 guard chans)
xx = load('/asl/data/airs/airs_freq.mat'); fairs=xx.freq; clear xx;                           % 2378
run '/home/chepplew/gitLib/asl_sno/data/airs_filled_gaps.m';                            % fa_gap [10 x 2]
[xa xc] = seq_match(sort(fairs),fcris);
fd    = fcris([1:1037 1158:1305]);                            % prior knowledge from Howards decon routine
[xd, xc] = seq_match(fd, fcris);                               % track indexes for both grids

bands = [640, 1100; 1200, 1620; 2170, 2400];

% prep the airs->cris AIRS gaps so they can be ommitted when plotting (LW and MW only)
fd_notgap = find(fcris > fcris(1) & fcris < fa_gap(1,1) -0.5);
for i = 1:8
  fd_notgap = [fd_notgap, find(fcris > fa_gap(i,2)+0.5 & fcris < fa_gap(i+1,1)-0.5 )];
end
 fd_notgap = [fd_notgap, find(fcris > fa_gap(2,2)+0.5 & fcris < fa_gap(3,1)-0.5 )];

 fcgd = fcris(~fd_gap);

% Get the AIRS detector module list:
FH = fopen('/home/chepplew/gitLib/asl_sno/data/AIRS_DetectorModule_Table.csv','r');
for i = 1:2 hdr = fgetl(FH); end
AMods = textscan(FH,'%f %s %f %f %f %f %f %f %f %f','Delimiter',',','EmptyValue',NaN);
fclose(FH);
ModEd = [652.16,  679.27,  712.06,  727.13,  731.91,  751.02,  790.58,  841.40, ...
         869.87,  890.29,  921.15,  952.30,  1029.47, 1113.44, 1118.19, 1122.45, ...
	 1253.88, 1611.43, 1384.30, 1397.51, 2218.02, 2560.46, 2342.81, 2642.94];
 ModEds = [AMods{9} AMods{10}];           % <- [34 x 2]
%
 x = load('../data/L2.chan_prop.2002.10.22.v9.5.3.mat');
cmods = unique(x.cmod);                   % <- 17 cell
for i = 1:length(cmods)
  j = ismember(x.cmod,cmods{i});            % provides index to use for valid cmods{i}
  ModEds(i,:) = [max(x.cfreq(j)) min(x.cfreq(j))]
end
% ----------------------------------
%      1. IASI CRIS Bias.
% -----------------------------------
dp = '/home/chepplew/gitLib/asl_sno/run/';
sfn = {'sno_IC_wmstats_2012x_chns400x01.mat','sno_IC_wmstats_2012x_chns400x02.mat',...
       'sno_IC_wmstats_2012x_chns400x03.mat'};

xcbm = []; xdbm = []; xwn = []; xbias = []; xser = []; xbinsz = []; xbinqa = [];

for i = 1:numel(sfn)
  g = load(strcat(dp,sfn{i}));
  
  [xm xn]  = size(g.wmstats.wn);
  xcbm   = [xcbm;  g.wmstats.cbm];
  xdbm   = [xdbm;  g.wmstats.dbm];
  xwn    = [xwn;   g.wmstats.wn'];
   junk  = g.wmstats.bias;
   junk(cellfun(@ischar,junk)) = {NaN};
   out   = cell2mat(junk);
   out   = reshape(out, 200, xn);          % <- 200 is the quantile profiler dimension
  xbias  = [xbias, out];
  clear junk out;
   junk  = g.wmstats.btser;
   junk(cellfun(@ischar,junk)) = {NaN};
   out   = cell2mat(junk);
   out   = reshape(out, 200, xn);          % <- 200 is the quantile profiler dimension
  xser   = [xser, out];   
  xbinsz = [xbinsz, g.wmstats.binsz];
  xbinqa = [xbinqa, g.wmstats.binqa];
end
whos xbias xcbm xdbm xwn xser xbinsz xbinqa;

[nm nn] = size(xbias);
clear aa ww xstd jtot;
xbinqa = reshape(cell2mat(xbinqa), 200, nn);
xbinsz = reshape(cell2mat(xbinsz), 200, nn);
jtot   = sum(xbinsz,1);
for i = 1:nn
  aa = squeeze(xser(:,i));   ww = squeeze(xbinsz(:,i));
  xstd(i) = wmean(aa, ww);
end

       
%{
% plotting options
figure(10);clf;plot(xwn, xdbm, 'b-');grid on;axis([bands(1,:) 200 255]);
figure(10);clf;plot(xwn, xcbm - xdbm,'.-',xwn, xstd,'.-');axis([bands(1,:) -0.45 0.45]);
  grid on;
  xlabel('Wavenumber cm^{-1}');ylabel('Bias (K)');title('2012-14 IASI CrIS SNO LW Bias'); 
  legend('Bias','Standard Error','location','South');
  % aslprint('./figs/IC_jplSNO_Bias_stErr_LW_spectrum.png');
figure(11);clf;h1=subplot(2,1,1);plot(xbinqa(:,403),xbias(:,403));hold on;grid on;
  plot(xbinqa(:,403),xser(:,403),xbinqa(:,403),-xser(:,403)); 
   axis([200 300 -1 1]);ylabel('Bias K');title('IASI CrIS Bias at 900 wn vs Scene');
   h2=subplot(2,1,2);semilogy(xbinqa(:,403),xbinsz(:,403));grid on;axis([200 300 100 10000]);
   xlabel('Scene BT (K)');ylabel('population');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])  
  % aslprint('./figs/IC_jplSNO_Bias_stErr_900wn_vsScene.png');

%}

% ----------------------------------
%       2. IASI AIRS Bias.
% ----------------------------------
dp  = '/home/chepplew/gitLib/asl_sno/run/';
yfn = {'sno_AI_wmstats_2011x_chns400x01.mat','sno_AI_wmstats_2011x_chns400x02.mat',...
       'sno_AI_wmstats_2011x_chns400x03.mat','sno_AI_wmstats_2011x_chns400x04.mat',...
       'sno_AI_wmstats_2011x_chns400x05.mat','sno_AI_wmstats_2011x_chns400x06.mat'};

yabm = []; ydbm = []; ywn = []; ybias = []; yser = []; ybinsz = []; ybinqa = [];
       
for i = 1:numel(yfn)

  g = load(strcat(dp,yfn{i}));
  [ym yn] = size(g.wmstats.wn);
  yabm  = [yabm;  g.wmstats.abm];
  ydbm  = [ydbm;  g.wmstats.i2abm];
  ywn   = [ywn;   g.wmstats.wn];
   junk = g.wmstats.bias;
   junk(cellfun(@ischar,junk)) = {NaN};
   out  = cell2mat(junk);
   out  = reshape(out, 200, ym);  
  ybias = [ybias, out];
  clear junk out;
   junk  = g.wmstats.btser;
   junk(cellfun(@ischar,junk)) = {NaN};
   out   = cell2mat(junk);
   out   = reshape(out, 200, ym);          % <- 200 is the quantile profiler dimension
  yser   = [yser, out];   
  ybinsz = [ybinsz, g.wmstats.binsz];
  ybinqa = [ybinqa, g.wmstats.binqa];
end       
whos yabm ydbm ywn ybias yser ybinsz ybinqa

[nm nn] = size(ybias);
ybinqa = reshape(cell2mat(ybinqa), 200, nn);
ybinsz = reshape(cell2mat(ybinsz), 200, nn);
clear aa ww ystd jtot;
jtot   = sum(ybinsz,1);
for i = 1:nn                              % some yser are filled with NaNs so wmean would fail
  aa = squeeze(yser(:,i));   ww = squeeze(ybinsz(:,i));
  if all(isnan(yser(:,i)))
    ystd(i) = NaN(1,1);
  else
    ystd(i) = wmean(aa, ww);
  end
end

%{
% plotting options
figure(10);clf;plot(ywn,ydbm,'-',ywn,yabm,'-');grid on;axis([bands(1,:) 200 255]);
figure(10);clf;plot(ywn, yabm - ydbm,'.-',ywn, ystd,'m-',ywn,-ystd,'m-');
  axis([bands(1,:) -0.75 0.75]);grid on;
  xlabel('Wavenumber cm^{-1}');ylabel('Bias (K)');title('2011-14 AIRS IASI SNO LW Bias'); 
  legend('Bias','Standard Error','location','South');
  % aslprint('./figs/IA_jplSNO_Bias_stErr_LW_spectrum.png');
%}
load('/home/chepplew/projects/airs/master_nig_01_2009.mat');   % nig [1 x 1535]


%{
% sanity check
figure(2);clf;plot(ywn(nig), yabm(nig) - ydbm(nig),'m-');grid on;axis([645 1100 -0.8 0.8]);
   title('AIRS - IASI');
%}

% -------------------------------
%      3.  AIRS CrIS Bias.
% -------------------------------
dp  = '/home/chepplew/gitLib/asl_sno/run/';
zfn = {'sno_AC_wmstats_2013_chns400x01_v2.mat','sno_AC_wmstats_2013_chns400x02_v2.mat',...
       'sno_AC_wmstats_2013_chns400x03_v2.mat'};

zcbm   = []; zdbm = []; zwn = []; zbias = []; zbinbt = []; zser = []; zbinsz = []; 
zbinqa = [];
       
for i = 1:numel(zfn)
  g = load(strcat(dp,zfn{i}));
  [zm zn] = size(g.wmstats.wn);
  zcbm  = [zcbm;  g.wmstats.cbm];
  zdbm  = [zdbm;  g.wmstats.dbm];
  zwn   = [zwn;   g.wmstats.wn'];
   junk = g.wmstats.bias;
   junk(cellfun(@ischar,junk)) = {NaN};
   out  = cell2mat(junk);
   out  = reshape(out, 200, zn);  
  zbias = [zbias, out];     clear junk out;
   junk = g.wmstats.binqa;
   junk(cellfun(@ischar,junk)) = {NaN};
   out  = cell2mat(junk);
   out  = reshape(out, 200, zn);
  zbinqa = [zbinqa, out];   clear junk out;  
   junk  = g.wmstats.btser;
   junk(cellfun(@ischar,junk)) = {NaN};
   out   = cell2mat(junk);
   out   = reshape(out, 200, zn);          % <- 200 is the quantile profiler dimension
  zser   = [zser, out];     clear junk out;
  zbinsz = [zbinsz, g.wmstats.binsz];
end       
whos zcbm zdbm zwn zbias zbinbt binsz zbinqa zser 

[nm nn] = size(zbias);
clear aa ww zstd jtot;
zbinqa = reshape(cell2mat(zbinqa), 200, nn);
zbinsz = reshape(cell2mat(zbinsz), 200, nn);
jtot   = sum(zbinsz,1);
for i = 1:nn
  aa = squeeze(zser(:,i));   ww = squeeze(zbinsz(:,i));
  zstd(i) = wmean(aa, ww);
end


%{
% sanity check
figure(12);clf;plot(zwn, zcbm, '.-', zwn, zdbm, '.-');grid on; axis([bands(2,:) 210 270]);
figure(12);clf;plot(zwn, zdbm - zcbm,'.-',zwn,zstd,'m-',zwn,-zstd,'m-');grid on;
  axis([bands(2,:) -0.8 0.8]);
  xlabel('Wavenumber cm^{-1}');ylabel('Bias (K)');title('2013 Stnd SNO AIRS CrIS LW bias');
  legend('Bias','Std Err','Location','North');
  % aslprint('./figs/2013_AC_jplSNO_std_MW_biasSte.png');
figure(12);clf;plot(zwn(nig), zcbm(nig) - zdbm(nig),'m-');grid on;title('AIRS - CrIS');
  axis([645 1100 -0.8 0.8]);
figure(13);clf;plot(zwn(fd_notgap),zcbm(fd_notgap) - zdbm(fd_notgap), 'c');grid on;title('AIRS - CrIS');
  axis([645 1100 -0.8 0.8]);
figure(12);clf;plot(zwn, zdbm - zcbm,'.-',zwn,zstd,'m-',zwn,-zstd,'m-');grid on;
  axis([bands(2,:) -0.8 0.8]);
 hold on; for i=1:length(ModEds) line([ModEds(i,1) ModEds(i,1)],[-1 1]); end %}

% ---------------------------------
%   All together
% ---------------------------------
% 3 lines - two panes: upper pane: AIRS (A2C) - CrIS from sno_airs_cris*
%                                  AIRS (A2C) - IASI (I2C) from sno_airs_iasi*.mat
%                      lower pane: Cris       - IASI (I2C) from sno_iasi_cris*.mat

% first load IASI CrIS (x*) from 1. & AIRS CrIS (z*) from 3. above

AI = load('AI_sno_A2C_I2C_plot_data.mat');     

figure(14);clf;
  h1=subplot(2,1,1);plot(AI.fa2c,AI.r.gavba2c - AI.r.gavbi2c(xc),'.-');axis([bands(2,:) -0.8 0.8])
    hold on; plot(fa2c,zdbm - zcbm,'.-'); ylabel('bias K') 
    grid on; legend('a2c - i2c','a2c - c','location','south');
  h2=subplot(2,1,2);plot(xwn,xcbm - xdbm, '.-');axis([bands(2,:) -0.8 0.8]); grid on;
    ylabel('bias K');xlabel('wavenumber cm^{-1}');legend('cris - iasi','Location','south')
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])  
  %aslprint('./figs/AC_AI_CI_biasBT_MW_legend.png');
  
figure(15);clf;
  h2=subplot(2,1,1);plot(ywn(nig), yabm(nig) - ydbm(nig),'.-');grid on;
    legend('AIRS - IASI','Location','southEast');axis([bands(1,:) -0.8 0.8]);

  h3=subplot(2,1,2);plot(zwn(fd_notgap),zdbm(fd_notgap) - zcbm(fd_notgap), '.-');grid on;legend('AIRS - CrIS');
    axis([bands(1,:) -0.8 0.8]);
  
  linkaxes([h1 h2 h3],'x');set(h1,'xticklabel','');set(h2,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.2 pp1(3) pp1(4)*1.2])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)-pp2(4)*0.1 pp2(3) pp2(4)*1.1])  
  pp3=get(h3,'position');set(h3,'position',[pp3(1) pp3(2)-pp3(4)*0.1 pp3(3) pp3(4)*1.1])  

  saveas(gcf,'./figs/AirsSensors_MW_meanBTspectrum.png','png');
  %aslprint('./figs/AirsSensors_MW_meanBTspectrum.png')   

% --------------------------------------
%     Scene dependency
% --------------------------------------

% AIRS CrIS

figure(10);clf;plot(zbins(:,350:400),zbias(:,350:400));axis([190 330 -1 1])

[XG,YG] = meshgrid(squeeze(zbins(:,375)),fd(350:400));
crange  = [-0.2:0.02:0.2];
figure(11);clf;contourf(XG,YG,zbias(:,350:400)',crange);                   
  xlabel('Scene BT K');ylabel('wn cm-1'); colorbar;
  %saveas(gcf,'./figs/sample_LW_window_bias_scene_contour.png','png')
  
  
