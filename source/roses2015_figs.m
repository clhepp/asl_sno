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
yfn = {'sno_AI_wmstats_2011x_chns400x01.mat','sno_AI_wmstats_2011x_chns400x02.mat',...
       'sno_AI_wmstats_2011x_chns400x03.mat','sno_AI_wmstats_2011x_chns400x04.mat',...
       'sno_AI_wmstats_2011x_chns400x05.mat','sno_AI_wmstats_2011x_chns400x06.mat'};

yabm = []; ydbm = []; ywn = []; ybias = [];
       
for i = 1:6

  g = load(strcat(dp,yfn{i}));
  [ym yn] = size(g.wmstats.wn);
  yabm  = [yabm;  g.wmstats.abm];
  ydbm  = [ydbm;  g.wmstats.dbm];
  ywn   = [ywn;   g.wmstats.wn];
   junk = g.wmstats.bias;
   junk(cellfun(@ischar,junk)) = {NaN};
   out  = cell2mat(junk);
   out  = reshape(out, 200, ym);  
  ybias = [ybias, out];
end       
whos yabm ydbm ywn ybias

load('/home/chepplew/projects/airs/master_nig_01_2009.mat');   % nig [1 x 1535]


%{
% sanity check
figure(11);clf;plot(ywn(nig), yabm(nig) - ydbm(nig),'m-');grid on;axis([645 1100 -0.8 0.8]);
   title('AIRS - IASI');
%}

% -------------------------------
%      3.  AIRS CrIS Bias.
% -------------------------------

zfn = {'sno_AC_wmstats_2013_chns400x01.mat','sno_AC_wmstats_2013_chns400x02.mat',...
       'sno_AC_wmstats_2013_chns400x03.mat'};

zcbm = []; zdbm = []; zwn = []; zbias = []; zbinbt = [];
       
for i = 1:3

  g = load(strcat(dp,zfn{i}));
  [zm zn] = size(g.wmstats.wn);
  zcbm  = [zcbm;  g.wmstats.cbm];
  zdbm  = [zdbm;  g.wmstats.dbm];
  zwn   = [zwn;   g.wmstats.wn'];
   junk = g.wmstats.bias;
   junk(cellfun(@ischar,junk)) = {NaN};
   out  = cell2mat(junk);
   out  = reshape(out, 200, zn);  
  zbias = [zbias, out];   clear junk out;
   junk = g.wmstats.binbt;
   junk(cellfun(@ischar,junk)) = {NaN};
   out  = cell2mat(junk);
   out  = reshape(out, 200, zn);
  zbinbt = [zbinbt, out];   clear junk out;  
end       
whos zcbm zdbm zwn zbias zbins

%{
% sanity check
figure(12);clf;plot(zwn, zdbm - zcbm,'.-');grid on;axis([bands(1,:) -1 01]);
figure(12);clf;plot(zwn(nig), zcbm(nig) - zdbm(nig),'m-');grid on;title('AIRS - CrIS');
  axis([645 1100 -0.8 0.8]);
figure(13);clf;plot(zwn(fd_notgap),zcbm(fd_notgap) - zdbm(fd_notgap), 'c');grid on;title('AIRS - CrIS');
  axis([645 1100 -0.8 0.8]);
%}

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
  
  
