function [] = trend_sno_airs_cris()

% INPUTS: sdate  string: 'dd/mm/yyyy' start date
%         edate  string: 'dd/mm/yyyy/ end date
%
%
%  choose CrIS channels that align with good AIRS channels: 
%    Use a2c_good_chans.m to identify an uniterrputed sequence of channels
%    e.g. 372:388 (880 to 890 cm-1)
%  
% channels of interest (see HM global map trends)
%  902.040, 902.387, 1231.330, 1613.862, 2499.533, 2500.601 cm-1
% For S.deSM trend retrievals:
%
% CLH: Aug 2019. Use version of load_sno_airs_cris_asl_mat.m without L1C status flags.
%

addpath /asl/matlib/aslutil
addpath /asl/matlib/plotutils
addpath /home/chepplew/gitLib/asl_sno/source

% Get the slurm array job and assign current job.
sindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));   % 0-19999
%sindx=0

thisJob = (sindex+1);
  disp(['this job: ' num2str(thisJob)]);

igrp = thisJob


% Use Job to select wavenumbers
all_wns = [650:1095 1210:1605 2183:2550];

% Limit max group
if( igrp > length(all_wns)/5 ); 
  error('channel selection out of range'); return; 
end

% Extract wavenumbers to load:
ixwns = [(igrp - 1)*5 + 1 : igrp*5];
xwns  = all_wns(ixwns);

% Complete initialization
phome = '/home/chepplew/projects/sno/airs_cris/trending/figs/';

% load cris FOV colors
addpath /home/chepplew/gitLib/asl_sno/source/    % fovcolors.m'
cc = fovcolors;

z = load('/home/chepplew/myLib/data/cris_freq_2grd.mat');    % NSR

% load AIRS good channel list for reference (nig [1x1535] and the AIRS channel list
load('/home/chepplew/myLib/data/master_nig_01_2009.mat');
xx = load('/home/chepplew/myLib/data/airs_f2378.mat');
%plot([1:2378], xx.afchan,'o'); hold on; plot(nig, xx.afchan(nig),'r.');

%{
% wavnumbers to load
xwns  = [2540:2550];
xwns  = [1489:1496];     % MW water
xwns  = [1228:1234];     % MW window
xwns  = [896:903];       % LW window
xwns  = [697:704];       % LW CO2
xchns = [379:388];       % 885.000 to 890.6250 cm-1
xchns = [1200:1209];     %
xchns = [find(z.vchan > xwns(1),1) : find(z.vchan > xwns(end),1)];

% CHeck input parameters
% check dates
try
   D1 = datenum(sdate,'yyyy/mm/dd');
   D2 = datenum(edate,'yyyy/mm/dd');
catch
   error('Incorrect Date Format')
   return
end
[nyr1 nmn1 ndy1] = datevec(D1);
[nyr2 nmn2 ndy2] = datevec(D2);
cyr1   = sdate(1:4);     cmn1 = sdate(6:7);     cdy1 = sdate(9:10);
cyr2   = edate(1:4);     cmn2 = edate(6:7);     cdy2 = edate(9:10);

DX = D1;

st = [];
% Time window/interval to average (16 days)
for i = 1:ceil((D2-D1)/16)
  fprintf(1,'%d',i)
  DY = DX+16;
  if(DY >= D2); error('Reached End of Period'); return; end

  junk1 = datevec(DX);
  junk2 = datevec(DY);
  dstr1 = [sprintf('%4d',junk1(1)) '/' sprintf('%02d',junk1(2)) '/' ...
           sprintf('%02d',junk1(3))];
  dstr2 = [sprintf('%4d',junk2(1)) '/' sprintf('%02d',junk2(2)) '/' ...
           sprintf('%02d',junk2(3))];
  DX = DY;

  sx = load_sno_airs_cris_asl_mat(dstr1,dstr2, xchns, 'low',1,'v20a');

  st = [st sx];

end
%}

disp(['Loading AIRS:CrIS.1 SNO selected wavenumbers channels']);
s1 = load_sno_airs_cris_asl_mat('2012/05/04','2012/12/31',xwns, 'low',1,'v20a');
s2 = load_sno_airs_cris_asl_mat('2013/01/01','2013/12/31',xwns, 'low',1,'v20a');
s3 = load_sno_airs_cris_asl_mat('2014/01/01','2014/12/02',xwns, 'low',1,'v20a');
s4 = load_sno_airs_cris_asl_mat('2015/01/01','2015/12/31',xwns, 'low',1,'v20a');
s5 = load_sno_airs_cris_asl_mat('2016/01/01','2016/12/31',xwns, 'low',1,'v20a');
s6 = load_sno_airs_cris_asl_mat('2017/01/01','2017/12/31',xwns, 'low',1,'v20a');
s7 = load_sno_airs_cris_asl_mat('2018/01/01','2018/12/31',xwns, 'low',1,'v20a');
s8 = load_sno_airs_cris_asl_mat('2019/01/01','2019/12/31',xwns, 'low',1,'v20a');

s = [s1, s2, s3, s4, s5, s6, s7, s8];

%{
subset = 'south';
r1 = stats_sno_airs_cris_asl_mat(s1,'MW', subset);
r2 = stats_sno_airs_cris_asl_mat(s2,'MW', subset);
r3 = stats_sno_airs_cris_asl_mat(s3,'MW', subset);
r4 = stats_sno_airs_cris_asl_mat(s4,'MW', subset);
r5 = stats_sno_airs_cris_asl_mat(s5,'MW', subset);
r6 = stats_sno_airs_cris_asl_mat(s6,'MW', subset);
r7 = stats_sno_airs_cris_asl_mat(s7,'MW', subset);
r8 = stats_sno_airs_cris_asl_mat(s8,'MW', subset);


bias_mn = [r1.bias_mn, r2.bias_mn, r3.bias_mn, r4.bias_mn, r5.bias_mn,r6.bias_mn, r7.bias_mn];
bias_se = [r1.bias_sd/sqrt(r1.nsam), r2.bias_sd/sqrt(r2.nsam), r3.bias_sd/sqrt(r3.nsam), ... 
           r4.bias_sd/sqrt(r4.nsam), r5.bias_sd/sqrt(r5.nsam), r6.bias_sd/sqrt(r6.nsam), ...
	   r7.bias_sd/sqrt(r7.nsam)];

figure(3);clf;plot(xscale, bias_mn,'.-')
   hold on; plot(xscale, bias_se,'.-');
   xlabel('year of mean');ylabel('bias and std.err'); grid on;
   legend(num2str(xwns)); title(['AIRS:CrIS.1 Bias from ' subset ' SNOs'])
%}

% SCreen out bad data (as determined dueing the loading)
ss.clat = []; ss.clon = []; ss.ctim = []; ss.ra2c = []; ss.rc = [];
ss.asc  = []; ss.td   = []; ss.dist = []; ss.cfov = [];
for i=1:length(s) 
  ss.clat = [ss.clat; s(i).cLat(s(i).ig)];
  ss.clon = [ss.clon; s(i).cLon(s(i).ig)];
  ss.ctim = [ss.ctim; s(i).cTime(s(i).ig)];
  ss.ra2c = [ss.ra2c, s(i).ra2c(:,s(i).ig)];
  ss.rc   = [ss.rc,   s(i).rc(:,s(i).ig)];
  ss.asc  = [ss.asc;  single(s(i).aAsc(s(i).ig))];
  ss.td   = [ss.td;   s(i).tdiff(s(i).ig)];
  ss.dist = [ss.dist; s(i).dist(s(i).ig)];
  ss.cfov = [ss.cfov; single(s(i).cFov(s(i).ig))]; 
end

wnums = s1.fc;
savdr = '/home/chepplew/data/sno/airs_cris/trending/';
savfn = [savdr 'ac_sno_' sprintf('%d_%d',xwns(1),xwns(end)) 'cm_2012_2019.mat'];
save(savfn, 'ss', 'wnums','-v7.3');

% ------------------------------------
%{
c1 = struct;
d1 = struct;

iok   = [s2(:).ig];
c1.ctime = []; c1.r = []; d1.r = [];
for i=1:length(s) 
  c2.ctime  = [c2.ctime s3.cTime(s3.ig)']; 
  c2.r      = [c2.r     s3.rc(:,s3.ig)];
  d2.r      = [c2.r     s3.ra2c(:,s3.ig)];
end
  whos rc rd s*; 
  
% Simple Bias from rms of channel group
c.dtime   = datetime(c.ctime,'convertfrom','datenum');
c.rad_rms = rms(rc,1);
d.rad_rms = rms(rd,1);

c1.bt_rms = rad2bt(mean(s(1).fc(s(1).cchns)), c.rad_rms);
d1.bt_rms = rad2bt(mean(s(1).fd(s(1).dchns)), d.rad_rms);

% calc PDFs
btbins  = [190:2:310];   btcens  = [191:2:309];
dbtbins = [-20:0.2:20];  dbtcens = [-19.9:0.2:19.9];
c.pdf  = histcounts(c.bt_rms, btbins);
d.pdf  = histcounts(d.bt_rms, btbins);
c.pdf_bias = histcounts(c.bt_rms - d.bt_rms, dbtbins);

% rolling mean bias
M = movmean(c.bt_rms - d.bt_rms,100);

% trend line (use datenum not datetime and center on mid-point)
junk   = (c.ctime(1));
c.rtim = c.ctime - junk;
%b1 = c.rtim\(c.bt_rms - d.bt_rms);     % not suitabel for large vectors
% polyfit returns coefficients in descending order of powers (last is intercept)
[p1,S1,mu1] = polyfit(c.rtim,(c.bt_rms - d.bt_rms),1);
[p2,S2,mu2] = polyfit(c.rtim,(c.bt_rms - d.bt_rms),2);

yval_1 = polyval(p1,c.rtim,[],mu1);
yval_2 = polyval(p2,c.rtim,[],mu2);

% Save data
savfn = '/home/chepplew/data/sno/airs_cris/stats/sno_ac1_lr_lw_window_trend.mat';
savvn = {'s','c','d','M','p*','S*','mu*','yval*'};
disp(['Saving data to: ' savfn]);
save(savfn,savvn{:},'-v7.3');

%}
% ----------------- LLS smoothing method ---------------
%{
%load ac_sno_mw_window4_2012_2018.mat
%load /home/chepplew/data/sno/airs_cris/ac_sno_mw_window4_2012_2019.mat

% Subset here
clear iifov
for i=1:9
  iifov{i} = find(ss.cfov == i); 
end
% Subset by Latitude
iie = find(ss.clat < 40 & ss.clat > -40);
iie = find(ss.clat < -40);

% select which subset to apply
isub = ':';           % no subset
ifv = 9;
%isub = iifov{ifv};
%isub = iie;

% Select which channel of collection to use
nchan = length(wnums);

for ich = 1:nchan

  % Get selected channel
  rc   = real(ss.rc(ich,isub));
  ra2c = real(ss.ra2c(ich,isub));

  % There is a stray zero in one of these
  k1 = find(rc > 0 & ra2c > 0);

  % datetime format
  amtime = datetime(ss.ctim(isub(k1)),'convertfrom','datenum');

  % Convert to BT (OK to use generic nu?)
  bta2c = rad2bt(wnums(ich),ra2c(k1));
  btc   = rad2bt(wnums(ich),rc(k1));

  btd   = btc - bta2c;

  % Find the cuts where there is no data for 2-3 days
  a = find(hours(diff(amtime)) > 49);

  clear aa astd atime;
  for i=1:(length(a)-1)
     aa(i)    = nanmean(btd([a(i):a(i+1)]));
     astd(i)  = nanstd(btd([a(i):a(i+1)]));
     atime(i) = mean(amtime([a(i):a(i+1)]));
  end

  % Get rid of bad days (empirical)
  k2 = find(astd < 2.5);

  % Mean time between averages
  mdd = nanmean(days(diff(atime(k2))));

  % Number of mdd's per year
  yearn = round(365/mdd);
  monthn = round(yearn/12);

  y  = aa(k2);
  ys1 = smooth(y,monthn,'loess');

  plot(atime(k2),ys1,'linewidth',1,'color','b')
  hold on;grid

  ys2 = smooth(y,yearn,'loess');
  plot(atime(k2),ys2,'linewidth',2,'color',cc(ifv,:))

  ylabel('CrIS - A2Cris in K')
  hl = legend('1-month Loess Filter','1-year Loess Filter');
  title(['AIRS:CrIS.1 bias from SNO ' sprintf('%7.3f',wnums(ich)) ' cm-1'])
%}
end
%{
  aslprint([phome 'ac_sno_bias_2012_2019_698p750cm_no_subset.png'])
  aslprint([phome 'ac_sno_bias_2012_2019_698p750cm_cfov1_3.png']) 

%}
%{
% Regression for slope
  dtime = datenum(atime(k));
  xtime = [ones(length(dtime),1) dtime'-dtime(1)];

  z = double(aa(k))';
  [b bint r rint stats] = regress(z,xtime);

  lag = xcorr(real(z),1,'coeff');
  lag = lag(1);
  error_mult = sqrt((1+lag)/(1-lag));

  slope = b(2)*365;  % K/year
  slope_unc = error_mult*(diff(bint(2,:)*365)/2);

  [slope slope_unc]

%}
% ---------------------- END LLS smooth ----------------
%{
% Select a three month window to average:
datn1 = datenum('2019/09/01','yyyy/mm/dd');
datn2 = datenum('2019/12/31','yyyy/mm/dd');
iit1  = find(ss.ctim > datn1 & ss.ctim < datn2);

% Convert to BT (OK to use generic nu?)
bta2c = rad2bt(wnums,real(ss.ra2c));
btc   = rad2bt(wnums,real(ss.rc));

btd   = btc(:,iit1) - bta2c(:,iit1);
btdmn = nanmean(btd,2);
btdsd = nanstd(btd,0,2);






ik    = find(abs(btd - btdmn) < 9*btdsd);

yz = btd(ik);
ysm2 = smooth(yz,2900,'loess');
B=smoothdata(ysm2,10000);
windowSize = 50000;
 b = (1/windowSize)*ones(1,windowSize);
 a = 1;
 B = filter(b,a,ysm2);
%}
%{
% Plotting section:

figure(7);clf;h1=subplot(211);plot(c.dtime, c.bt_rms,'.', c.dtime, d.bt_rms,'.');
   ylim([190 320]); grid on;
  h2=subplot(212);plot(c.dtime, c.bt_rms - d.bt_rms,'.');
   ylim([-20 20]); grid on;
   hold on; plot(c.dtime, M, 'r-');
   plot(c.dtime, yval_1,'k-','linewidth',1.5);
   
figure(8);clf;plot(btcens,c.pdf,'.-', btcens,d.pdf,'.-');legend('CrIS','IASI');
  grid on;
figure(8);clf;plot(dbtcens, c.pdf_bias, '.-');grid on;

fh8=figure(8);clf;set(gcf,'Resize','off');set(fh8,'Position',fh8.Position+[0 0 240 0]);
  h1=subplot(211);plot(c.dtime, c.bt_rms,'.', c.dtime, d.bt_rms,'.');grid on;
    ylim([180 335]);title('AIRS:CrIS.1 SNO LW window channel trend');ylabel('BT (K)');
  h2=subplot(212);plot(c.dtime, c.bt_rms - d.bt_rms,'.');ylabel('d(BT) K');
   ylim([-23 23]); grid on;
   hold on; plot(c.dtime, M, 'r-');
  phome='/home/chepplew/projects/sno/airs_cris/LR/figs/';
  % saveas(gcf,[phome 'sno_ac1_lr_lw_window_trend_2012_18.png'],'png');  
%}
%{
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])  
%}
  
