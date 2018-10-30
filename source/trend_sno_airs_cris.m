trend_sno_airs_cris

%
%  choose CrIS channels that align with good AIRS channels: 
%    Use a2c_good_chans.m to identify an uniterrputed sequence of channels
%    e.g. 372:388 (880 to 890 cm-1)
%  
% channels of interest (see HM global map trends)
%  902.040, 902.387, 1231.330, 1613.862, 2499.533, 2500.601 cm-1
% For S.deSM trend retrievals:
%
%
%

addpath /home/chepplew/gitLib/asl_sno/source

disp(['Loading AIRS:CrIS.1 SNO for 2017, LW window channels']);
s1 = load_sno_airs_cris_asl_mat('2012/01/01','2012/12/31',[379:388], 'low',1,'v20a');
s2 = load_sno_airs_cris_asl_mat('2013/01/01','2013/12/31',[379:388], 'low',1,'v20a');
s3 = load_sno_airs_cris_asl_mat('2014/01/01','2014/12/31',[379:388], 'low',1,'v20a');
s4 = load_sno_airs_cris_asl_mat('2015/01/01','2015/12/31',[379:388], 'low',1,'v20a');
s5 = load_sno_airs_cris_asl_mat('2016/01/01','2016/12/31',[379:388], 'low',1,'v20a');
s6 = load_sno_airs_cris_asl_mat('2017/01/01','2017/12/31',[379:388], 'low',1,'v20a');
s7 = load_sno_airs_cris_asl_mat('2018/01/01','2018/06/30',[379:388], 'low',1,'v20a');

s = [s1, s2, s3, s4, s5, s6, s7];

c = struct;
d = struct;

iok   = [s(:).ig];
c.ctime = []; rc = []; rd = [];
for i=1:length(s) 
  c.ctime = [c.ctime s(i).cTime(s(i).ig)']; 
  rc      = [rc    s(i).rc(:,s(i).ig)];
  rd      = [rd    s(i).ra2c(:,s(i).ig)];
end
  whos rc rd s*; 
  
% Simple Bias from rms of channel group
c.dtime   = datetime(c.ctime,'convertfrom','datenum');
c.rad_rms = rms(rc,1);
d.rad_rms = rms(rd,1);

c.bt_rms = rad2bt(mean(s(1).fc(s(1).cchns)), c.rad_rms);
d.bt_rms = rad2bt(mean(s(1).fd(s(1).dchns)), d.rad_rms);

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

%{
% Plotting section:

figure(1);clf;h1=subplot(211);plot(c.dtime, c.bt_rms,'.', c.dtime, d.bt_rms,'.');
   ylim([190 320]); grid on;
  h2=subplot(212);plot(c.dtime, c.bt_rms - d.bt_rms,'.');
   ylim([-20 20]); grid on;
   hold on; plot(c.dtime, M, 'r-');
   plot(c.dtime, yval_1,'k-','linewidth',1.5);
   
figure(2);clf;plot(btcens,c.pdf,'.-', btcens,d.pdf,'.-');legend('CrIS','IASI');
  grid on;
figure(2);clf;plot(dbtcens, c.pdf_bias, '.-');grid on;

fh3=figure(3);clf;set(gcf,'Resize','off');set(fh3,'Position',fh3.Position+[0 0 240 0]);
  h1=subplot(211);plot(c.dtime, c.bt_rms,'.', c.dtime, d.bt_rms,'.');grid on;
    ylim([180 335]);title('AIRS:CrIS.1 SNO LW window channel trend');ylabel('BT (K)');
  h2=subplot(212);plot(c.dtime, c.bt_rms - d.bt_rms,'.');ylabel('d(BT) K');
   ylim([-23 23]); grid on;
   hold on; plot(c.dtime, M, 'r-');
  phome='/home/chepplew/projects/sno/airs_cris/LR/figs/';
  % saveas(gcf,[phome 'sno_ac1_lr_lw_window_trend_2012_18.png'],'png');  
%{
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])  
%}
  
