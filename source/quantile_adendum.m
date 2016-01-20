% quantile_addendum.m

%
% multiple channel quantile/scene dependent summary: IASI.CrIS SNO channels 250:410 window
%      exhibit a very similar bias trend w/scene brightness temperature, so make a summary
%
% start with wmstats loaded for the LW band using 
%  roses2015_figs.m section 1. IASI CRIS Bias.
%
%
% take LW channels 
bias120 = nanmean(xbias(:,250:410),2);
% take window channels 250:410 (160 LW window channels)
bias160 = nanmean(xbias(:,250:410),2);
stdv160 = nanstd(xbias(:,250:410),1,2);

 figure(2);clf;plot(xbinqa(:,403),bias160,xbinqa(:,403),stdv160,'m',xbinqa(:,403),-stdv160,'m');
   grid on;axis([200 290 -1 1]); 
   xlabel('Scene BT (K)');ylabel('Bias (K)');title('IASI - CrIS 800 to 900 cm^{-1} Bias vs Scene')
   legend('Mean Bias','Std Deviation','Location','South');
   % aslprint('./figs/IC_jplSNO_bias_std_LWchns_vsScene.png');
