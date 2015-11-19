function spectrum_tsno_airs_cris_jpl(sd,ans)

% Computes spectral average and standard deviation from all samples supplied
%   from read_tsno_airs_cris_jpl_nc().
%
% Synposis: takes existing structure returned from 'read_tsno_airs_cris_jpl_nc'
%     and plot request affirmation string.
%
% Dependencies: rad2bt.m, drdbt.m, Nominal frequency grids for Airs, CrIS and
%    deconvolved Airs->CrIS.
%
% Author: C. L. Hepplewhite. UMBC/JCET
%
% Version: Initial 02-May-2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General initialize
addpath /home/chepplew/myLib/matlib/aslutil/           % rad2bt
addpath /asl/matlab2012/aslutil/                       % drdbt.m

% check data structure exists
if(~exist('sd','var')) fprintf('Error: variable structure absent\n'); exit; end

% determine of want to create and save plots
if(lower(ans(1)) == 'y') PLOT=1; 
elseif(lower(ans(1)) == 'n') PLOT=0; 
else fprintf('Error, plz enter y(yes) or n(no) for plotting option\n'); exit; end

%{
load('~strow/Matlab/Airs/airs_f'); fa=f;
load('~strow/Matlab/Cris/Data/freq_cris');             % 1305 chns (no guard chans)
load('~chepplew/cris_f1317.mat');  fc = cfchan;        % 1317 chns (12 guard chans)
%}
xx=load('/asl/data/airs/airs_freq.mat'); fa=xx.freq; clear xx;
xx=load('cris_freq_nogrd.mat'); f_cris=xx.vchan; clear xx;  % 1305 chns (no guard chans)
xx=load('cris_freq_2grd.mat');  fc = xx.vchan; clear xx;    % 1317 chns (12 guard chans)
fd = f_cris([1:1037 1158:1305]);                        % prior knowledge from Howards decon routine

ratpm=0; rctpm=0; rdtpm=0; ratps=0; rctps=0; rdtps=0;

for i = 1:numel(sd.nSam)
  ratpm = ratpm + sd.avra(:,i).*sd.nSam(i);    rctpm = rctpm + sd.avrc(:,i).*sd.nSam(i);
  rdtpm = rdtpm + sd.avrd(:,i).*sd.nSam(i);
  ratps = ratps + ( sd.sdra(:,i).*sd.sdra(:,i) + sd.avra(:,i).*sd.avra(:,i) )*sd.nSam(i);
  rctps = rctps + ( sd.sdrc(:,i).*sd.sdrc(:,i) + sd.avrc(:,i).*sd.avrc(:,i) )*sd.nSam(i);
  rdtps = rdtps + ( sd.sdrd(:,i).*sd.sdrd(:,i) + sd.avrd(:,i).*sd.avrd(:,i) )*sd.nSam(i);
end
gatrm = ratpm/sum(sd.nSam);   gctrm = rctpm/sum(sd.nSam);   gdtrm = rdtpm/sum(sd.nSam);

garsd = real(sqrt( ratps/sum(sd.nSam) - gatrm.*gatrm ));
gcrsd = real(sqrt( rctps/sum(sd.nSam) - gctrm.*gctrm ));
gdrsd = real(sqrt( rdtps/sum(sd.nSam) - gdtrm.*gdtrm ));
garse = garsd/sqrt(sum(sd.nSam));   gcrse = gcrsd/sqrt(sum(sd.nSam));
gdrse = gdrsd/sqrt(sum(sd.nSam));

gatbm = real(rad2bt(fairs,gatrm)); gctbm = real(rad2bt(fc,gctrm)); gdtbm = real(rad2bt(fd,gdtrm));

incd  = find(ismember(fc, fd));
biasg = gdtbm - gctbm(incd);

btm   = 0.5*(gdtbm + gctbm(incd));
mdr   = 1E-3*(1./drdbt(fd,btm) );
drse  = sqrt((gdrsd.^2 + gcrsd(incd).^2))/sqrt(sum(sd.nSam));
dbse  = mdr.*drse;

whos dbse biasg gatbm gctbm gdtbm

if(PLOT)
figure(1);clf;semilogy(fairs,gatrm,'b',fairs,garsd,'b--',fairs,garse,'b:'); grid;hold on;
    semilogy(fc,gctrm,'m',fc,gcrsd,'m--',fc,gcrse,'m:');
    semilogy(fd,gdtrm,'c',fd,gdrsd,'c--',fd,gdrse,'c:');axis([500 2700 10^-5 10^3]);
    xlabel('wn (cm-1)');ylabel('radiance');title('Airs CrIS SNO All samples Spectrum');

ax1=([600 2700 -Inf Inf]); ax2=([600 2700 -1 1]);
figure(2);clf;h1=subplot(2,1,1);plot(fairs,gatbm,'b',fc,gctbm,'g',fd,gdtbm,'r');grid;axis(ax1);
  ylabel('BT (K)');title('TSNO 2014Apr AIRS l1b (b) CRIS (g) AtoC (r)'); 
  h2=subplot(2,1,2);plot(fd,biasg,'k-',fd,dbse,'m-',fd,-1*dbse,'m-');grid;axis(ax2);
  ylabel('BT (K)');xlabel('wn (cm-1)');
    linkaxes([h1 h2],'x');set(h1,'xticklabel','');p=get(h1,'position');
    set(h1,'position',[p(1) p(2)-p(4)*0.1 p(3) p(4)*1.1])
    p=get(h2,'position'); set(h2,'position',[p(1) p(2) p(3) p(4)*1.1]);
    saveas(gcf,'spectrum_tsno_airs_cris_jpl.pdf','pdf');
    saveas(gcf,'spectrum_tsno_airs_cris_jpl.fig','fig');
end

end
