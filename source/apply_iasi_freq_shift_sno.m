function apply_iasi_freq_shift_sno()

% SYNPOSIS:
%
% USEAGE:
%
% INPUTS:
%
% OUTPUT:
%
% DEPENDENCIES:
%
% NOTES:
%      Take the fixed percentage frequency shifts from do_iasi_speccal for each FOV
%   and each band of IASI, and compute new full spectrum frequency grid for each FOV.
%   Take the IASI radiance Obs on the new shifted frequency grid and use Howard's decon
%   to generate translated IASI to AIRS spectrum on the baseline AIRS spectral grid. This
%   is applied to each subset FOV obs in turn and then the full data set is re-assembled.
%
% AUTHOR:
%
% REVISION:
%

% for now work in this directory:
cd /home/chepplew/projects/sno/airs_iasi/

load('/asl/data/iremis/danz/iasi_f.mat');          % fiasi [8461 x 1]
load('/asl/data/airs/airs_freq.mat'); fa=freq; clear freq; % fa    [2378x1]

% table of fixed percent frequency shifts by band and FOV (from iasi_freqCal.org)
% rows: FOV, column: band
tabShift = {
{ 0,    'LW',    'MW',    'SW'; }
{ 1,   0.3625,  1.5239,  0.9812;}
{ 2,  -0.1782,  1.0814,  0.6456;}
{ 3,   0.4124,  1.1049,  0.7334;}
{ 4,   1.2157,  1.7141,  1.0934;}
};

% Load the relative frequency shift corresponding to 1% for all channels:
h = load('/home/chepplew/projects/sno/airs_iasi/iasi_relFreqShift_vs_channel.mat');

% The IASI bands are defined with frequency ranges as follows:
ibands = {
{'LW',  645,  1210;}
{'MW', 1210,  2000;}
{'SW', 2000,  2760;}
};
% convert the band frequencies into band channel numbers and append to ibands:
ibands = [ibands; {...
{'LW',    1,  2261;}
{'MW', 2262,  5421;}
{'SW', 5422,  8461;}
}];

% adjust the nominal IASI frequency grid 'fiasi' onto new grid 'fsi' for each FOV.
fsi   = struct('A',[],'B',[],'C',[],'D',[]);
fsNam = fieldnames(fsi);
% LW
for j = 1:2261
  fsi.A(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{2}(2)) * h.ppm(j)) );
  fsi.B(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{3}(2)) * h.ppm(j)) );
  fsi.C(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{4}(2)) * h.ppm(j)) );
  fsi.D(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{5}(2)) * h.ppm(j)) );
end
% MW
for j = 2262:5421
  fsi.A(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{2}(3)) * h.ppm(j)) );
  fsi.B(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{3}(3)) * h.ppm(j)) );
  fsi.C(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{4}(3)) * h.ppm(j)) );
  fsi.D(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{5}(3)) * h.ppm(j)) );
end
% SW
for j = 5422:8461
  fsi.A(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{2}(4)) * h.ppm(j)) );
  fsi.B(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{3}(4)) * h.ppm(j)) );
  fsi.C(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{4}(4)) * h.ppm(j)) );
  fsi.D(j) = fiasi(j) * ( 1 + (cell2mat(tabShift{5}(4)) * h.ppm(j)) );
end

% figure(2);clf;plot(fiasi,(fsi.A'-fiasi)./fiasi,'g',fiasi,(fsi.B'-fiasi)./fiasi,'b'); hold on;
%   plot(fiasi,(fsi.C'-fiasi)./fiasi,'r',fiasi,(fsi.D'-fiasi)./fiasi,'c');grid;

% Get the sno data, one file at a time.
dp     = '/asl/s1/chepplew/projects/sno/airs_iasi/JPL/'; % standard/';
snoLst = dir(strcat(dp,'sno_airs_iasi_*.mat'));
fprintf(1,'Found %d SNO files\n',numel(snoLst));

ifn = 7;                                 % from August 2007 onwards.
  g = load(strcat(dp,snoLst(ifn).name));
  
% IASI to AIRS translation:
% set paths to asl libs
addpath /asl/matlib/h4tools
addpath /asl/packages/iasi_decon

% specify an SRF tabulation file
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
  
% Include screening on IASI radiance as some spectral values are (close to) zero.
% TBD

% do the translation, each FOV at a time. arad3 holds the 'iasi_frqCal_to_airs'
clear arad3;
[m n] = size(g.ra);
arad3 = single(zeros(m,n));            % to hold the shifted IASI-> AIRS spectra.
tic
for j = 1:4
  clear indF fiaX iradX aradX afrqX;
  indF  = find(g.iifov == j);
  fiaX  = fsi.(fsNam{j})';
  iradX = g.ri(:,indF);
  [aradX, afrqX] = iasi2airs(iradX, fiaX, sfile, fa);

  % timing report
  [m,n] = size(g.ri(:,indF));
  fprintf(1, 'translated %d obs in %d seconds\n', n, toc)

% re-assemble each FOV subset
  arad3(:,indF) = aradX;

end    % of FOV loop

% convert to BT (airs; original iasi; shifted iasi; iasi2airs; freqCal_iasi2airs)
bta    = real(rad2bt(fa,g.ra));              btam    = nanmean(bta,2);
%
bti    = real(rad2bt(fiasi,g.ri));           btim    = nanmean(bti,2);
%
clear btsi;  btsi  = struct('A',[],'B',[],'C',[],'D',[]);
clear btsim; btsim = struct('A',[],'B',[],'C',[],'D',[]);
for j = 1:4
  indF = find(g.iifov == j);
  btsi.(fsNam{j})   = real( rad2bt(fsi.(fsNam{j}), real(g.ri(:,indF))) );   
  btsim.(fsNam{j})  = nanmean(btsi.(fsNam{j}),2);
  % reassemble each FOV subset
  btsiA(:,indF) = btsi.(fsNam{j});
end
  btsiAm = nanmean(btsiA,2);
%
bti2a  = real(rad2bt(fa, g.i2ra));           bti2am  = nanmean(bti2a,2);
%
bti3a  = real(rad2bt(fa, arad3));            bti3am  = nanmean(bti3a,2);

%{
  figure(1);clf;plot(fsi.A,btsim.A,'b');grid;
%}

fDiffppm = 1.e6*(fa - afrqX)./fa; 

%{
  figure(3);clf;plot(fa,fDiffppm,'m.');grid;

  figure(3);clf;semilogy(fiasi,g.ri(:,indF(1)),'b',fiaX,g.ri(:,indF(1)),'r');grid;

  figure(3);clf;plot(fa,bti2ra(:,1),'b',fa,bta(:,1),'g');grid;

  figure(3);clf;plot(fa,bti2ra(:,998),'b',afrqX,bti3ra(:,998),'g');grid;

  figure(3),clf;plot(fa,btam,'b',fiasi,btim,'g',fa,bti2am,'c',afrqX,bti3am,'k');grid;
  figure(3);clf;plot(afrqX,bti2am-bti3am,'m');grid;
     xlabel('nominal frequency (wn)');ylabel('BT (K)'); title('IASI.2Airs - IASI.freqCal.2Airs SNO mean bias');
  % saveas(gcf,'./figs/iasi2airs_iasi_vs_iasiFreqCal_bias_btSpectrum.png','png');
%}


% Examine each FOV separately
figure(3);clf;ax=[600 2800 -0.3 0.3];
for j = 1:4;
  fiDiff.(fsNam{j}) = 1.e6*(fsi.(fsNam{j})' - fiasi)./fiasi;
  indF = find(g.iifov == j);
  ir.(fsNam{j})     = g.ri(:,indF);                         % original IASI
   ib.(fsNam{j})    = real( rad2bt(fiasi,ir.(fsNam{j})) );
   ibm.(fsNam{j})   = nanmean(ib.(fsNam{j}),2);

  sb.(fsNam{j})     = real( rad2bt(fsi.(fsNam{j})',ir.(fsNam{j})) );   % shifted IASI
   sbm.(fsNam{j})   = nanmean(sb.(fsNam{j}),2);
      
  i3ra.(fsNam{j})   = arad3(:,indF);                         % shifted IASI -> AIRS
   i3ba.(fsNam{j})  = real( rad2bt(fa,i3ra.(fsNam{j})) );
   i3bam.(fsNam{j}) = nanmean(i3ba.(fsNam{j}),2);

  i2ra.(fsNam{j})   = g.i2ra(:,indF);                        % original IASI -> AIRS
   i2ba.(fsNam{j})  = real( rad2bt(fa,i2ra.(fsNam{j})) );
   i2bam.(fsNam{j}) = nanmean(i2ba.(fsNam{j}),2);

  intrpsbm.(fsNam{j}) = interp1(fsi.(fsNam{j}),sbm.(fsNam{j}),fiasi,'spline','extrap');
    
  subplot(2,2,j);plot(fiasi,ibm.(fsNam{j})-intrpsbm.(fsNam{j}),'m');grid;axis(ax);
    title(sprintf('FOV %d',j));
end
% saveas(gcf,'./figs/iasi_freqCal_FOVs_vs_orig_meanBTSpectrum.png','png');
  
  figure(3);clf;plot(fiasi,fiDiff,'k');grid;
  figure(3);clf;plot(fiasi,ibm.D,'b',fsi.D,sbm.D,'g');grid;
  
