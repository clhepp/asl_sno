function [] = get_airs_l1c_qual_flags_sno(sdate)
%
% function [] = get_airs_l1c_qual_flags_sno(sdate)
%
% INPUT: date. string format: 'yyyy/MM/dd'
%
%
% L1cProc:
% Zero means the channel was unchanged in Level-1C.;
% Bit 7 (MSB, value 128): This is a synthesized fill channel where the AIRS instrument 
%        does not have a detector;
% Bit 6: (value 64) Cleaned. See L1cCleanReason for the cause;
% Bit 5: (value 32) Shifted frequency (not used in release 6.0);
% Bit 4: (value 16) radiometric correction applied (not used in release 6.0);
% Bit 3: (value 8) unused/reserved (value 0);
% Bit 2: (value 4) unused/reserved (value 0);
% Bit 1: (value 2) unused/reserved (value 0);
% Bit 0: (LSB, value 1) Output value is a dummy/filler value because data is missing 
%        or otherwise could not be processed.
%
% L1cSynthReason:
% 0: value is preserved from Level-1B;
% 1: Filled because this channel falls in a gap between AIRS instrument modules;
% 2: Cleaned because this channel is known to be of low quality;
% 3: Cleaned because of bad (-9999.0) Level-1B radiance value;
% 4: Cleaned because of high Level-1B NeN noise measurement;
% 5: Cleaned because Level-1B reported a zero or negative value in the NeN noise measurement 
%   indicating that the channel is in too poor a state for noise level to be measured 
%   effectively;
% 6: Cleaned because the telemetry, gain, offset, or pop flag bits were set in Level-1B 
%   CalFlag (not used);
% 7: Cleaned because Level-1B radiance is unphysically hot;
% 8: Cleaned because Level-1B radiance is unphysically cold;
% 9: Cleaned because Level-1B radiance is hotter than expected based on the radiances of 
%   correlated channels;
% 10: Cleaned because Level-1B radiance is colder than expected based on the radiances of 
%    correlated channels;
% 11: Cleaned because Level-1B radiance is significantly increased by scene spatial 
%    inhomgeneity;
% 12: Cleaned because Level-1B radiance is significantly decreased by scene spatial 
%    inhomgeneity
% 100: Cleaned by runtime user command (Test mode only)
%

addpath /asl/packages/airs_decon/source               % seq_match.m

% directory to dump plots
phome = '/home/chepplew/projects/sno/airs_cris/figs/';

% get latest AIRS good channel lists, use 'kg' list
load /asl/matlib/airs/good_chans_2016

% match these to the 2645 L1c channel list, use fairs [2645], good channels are: fairs(xj)
load /home/chepplew/projects/airs/airs_f
[xi xj] = seq_match(sort(f(kg)), fairs);

% Check requested date
try 
   D1 = datenum(sdate,'yyyy/mm/dd');
catch
   error('Incorrect Date Format')
   return
end
[nyr1 nmn1 ndy1] = datevec(D1);
cyear = sprintf('%4d',nyr1);

sd = ['/home/chepplew/data/sno/airs_cris/ASL/HR/' cyear '/'];
ad = dir(strcat(sd,'sno_airs_cris_asl_wngbr_*_frmL1c.mat'));

ifn = 1;
g = load([ad(ifn).folder,'/',ad(ifn).name],'l1cProc','l1cSynthReason');

% count the occurances of modified values
tally = [];
for i = 1:length(xj)
  tally(i) = length(find(g.l1cProc(i,:) ~=0 ));
end

% Some plots
figure(1);clf;imagesc([1:1972],fairs(xj(426:500)),g.l1cProc(xj(426:500),:))
  xlabel('sample no.');ylabel('wavenumber cm^{-1}');title('2017d001 AC SNO AIRS l1cProc');
  saveas(gcf,[phome '2017d001_ac_hr_l1cproc_sample.png'],'png')  
figure(1);clf;pcolor([1:1972], fairs(xj(426:500)),g.l1cProc(xj(426:500),:) )
  shading flat 
% check alignment of good channels between the 2378 and 2645 grids
figure(2);clf;plot([1:1842],sort(f(kg(xi))),'.', [1:1842],fairs(xj),'o') 
  
