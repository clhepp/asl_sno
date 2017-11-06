function [s] = load_sno_airs_iasi_asl_mat_v2(sdate1, sdate2, xchns, src)
%
% function load_sno_airs_iasi_asl_mat() loads up radiances for a selected number
%   of channels, specified by AIRS channel number, from the ASL SNO mat files 
%   and for the specified year and months. 
%   Unlike the sister function 'read_sno...'
%   this function cacluates statistics during load and subsets of CrIS FOV.
%
% Synopsis: load_sno_airs_iasi_asl_mat('date1','date2',[chan1...chan10]);
%           date1: first month as string: 'YYYY/MM/DD'
%           date2: last month as string:  'YYYY/MM/DD'
%           N.B. Only accepts the same year.
%           xchns: numeric IDs of CrIS channels to load based on NO guard channel list. 
%           (max 10).
%           eg [403 499 737 884 905 998 1021 1297] or [566:824];
%           [   1:1269] (645  - 1100) cm-1
%           [1270:2160] (1100 - 1615) cm-1
%
% Output:  Two structures of arrays. 
%             s: the SNO single fields.
%             a: whole spectrum averages and first moment.
%
%
% Notes: If the selected channel is associated with a bad AIRS channel or has been
%        modified by AIRS L1C (clean and fill) then the next nearest good channel
%        is substituted.
%        The IASI spectra are apodized.
%
% Dependencies: i) AIRS good channel list; ii) nominal AIRS and CrIS frequency grids.
%    iii) fixed path and file name syntax for SNO files.
%
% Notes: i) No QA is applied. ii) time separation of SNO pairs from file is positive-only
%    so is recomputed here.
%
% Author: C. L. Hepplewhite, UMBC/JCET
%
% Version: 04-October-2017
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /home/chepplew/gitLib/asl_sno/run

addpath /asl/packages/airs_decon/source             % hamm_app.m
addpath /asl/matlib/aslutil                         % rad2bt.m
addpath /home/chepplew/gitLib/asl_sno/source/
addpath /home/strow/Git/breno_matlab/Math           % Math_bin.m
addpath /home/chepplew/gitLib/airs_deconv/source    % seq_match.m
addpath /home/chepplew/myLib/matlib/math            % remove_6sigma 

s = struct;

junk = [-5:.05:5]; y0 = normpdf(junk,0,1); yp = cumsum(y0)./20.0; clear junk y0;
% Choose which profiler to use (goes in prf)
s.prf  = yp;

% Process and check the date strings
posYrs = [2002:2017];
posMns = [1:12];
whos sdate1; disp([sdate1 ' to ' sdate2]); fprintf('\n');
try 
   D1 = datenum(sdate1,'yyyy/mm/dd');
   D2 = datenum(sdate2,'yyyy/mm/dd');
catch
   error('Incorrect Date Format')
   return
end
s.dtime1 = datetime(D1,'convertFrom','datenum');
[nYr1 nMn1 nDy1] = datevec(D1);
[nYr2 nMn2 nDy2] = datevec(D2);
if(nYr1 ~= nYr2) error('Use same year only'); return; end
cYr1   = sdate1(1:4);     cMn1 = sdate1(6:7);     cDy1 = sdate1(9:10);
cYr2   = sdate2(1:4);     cMn1 = sdate2(6:7);     cDy1 = sdate2(9:10);

  junk = sprintf('%4d/%02d/%02d',nYr1-1,12,31);
jdy1   = datenum(sdate1)-datenum(junk);  clear junk;           % needed for data directory
  junk = sprintf('%4d/%02d/%02d',nYr2-1,12,31);
jdy2   = datenum(sdate2)-datenum(junk);  clear junk;           % needed for data directory

% Check channel numbers entered correctly
if(length(xchns) > 20 || length(xchns) < 1 ) fprintf(1,'Wrong number channels\n'); end
if(min(xchns) < 1 || max(xchns) > 1317 ) fprintf(1,'Wrong channel numbers used\n'); end

% load IASI and AIRS channels & good AIRS channels (nig) to use, & bad (nib) to avoid
load('/home/chepplew/projects/iasi/f_iasi.mat');               % f_iasi [8641 x 1]
load('/home/chepplew/projects/airs/airs_f.mat'); f2645 = fairs;  f2378 = f;
load('/home/chepplew/projects/airs/master_nig_01_2009.mat');   % nig [1 x 1535]
% nig = importdata('/home/strow/Work/Airs/good_chan_list');
junk = ismember([1:2378], nig);  nib = find(junk == 0);  clear junk;

% Get over-lapping IASI channels
[zi ichns] = seq_match(f2645(xchns),f_iasi);
   
% Screen the channel selection for AIRS bad channels and report any:
achns  = xchns;
%{
agood  = logical(zeros(1,numel(achns)));         % set all to good.
for i=1:numel(cWavs)
  achns(i)  = find(f2645  > cWavs(i)-0.1,1);
  agood(i)  = ~isempty( find( ismember(f2645(zj), f2645(achns(i))) ) ); 
end
for i=1:numel(cWavs) sWavs{i}  = sprintf('%6.2f',cWavs(i)); end
aWavs = f2645(achns);
for i = 1:numel(aWavs)
  if(agood(i)) fprintf(1,'Good AIRS l1b channel: %d, %7.3f\n', i,aWavs(i)); end
end
%}
% ************* load up SNO data ********************

dp = ['/home/chepplew/data/sno/airs_iasi/ASL/' cYr1 '/'];
% snoLst = dir(strcat(dp, 'sno_airs_crisLR_clh_*_018d600s.mat'));
snoLst = dir(strcat(dp, 'sno_airs_iasi_*_frmL1c.mat'));
fprintf(1,'Found %d total SNO files\n',numel(snoLst));

%{
% subset range by date as requested:
dstart = datenum([nyr1 nmn1 ndy1]);
dlast  = datenum([nyr2 nmn2 ndy2]);
ifn1 = 1;             % default start with first file unless later.
for i=1:numel(snoLst)
  %junk = snoLst(i).name(15:22);               % specific file name only
  junk = regexp(snoLst(i).name,'(?<=_)[\d8]+(?=_018d)','match'); 
  thisdat = datenum( [str2num(junk{1}(1:4)) str2num(junk{1}(5:6)) str2num(junk{1}(7:8))] );
  if(thisdat <= dstart) ifn1 = i; end
  if(thisdat <= dlast) ifn2 = i; end
end
fprintf(1,'Processing SNO files from: %s to %s\n',snoLst(ifn1).name(21:28), ...
        snoLst(ifn2).name(21:28));
%}
%%%%%%%%%
s.tdiff = [];    s.ra = [];    s.ri = [];      s.rd = [];  s.itime = [];  s.atime = []; 
 s.alat = [];  s.alon = []; s.dist  = [];    s.ilat = [];   s.ilon = [];  s.csolz = [];  
s.iqual = []; s.alnfr = [];  s.ifov = []; 

for ifn = 1:numel(snoLst)  % ifn1:ifn2;
  vars = whos('-file',strcat(dp,snoLst(ifn).name));
  if( ismember('iObs', {vars.name}) & ismember('aObs', {vars.name}) & ...
      ismember('i2ra', {vars.name})  )  
    load(strcat(dp, snoLst(ifn).name));
    if  (size(iObs,1) == 8461 & size(aObs,1) == 2645 & size(i2ra,1) == 2645) 
      s.ra      = [s.ra, aObs(achns,:)];                   % 
      s.ri      = [s.ri, iObs(ichns,:)];      clear rc_ham;
      s.rd      = [s.rd, i2ra(achns,:)];              %
      s.atime   = [s.atime; sno.atim];
      s.itime   = [s.itime; sno.itim];
      s.alat    = [s.alat;  sno.alat];         s.alon = [s.alon;  sno.alon];
      s.ilat    = [s.ilat;  sno.ilat];         s.ilon = [s.ilon;  sno.ilon];
      s.ifov    = [s.ifov;  sno.ifov];
      s.iqual   = [s.iqual; sno.iqual];
      s.tdiff   = [s.tdiff; sno.tdiff];                       %
      s.dist    = [s.dist;  sno.dist'];
      s.alnfr   = [s.alnfr; sno.alandfr];
      %s.l1cr   = [s.l1cr, g.l1cReason'];
      %s.l1cp   = [s.l1cp, g.l1cProc'];
    end
  else
    disp(['Skipping: ' snoLst(ifn).name]);
  end           % if ismember(nbr_rLW)
  fprintf(1,'.');
end             % end for ifn
fprintf(1,'Loaded %d SNO pairs\n',size(s.ilat,1));

s.fa = f2645;   s.fi = f_iasi;  s.fd = fi2a;
s.achns  = achns;
s.ichns  = ichns;
s.dchns  = achns;
%s.aWavs  = aWavs;
%s.iWavs  = iWavs;
%s.agood  = agood;

% Check QA
iok  = find(s.iqual == 0);
ibad = find(s.iqual > 0);

% Remove 6-sigma
aibias = s.rd - s.ra;
disp(['Removing outliers']);
clear gx;
for i=1:length(achns)
   n  = single(remove_6sigma(aibias(i,:)));
   nn = single(remove_6sigma(aibias(i,n)));
   gx(i).n = n(nn);
end

% Now find unique set of bad SNO samples
ux = [];
[~, psz] = size(aibias);
for i=1:length(achns)
   ux = [ux setdiff(1:psz,gx(i).n)];
end
sbad   = unique(ux);
s.r6s  = single(setdiff(1:psz,ux));
disp(['  ' num2str(numel(ux)) ' outliers removed']);
clear gx n nn ux aibias;

% combine iqual and r6s
s.ibad   = sort(unique([ibad; sbad']));
s.iok    = setdiff(1:psz, s.ibad);

%{
% find highest l1cProc value for each channel 
% [0:unchanged, 64:cleaned, see l1cReason, 128:synthesized, 128+1:dummy fill]
% highest l1cReason: [0:preserved, 1:gap, 3, 4, 5, 8, 9, 10, 11, 12, 129:?]
for i = 1:2645
  chanProc(i) = max(s.l1cp(:,i));
  chanReas(i) = max(s.l1cr(:,i));
end
presChanID = find(chanReas == 0);    % what chance - have 2378 channels preserved

% Find the L1b channel IDs corresponding to these L1C preserved chan IDs.
b=sort(f2645(presChanID));
a=sort(fa);
[ai,bi]=seq_match(a,b);     %a(ai) are the L1b set that are preserved

% find which of these apply to the AIRS2CrIS subset.
c = sort(fd);
[ci,bi] = seq_match(c,b);


% Compute averages and standard deviations

ratpm=0; rctpm=0; rdtpm=0; ratps=0; rctps=0; rdtps=0; ratxs=0; rctxs=0; rdtxs=0;
for i = 1:numel(a.nSam)
  ratpm = ratpm + a.avra(:,i).*a.nSam(i);    rctpm = rctpm + a.avrc(:,i).*a.nSam(i);
  rdtpm = rdtpm + a.avrd(:,i).*a.nSam(i);
  ratps = ratps + ( a.sdra(:,i).*a.sdra(:,i) + a.avra(:,i).*a.avra(:,i) )*a.nSam(i);
  rctps = rctps + ( a.sdrc(:,i).*a.sdrc(:,i) + a.avrc(:,i).*a.avrc(:,i) )*a.nSam(i);
  rdtps = rdtps + ( a.sdrd(:,i).*a.sdrd(:,i) + a.avrd(:,i).*a.avrd(:,i) )*a.nSam(i);
end

a.gavrm = ratpm/sum(a.nSam);  a.gcvrm = rctpm/sum(a.nSam);  a.gdvrm = rdtpm/sum(a.nSam);
a.gadrs = real(sqrt( ratps/sum(a.nSam) - a.gavrm.*a.gavrm ));
a.gcdrs = real(sqrt( rctps/sum(a.nSam) - a.gcvrm.*a.gcvrm ));
a.gddrs = real(sqrt( rdtps/sum(a.nSam) - a.gdvrm.*a.gdvrm ));
a.garse = a.gadrs/sqrt(sum(a.nSam));   a.gcrse = a.gcdrs/sqrt(sum(a.nSam));
a.gdrse = a.gddrs/sqrt(sum(a.nSam));
%}

%{
% plot checks
phome   = '/home/chepplew/projects/sno/airs_iasi/figs/';
cyr     = num2str(year(s.dtime1));
[ns nz] = size(s.ra);
abt = real(rad2bt(s.fa(s.achns), s.ra(:,s.iok)));
ibt = real(rad2bt(s.fi(s.ichns), s.ri(:,s.iok)));
dbt = real(rad2bt(s.fa(s.achns), s.rd(:,s.iok)));
abm = nanmean(abt,2);
ibm = nanmean(ibt,2);
dbm = nanmean(dbt,2);
btmbias = nanmean(abt - dbt,2);
btmstd  = nanstd(abt - dbt, 0,2);

ach = find(s.fa(s.achns) > 900,1);
ich = find(s.fi(s.ichns) > 900,1);

dbtbin = 2; btbins = [190:dbtbin:300]; 
btcens = [btbins(1)+dbtbin/2:dbtbin:btbins(end)-dbtbin/2];
pdf_abt = histcounts(abt(ach,:),btbins);
pdf_ibt = histcounts(ibt(ich,:),btbins);
pdf_dbt = histcounts(dbt(ach,:),btbins);

pdf_btbias = histcounts(abt - dbt, [-20:0.2:20]);                  % All channels
pdf_btbias = histcounts(abt(ach,:) - dbt(ach,:), [-20:0.2:20]);    % Single channel

figure(1);clf;simplemap(s.ilat, s.ilon, s.tdiff*24*60);
  title([cyr ' AI ASL SNO Obs delay minutes']);
  %aslprint([phome cyr '_AI_l1c_ASL_SNO_map_delay.png']);
figure(1);clf;simplemap(s.ilat, s.ilon, s.dist); 
figure(1);clf;simplemap(s.ilat, s.ilon, ibt(ich,:)');
figure(1);clf;simplemap(s.ilat, s.ilon, (abt(ach,:) - dbt(ach,:))');

figure(2);plot(btcens,pdf_abt,'.-', btcens,pdf_ibt,'.-');legend('AIRS','IASI');
  xlabel('BT bin (K)');ylabel('bin population');grid on;
  title([cyr ' AI SNO 900wn channel population']);
  %aslprint([phome cyr '_AI_l1c_ASL_SNO_900wn_pdf.png']);

figure(3);clf;plot(s.fa(s.achns),abm,'-', s.fi(s.ichns),ibm,'-', ...
          s.fd(s.achns),dbm,'-');grid on;legend('AIRS','IASI','I2A')
  xlabel('wavenumber cm^{-1}'); ylabel('BT (K)'); xlim([640 1100])
  title([cyr ' AIRS:IASI SNO Mean BT (LW) ' num2str(nz) ' samples']);
  % aslprint([phome cyr '_AI_l1c_ASL_SNO_meanBT_LW.png']);

figure(4);clf;semilogy([-19.9:0.2:19.9], pdf_btbias,'.-');grid on;axis([-16 16 10 1E7]);
  xlabel('BT bin (K)');ylabel('Bin population');
  title([cyr ' SNO AIRS minus IASI all LW channels']);
  %aslprint([phome cyr '_AI_l1c_ASL_SNO_allchans_LW_pdf.png'])
  
figure(5);clf;  % plotxx(s.fa(s.achns),btmbias,'-',s.fa(s.achns),btmstd/sqrt(nz),'-');
[hax,hl1,hl2] = plotyy(s.fa(s.achns),btmbias,s.fa(s.achns),btmstd/sqrt(nz));grid on;
   hax(1).YTick = [-0.6:0.2:0.6];ylim(hax(2),[0 0.08]);
   xlim(hax(1),[640 1100]); xlim(hax(2),[640 1100]);
   ya2=ylabel(hax(2),'Std. error (K)');set(ya2, 'Units', 'Normalized', 'Position', [1.05, 0.7, 0]);
   ylabel(hax(1),'AIRS - IASI (K)');legend('Mean bias AIRS - IASI','Std.Err AIRS - IASI');
   
figure(5);clf;hax = gca;
  yyaxis left;  ha1 = plot(s.fa(s.achns),btmbias,'-');      ylabel('AIRS - IASI (K)');
  yyaxis right; ha2 = plot(s.fa(s.achns),btmstd/sqrt(nz));  ylabel('std. error (K)')
  yyaxis right; hax.YLim = ([0 0.08]);  hax.XLim=([1100 1650]);grid on;
  xlabel('wavenumber cm^{-1}');title([cyr ' AIRS IASI SNO Mean Bias, std.Err']);

   %aslprint([phome cyr '_AI_l1c_ASL_SNO_meanBias_stdErr_MW.png']);
%}
