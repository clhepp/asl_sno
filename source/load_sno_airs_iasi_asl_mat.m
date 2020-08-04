function [s] = load_sno_airs_iasi_asl_mat_v2(sdate, edate, xchns, src)
%
% function load_sno_airs_iasi_asl_mat() loads up radiances for a selected number
%   of channels, specified by AIRS channel number, from the ASL SNO mat files 
%   and for the specified year and months. 
%
% Synopsis: load_sno_airs_iasi_asl_mat('sdate','edate',[chan1...chan10], src);
%
% INPUTS    1. sdate: start date as string: 'YYYY/MM/DD'
%           2. edate: end date as string:   'YYYY/MM/DD'
%              N.B. Only accepts the same year.
%           3. xchns: numeric IDs of AIRS L1C channels to load. 
%             eg [403 499 737 884 905 998 1021 1297] or [566:824];
%                [   1:1269] (645  - 1100) cm-1
%                [1270:2162] (1100 - 1615) cm-1
%                [2163:2645] (2181 - 2665) cm-1.
%           4. src: IASI mission number [1, 2 or 3] IASI-1 or -2, -3
%
% Output:  structure of arrays. 
%             s: the SNO single fields.
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

cd /home/chepplew/projects/sno/airs_iasi/

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
posYrs = [2002:2019];
posMns = [1:12];
whos sdate; disp([sdate ' to ' edate]); fprintf('\n');
try 
   D1 = datenum(sdate,'yyyy/mm/dd');
   D2 = datenum(edate,'yyyy/mm/dd');
catch
   error('Incorrect Date Format')
   return
end
s.sdtime = datetime(D1,'convertFrom','datenum');
[nYr1 nMn1 nDy1] = datevec(D1);
[nYr2 nMn2 nDy2] = datevec(D2);
if(nYr1 ~= nYr2) error('Use same year only'); return; end
cYr1   = sdate(1:4);     cMn1 = sdate(6:7);     cDy1 = sdate(9:10);
cYr2   = edate(1:4);     cMn2 = edate(6:7);     cDy2 = sdate(9:10);

  junk = sprintf('%4d/%02d/%02d',nYr1-1,12,31);
jdy1   = datenum(sdate)-datenum(junk);  clear junk;           % needed for data directory
  junk = sprintf('%4d/%02d/%02d',nYr2-1,12,31);
jdy2   = datenum(edate)-datenum(junk);  clear junk;           % needed for data directory
s.sdate = sdate;
s.edate = edate;

% Check channel numbers entered correctly
if(length(xchns) > 20 || length(xchns) < 1 ) fprintf(1,'Wrong number channels\n'); end
if(min(xchns) < 1 || max(xchns) > 1317 ) fprintf(1,'Wrong channel numbers used\n'); end

% Check IASI Mission source
if(~ismember(src,[1,2,3])) error('Invalid IASAI mission number [1 or 2]'); return; end
if(src == 1) IX = '';  IR = 'M02'; end
if(src == 2) IX = '2'; IR = 'M01'; end
if(src == 3) IX = '3'; IR = 'M03'; end
s.src = src;

% load IASI and AIRS channels & good AIRS channels (nig) to use, & bad (nib) to avoid
load('/home/chepplew/projects/iasi/f_iasi.mat');               % f_iasi [8641 x 1]
load('/home/chepplew/projects/airs/airs_f.mat'); f2645 = fairs;  f2378 = f;
load('/home/chepplew/projects/airs/master_nig_01_2009.mat');   % nig [1 x 1535]
% nig = importdata('/home/strow/Work/Airs/good_chan_list');
junk = ismember([1:2378], nig);  nib = find(junk == 0);  clear junk;

% Get over-lapping IASI channels
% Get IASI channels that cover the AIRS range
%[zi ichns] = seq_match(f2645(xchns),f_iasi);
ichns = [find(f_iasi >= f2645(xchns(1)),1):find(f_iasi >= f2645(xchns(end)),1)-1];
% if all channels are required (DON'T RUN OUT OF MEMORY!) then:
if(length(xchns) == 2645) ichns = [1:8461]; end
   
% Screen the channel selection for AIRS bad channels and report any:
achns  = xchns';

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
% ************* Get the SNO data for dates range ********************

dp     = ['/home/chepplew/data/sno/airs_iasi' IX '/ASL/' cYr1 '/'];
snoLst = dir(strcat(dp, 'sno_airs_iasi_*_frmL1C.mat'));
fprintf(1,'Found %d total SNO files\n',numel(snoLst));

ifn1 = 1;             % default start with first file unless later.
for i=1:numel(snoLst)
  junk = regexp(snoLst(i).name,'[0-9]','match');
  junk = cell2mat(junk(1:end));                   % no. 4 appears before date
  thisdat = datenum(junk,'yyyymmdd');
  if(thisdat < D1)  ifn1 = i+1; end
  if(thisdat <= D2) ifn2 = i; end
end
disp(['Source dir: ' dp]);
fprintf(1,'Loading %d SNO files from: %s to %s\n',(ifn2-ifn1+1),snoLst(ifn1).name, ...
        snoLst(ifn2).name);
s.dp    = dp;
s.flist = snoLst(ifn1:ifn2);

% -------------------------------------------------------------------------
s.tdiff = [];    s.ra = [];    s.ri = [];    s.ri2a = [];  s.itime = [];  s.atime = []; 
 s.alat = [];  s.alon = []; s.dist  = [];    s.ilat = [];   s.ilon = [];  s.isolz = [];
s.asolz = [];  
s.iqual = []; s.alnfr = [];  s.ifov = [];    s.l1cproc = []; s.l1creas = [];

for ifn = ifn1:ifn2;
  vars = whos('-file',strcat(dp,snoLst(ifn).name));
  if( ismember('ri', {vars.name}) & ismember('ra', {vars.name}) & ...
      ismember('ri2a', {vars.name})  )  
    load(strcat(dp, snoLst(ifn).name));
    if(size(ri,1) == 8461 & size(ra,1) == 2645 & size(ri2a,1) == 2645)
    if(size(ra,2) ~= size(ri2a,2)) disp(['unequal samples: ' num2str(ifn)]); end
      s.ra      = [s.ra,   ra(achns,:)];                   % 
      s.ri      = [s.ri,   ri(ichns,:)];      clear rc_ham;
      s.ri2a    = [s.ri2a, ri2a(achns,:)];              %
      s.atime   = [s.atime; sno.atim];
      s.itime   = [s.itime; sno.itim];
      s.alat    = [s.alat;  sno.alat];         s.alon = [s.alon;  sno.alon];
      s.ilat    = [s.ilat;  sno.ilat];         s.ilon = [s.ilon;  sno.ilon];
      s.ifov    = [s.ifov;  sno.ifov];
      s.iqual   = [s.iqual; sno.iqual];
      s.tdiff   = [s.tdiff; sno.tdiff];                       %
      s.dist    = [s.dist;  sno.dist'];
      s.alnfr   = [s.alnfr; sno.alandfr];
      s.asolz   = [s.asolz; sno.asolzen];    s.isolz = [s.isolz; sno.isolzen];
      s.l1creas = [s.l1creas, l1cSynthReason(achns,:)];
      s.l1cproc = [s.l1cproc, l1cProc(achns,:)];
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

% Check QA and non-physical radiances
clear izx;
iok   = find(s.iqual == 0);
ibad  = find(s.iqual > 0);
for i=1:length(ichns)
  izx(i).n = find(s.ri(i,:) <= 0 );
end

% Remove 6-sigma
aibias = s.ri2a - s.ra;
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


% ******************* END of routine ***********************



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
