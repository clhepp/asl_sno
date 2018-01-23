function [s] = load_sno_iasi_cris_asl_mat(sdate, edate, xchns, src)
%
% function load_sno_iasi_cris_asl_mat() loads up radiances for a selected number
%   of channels, specified by CrIS channel number, from the ASL SNO mat files 
%   and for the specified year and months. 
%   Unlike the sister function 'read_sno...'
%   this function cacluates statistics during load and subsets of CrIS FOV.
%
% Synopsis: load_sno_iasi_cris_asl_mat('date1','date2',[chan1...chan10]);
%           sdate: start date as string: 'YYYY/MM/DD'
%           edate: end date as string:   'YYYY/MM/DD'
%           N.B. Only accepts the same year.
%           xchns: numeric IDs of CrIS channels to load based on NO guard channel list. 
%           (max 10).
%           eg [403 499 737 884 905 998 1021 1297] or [566:824]; 
%              LW: [1:717], MW:[718:1154], SW:[1155:1317];
%                  [1:1269] (645  - 1100) cm-1
%                  [1270:2160] (1100 - 1615) cm-1
%           src:  [M1, M2]. mission numbers: iasi-1 or iasi-2 (MetOp-A or -B) 
%                 and cris-1 or cris-2 (NPP or JPSS-1), 
%
% Output:  Two structures of arrays. 
%             s: the SNO geo, time and related vector fields.
%             a: whole spectrum averages and first moment.
%
%
% Notes: 
%        The CrIS CCAST spectral resolution is hard-wired for LR (low-res).
%        The IASI spectra are apodized.
%
% Dependencies: i) nominal IASI and CrIS w/2 guard channels per edge frequency grids.
%    iii) fixed path and file name syntax for SNO files.
%
% Notes: i) No QA is applied. ii) time separation of SNO pairs from file is positive-only
%    so is recomputed here.
%
% Author: C. L. Hepplewhite, UMBC/JCET
%
% Version: 24-October-2017
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

% Check number of input arguments
if(nargin ~= 4) error('Please enter all 4 input arguments'); return; end

% Check CrIS resolution (high or low)
cris_res='high';
cris_res = upper(cris_res);
if(strcmp(cris_res,'HIGH')) CR='HR'; ncc=2223; end
if(strcmp(cris_res,'LOW'))  CR='LR'; ncc=1317; end

% Check mission numbers
if(length(src) ~=2) error('Need IASI and CRIS mission numbers'); return; end
junk = ismember(src,[1,2]);
if(~all(junk)) error('Mission numbers can only be 1 or 2 for now'); return; end
disp(['you have selected IASI-' num2str(src(1)) ' and CRIS-' num2str(src(2))]);
if(src(1) == 1) IX = '';  end
if(src(1) == 2) IX = '2'; end
if(src(2) == 1) CX = '';  end
if(src(2) == 2) CX = '2'; end

% Process and check the date strings
posYrs = [2002:2018];
posMns = [1:12];
whos sdate; disp([sdate ' to ' edate]); fprintf('\n');
try 
   D1 = datenum(sdate,'yyyy/mm/dd');
   D2 = datenum(edate,'yyyy/mm/dd');
catch
   error('Incorrect Date Format')
   return
end
s.dtime1 = datetime(D1,'convertFrom','datenum');
[nYr1 nMn1 nDy1] = datevec(D1);
[nYr2 nMn2 nDy2] = datevec(D2);
if(nYr1 ~= nYr2) error('Use same year only'); return; end
cYr1   = sdate(1:4);     cMn1 = sdate(6:7);     cDy1 = sdate(9:10);
cYr2   = edate(1:4);     cMn1 = edate(6:7);     cDy1 = edate(9:10);

  junk = sprintf('%4d/%02d/%02d',nYr1-1,12,31);
jdy1   = datenum(sdate)-datenum(junk);  clear junk;           % needed for data directory
  junk = sprintf('%4d/%02d/%02d',nYr2-1,12,31);
jdy2   = datenum(edate)-datenum(junk);  clear junk;           % needed for data directory

% Check channel numbers entered correctly
if(length(xchns) > 20 || length(xchns) < 1 ) fprintf(1,'Wrong number channels\n'); end
if(min(xchns) < 1 || max(xchns) > 1317 ) fprintf(1,'Wrong channel numbers used\n'); end

% load IASI and AIRS channels & good AIRS channels (nig) to use, & bad (nib) to avoid
load('/home/chepplew/projects/iasi/f_iasi.mat');               % f_iasi [8641 x 1]
load('/home/chepplew/projects/cris/cris_freq_2grd.mat'); fc = vchan;

% Get over-lapping IASI channels
%[zi ichns] = seq_match(fc(xchns),f_iasi);
ichns = [find(f_iasi >= fc(xchns(1)),1): find(f_iasi >= fc(xchns(end)),1)]; 
cchns = xchns;
   
% ************* load up SNO data ********************

dp     = ['/home/chepplew/data/sno/iasi' IX '_cris' CX '/ASL/' CR '/' cYr1 '/'];
snoLst = dir(strcat(dp, 'sno_iasi_cris_asl_*.mat'));
fprintf(1,'Found %d total SNO files in %s\n',numel(snoLst), dp);

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
s.tdiff = [];    s.rc = [];    s.ri = [];      s.rd = [];  s.itime = [];  s.ctime = []; 
 s.clat = [];  s.clon = []; s.dist  = [];    s.ilat = [];   s.ilon = [];  s.csolz = [];  
s.iqual = []; s.clnfr = [];  s.ifov = [];    s.cfov = []; s.prcver = [];

for ifn = 1:numel(snoLst)  % ifn1:ifn2;
  vars = whos('-file',strcat(dp,snoLst(ifn).name));
  if( ismember('ri', {vars.name}) & ismember('rc', {vars.name}) & ...
      ismember('ri2c', {vars.name})  )  
    g=load(strcat(dp, snoLst(ifn).name));
    if  (size(g.ri,2) == 8461 & size(g.rc,2) == ncc & size(g.ri2c,1) == ncc) 
        %junk    = single(hamm_app(double(g.rc(:,cchns))'));
      s.rc      = [s.rc, g.rc(:,cchns)'];               % 
      s.ri      = [s.ri, g.ri(:,ichns)'];               %
      s.rd      = [s.rd, g.ri2c(cchns,:)];              %
      s.ctime   = [s.ctime; g.ctime];
      s.itime   = [s.itime; g.itime];
      s.clat    = [s.clat;  g.clat];         s.clon = [s.clon;  g.clon];
      s.ilat    = [s.ilat;  g.ilat];         s.ilon = [s.ilon;  g.ilon];
      s.ifov    = [s.ifov;  g.ifov];
      s.cfov    = [s.cfov;  g.cfov];
      s.csolz   = [s.csolz; g.csolzen];
      s.tdiff   = [s.tdiff; g.tdiff];                       %
      s.dist    = [s.dist;  g.dist];
      %s.prcver  = [s.prcver; g.process_version];
      %s.alnfr   = [s.alnfr; sno.alandfr];
      %s.l1cr    = [s.l1cr, g.l1cReason'];
      %s.l1cp    = [s.l1cp, g.l1cProc'];
    end
  else
    disp(['Skipping: ' snoLst(ifn).name]);
  end           % if ismember(nbr_rLW)
  fprintf(1,'.');
end             % end for ifn
fprintf(1,'Loaded %d SNO pairs\n',size(s.ilat,1));

s.fi = f_iasi;  s.fd = g.fi2c;  s.fc = g.fc;
s.cchns  = cchns;
s.ichns  = ichns';
s.dchns  = cchns;

% Check QA
iok  = ':'; % find(s.iqual == 0);
ibad = [];  % find(s.iqual > 0);

% Remove 6-sigma
icbias = s.rd - s.rc;
disp(['Removing outliers']);
clear gx;
for i=1:length(cchns)
   n  = single(remove_6sigma(icbias(i,:)));
   nn = single(remove_6sigma(icbias(i,n)));
   gx(i).n = n(nn);
end

% Now find unique set of bad SNO samples
ux = [];
[~, psz] = size(icbias);
for i=1:length(cchns)
   ux = [ux setdiff(1:psz,gx(i).n)];
end
sbad   = unique(ux);
s.r6s  = single(setdiff(1:psz,ux));
disp(['  ' num2str(numel(ux)) ' outliers removed']);
clear gx n nn ux icbias;

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
phome   = '/home/chepplew/projects/sno/iasi_cris/figs/';
cyr     = num2str(year(s.dtime1));
[ns nz] = size(s.rc);
cbt = real(rad2bt(s.fc(s.cchns), s.rc(:,s.iok)));
ibt = real(rad2bt(s.fi(s.ichns), s.ri(:,s.iok)));
dbt = real(rad2bt(s.fc(s.cchns), s.rd(:,s.iok)));
cbm = nanmean(cbt,2);
ibm = nanmean(ibt,2);
dbm = nanmean(dbt,2);
btmbias = nanmean(cbt - dbt,2);
btmstd  = nanstd(cbt - dbt, 0,2);
  whos cbt ibt dbt cbm ibm dbm btmbias btmstd;

% in case of LW band:
cch = find(s.fc(s.cchns) > 900,1);
ich = find(s.fi(s.ichns) > 900,1);

dbtbin = 2; btbins = [190:dbtbin:300]; 
btcens = [btbins(1)+dbtbin/2:dbtbin:btbins(end)-dbtbin/2];
pdf_cbt = histcounts(cbt(cch,:),btbins);
pdf_ibt = histcounts(ibt(ich,:),btbins);
pdf_dbt = histcounts(dbt(cch,:),btbins);

pdf_btbias = histcounts(cbt - dbt, [-20:0.2:20]);                  % All channels
pdf_btbias = histcounts(cbt(cch,:) - dbt(cch,:), [-20:0.2:20]);    % Single channel

figure(1);clf;simplemap(s.ilat, s.ilon, s.tdiff*24*60);
  title([cyr ' I1:C ASL SNO Obs delay minutes']);
  %aslprint([phome cyr '_I1C_ASL_SNO_map_delay.png']);
figure(1);clf;simplemap(s.ilat, s.ilon, s.dist); 
figure(1);clf;simplemap(s.ilat, s.ilon, ibt(ich,:)');
figure(1);clf;simplemap(s.ilat, s.ilon, (cbt(cch,:) - dbt(cch,:))');

figure(2);plot(btcens,pdf_cbt,'.-', btcens,pdf_ibt,'.-');legend('CrIS','IASI');
  xlabel('BT bin (K)');ylabel('bin population');grid on;
  title([cyr ' I1:C SNO 900wn channel population']);
  %aslprint([phome cyr '_I1C_ASL_SNO_900wn_pdf.png']);

figure(3);clf;plot(s.fc(s.cchns),cbm,'-', s.fi(s.ichns),ibm,'-', ...
          s.fd(s.cchns),dbm,'-');grid on;legend('CrIS','IASI','I2C')
  xlabel('wavenumber cm^{-1}'); ylabel('BT (K)'); xlim([640 1100]); 
  % xlim([1200 1780]); xlim([2140 2560])
  title([cyr ' IASI1:CrIS SNO Mean BT (LW) ' num2str(nz) ' samples']);
  % aslprint([phome cyr '_I1C_ASL_SNO_meanBT_LW.png']);

figure(4);clf;semilogy([-19.9:0.2:19.9], pdf_btbias,'.-');grid on;axis([-16 16 10 1E7]);
  xlabel('BT bin (K)');ylabel('Bin population');
  title([cyr ' SNO CrIS minus IASI 900 cm-1 channel']);
  %aslprint([phome cyr '_I1C_ASL_SNO_900wn_bias_pdf.png'])
  
figure(5);clf;  % plotxx(s.fa(s.achns),btmbias,'-',s.fa(s.achns),btmstd/sqrt(nz),'-');
[hax,hl1,hl2] = plotyy(s.fc(s.cchns),btmbias,s.fc(s.cchns),btmstd/sqrt(nz));grid on;
   hax(1).YTick = [-0.6:0.2:0.6]; xlim(hax(1),[640 1100]); ylim(hax(1),[-1 1]);
   xlim(hax(2),[640 1100]); ylim(hax(2),[0 0.08]);
   ya2=ylabel(hax(2),'Std. error (K)');set(ya2, 'Units', 'Normalized', 'Position', [1.05, 0.7, 0]);
   xlabel(hax(1),'Wavenumber (cm^{-1})');
   ylabel(hax(1),'CrIS - IASI (K)');legend('Mean bias CrIS-IASI','Std.Err CrIS-IASI');
   title([cyr ' CrIS - IASI SNO Mean Bias, std.Err']);
figure(5);clf;hax = gca;
  yyaxis left;  ha1 = plot(s.fc(s.cchns),btmbias,'-');      ylabel('CrIS - IASI (K)');
  hax.YLim=([-0.8 0.8]);
  yyaxis right; ha2 = plot(s.fc(s.cchns),btmstd/sqrt(nz));  ylabel('std. error (K)')
  yyaxis right; hax.YLim = ([0 0.08]);  hax.XLim=([640 1100]);grid on;
  xlabel('wavenumber cm^{-1}');title([cyr ' CrIS - IASI SNO Mean Bias, std.Err']);
   %aslprint([phome cyr '_I1C_ASL_SNO_meanBias_stdErr_LW.png']);
%}
