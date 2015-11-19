function [s a] = load_sno_airs_cris_asl_mat(sdate1, sdate2, cchns, src)
%
% function load_sno_airs_cris_asl_mat() loads up radiances for a selected number
%   of channels, specified by CrIs channel number, from the JCET SNO mat files 
%   and for the specified year and months. 
%   Options for CCAST LR, CCAST HR, IDPS LR.
%   Unlike the sister function 'read_sno...'
%   this function cacluates statistics during load and subsets of CrIS FOV.
%
% Synopsis: load_sno_airs_cris_jct_mat('date1','date2',[chan1...chan10]);
%           date1: first month as string: 'YYYY/MM/DD'
%           date2: last month as string: 'YYYY/MM/DD'
%           [chan1, ...]: numeric list of CrIS channels to load (max 10).
%           eg [403 499 737 884 905 998 1021 1297]
%           src: source. one of: ' '
%
% Output:  Two structures of arrays. 
%             s: the SNO single fields.
%             a: whole spectrum averages and first moment.
%
%
% Notes: If the CrIS channel is associated with a bad AIRS channel or has been
%        modified by AIRS L1C (clean and fill) then the next nearest good channel
%        is substituted.
%        The CrIS spectra are from the CCAST production and are Sinc (unapodized).
%
% Dependencies: i) AIRS good channel list; ii) nominal AIRS and CrIS frequency grids.
%    iii) fixed path and file name syntax for SNO files.
%
% Notes: i) No QA is applied. ii) time separation of SNO pairs from file is positive-only
%    so is recomputed here.
%
% Author: C. L. Hepplewhite, UMBC/JCET
%
% Version: Initial draft: 02-May-2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/packages/airs_decon/source             % hamm_app.m
addpath /asl/matlab2012/aslutil                     % rad2bt.m

cd /home/chepplew/projects/sno/sno_git_repo/source/


% Process the date strings
posYrs = [2002:2015];
posMns = [1:12];
% check dates are entered correctly:
if (length(sdate1) ~= 10 || length(sdate2) ~= 10) fprintf(1,'Error in dates\n'); exit; end
syr1 = sdate1(1:4);  smn1 = sdate1(6:7);  sdy1 = sdate1(9:10);
syr2 = sdate2(1:4);  smn2 = sdate2(6:7);  sdy2 = sdate2(9:10);
junk = ismember(posYrs, str2num(syr1)); if(isempty(~find(junk))) fprintf('invalid year\n'); end
junk = ismember(posMns, str2num(smn1)); if(isempty(~find(junk))) fprintf('invalid month\n'); end
junk = ismember([1:31], str2num(sdy1)); if(isempty(~find(junk))) fprintf('invalid day\n'); end
junk = ismember(posYrs, str2num(syr2)); if(isempty(~find(junk))) fprintf('invalid year\n'); end
junk = ismember(posMns, str2num(smn2)); if(isempty(~find(junk))) fprintf('invalid month\n'); end
junk = ismember([1:31], str2num(sdy2)); if(isempty(~find(junk))) fprintf('invalid day\n'); end
nyr1 = str2num(syr1); nmn1 = str2num(smn1);  ndy1 = str2num(sdy1);
nyr2 = str2num(syr2); nmn2 = str2num(smn2);  ndy2 = str2num(sdy2);
junk = sprintf('%4d/%02d/%02d',nyr1-1,12,31);
jdy1 = datenum(sdate1)-datenum(junk);  clear junk;           % needed for data directory
junk = sprintf('%4d/%02d/%02d',nyr2-1,12,31);
jdy2 = datenum(sdate2)-datenum(junk);  clear junk;           % needed for data directory

% Check channel numbers entered correctly
if(length(cchns) > 10 || length(cchns) < 1 ) fprintf(1,'Wrong number channels\n'); end
if(min(cchns) < 1 || max(cchns) > 1317 ) fprintf(1,'Wrong channel numbers used\n'); end

% get list of good AIRS channels (nig) to use, & bad (nib) to avoid
load('/home/chepplew/projects/airs/master_nig_01_2009.mat');   % nig [1 x 1535]
% nig = importdata('/home/strow/Work/Airs/good_chan_list');
junk = ismember([1:2378], nig);  nib = find(junk == 0);  clear junk;

load('/asl/data/airs/airs_freq.mat'); fa=freq;  clear freq;
xx = load('cris_freq_nogrd.mat'); fcris = xx.vchan;    % 1305 chns (no guard chans)
xx = load('cris_freq_2grd.mat');  fc = xx.vchan;       % 1317 chns (12 guard chans)
fd = fcris([1:1037 1158:1305]);                        % prior knowledge from Howards decon routine
FH = fopen('freq2645.txt','r');
  junk = textscan(FH,'%f');
  f2645 = single(cell2mat(junk));
fclose(FH);
  
% Screen the channel selection for AIRS bad channels and update if necessary:
cWavs = fc(cchns);
for i=1:numel(cWavs)
  tmp     = find(fa  > cWavs(i)-0.1,1);
  aref    = find(nig > tmp,1);  achn(i) = nig(aref);
end
cWavs = fa(achn);
for i=1:numel(cWavs)
  cchn(i) = find(fc  > cWavs(i)-0.25, 1);
  dchn(i) = find(fd  > cWavs(i)-0.25, 1);
end
for i=1:numel(cWavs) sWavs{i}  = sprintf('%6.2f',cWavs(i)); end

% ************* load up SNO data ********************

dp = '/asl/s1/chepplew/projects/sno/airs_cris/LR/';
snoLst = dir(strcat(dp, 'sno_airs_crisLR_clh_*_018d600s.mat'));
fprintf(1,'Found %d total SNO files\n',numel(snoLst));

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

%%%%%%%%%
s.td    = [];  s.arad = [;]; s.crad = [;]; s.drad = [;]; s.ctim = [];  s.atim = []; 
s.arlat = []; s.arlon = [];  s.dsn  = []; s.crlat = []; s.crlon = []; s.csolz = [];  
s.alnfr = []; s.clnfr = [];  s.cifv = [];  s.l1cr = [;]; s.l1cp = [;];
a.nSam  = [];  a.avrd = [;]; a.avra = [;]; a.avrc = [;]; a.sdra = [;]; a.sdrc = [;]; 
a.sdrd  = [;]; 

for ifnum = ifn1:ifn2;
  vars = whos('-file',strcat(dp,snoLst(ifnum).name));
  if( ismember('ra2c', {vars.name}) )              % AIRS->CrIS is present
    g = load(strcat(dp, snoLst(ifnum).name));
    %if(size(g.alat,1) < size(g.alat,2)) fprintf(1,'Unexpected row/column swapped\n'); end

      s.arad  = [s.arad, g.ra(achn,:)];                % [arad, [ra(achn,:); avaw]]; etc
         junk = double(g.rc(cchn,:));                  % unapodized, double for hamm_app
      rc_hamm = single(hamm_app(junk)); clear junk;
      s.crad  = [s.crad, rc_hamm];                     % 1317 chns (12 guard chans)
      s.drad  = [s.drad, g.ra2c(dchn,:)];              %
         clear rc_hamm;
      s.atim  = [s.atim, g.atime];
      s.ctim  = [s.ctim, g.ctime];
      s.arlat = [s.arlat,g.alat];         s.arlon = [s.arlon,g.alon];
      s.crlat = [s.crlat,g.clat];         s.crlon = [s.crlon,g.clon];
      s.cifv  = [s.cifv, g.cifov];
      s.td    = [s.td,   g.tdiff];                       %
      s.dsn   = [s.dsn,  g.dist];
      s.l1cr  = [s.l1cr, g.l1cReason'];
      s.l1cp  = [s.l1cp, g.l1cProc'];
            
      junk    = double(g.rc);
      rc_hamm = hamm_app(junk);
      a.nSam  = [a.nSam,single(size(g.ra,2))];
      a.avra  = [a.avra,nanmean(g.ra,2)];       a.sdra = [a.sdra,nanstd(g.ra,1,2)];  
      a.avrc  = [a.avrc,nanmean(rc_hamm,2)];    a.sdrc = [a.sdrc,nanstd(rc_hamm,1,2)];
      a.avrd  = [a.avrd,nanmean(g.ra2c,2)];     a.sdrd = [a.sdrd,nanstd(g.ra2c,1,2)];

         clear rc_hamm;
  end    % if ismember(ra2c)
  fprintf(1,'.');
end

% find highest l1cProc value for each channel 
% [0:unchanged, 64:cleaned, see l1cReason, 128:synthesized, 128+1:dummy fill]
% highest l1cReason: [0:preserved, 1:gap, 3, 4, 5, 8, 9, 10, 11, 12, 129:?]
for i = 1:2645
  chanProc(i) = max(s.l1cp(:,i));
  chanReas(i) = max(s.l1cr(:,i));
end
presChanID = find(chanReas == 0);    % what chance - have 2378 channels preserved

%{
% Find the L1b channel IDs corresponding to these L1C preserved chan IDs.
b=sort(f2645(presChanID));
a=sort(fa);
[ai,bi]=seq_match(a,b);     %a(ai) are the L1b set that are preserved

% find which of these apply to the AIRS2CrIS subset.
c = sort(fd);
[ci,bi] = seq_match(c,b);
%}

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




%{
% plot options

abtm  = real(rad2bt(fa,a.avra));
cbtm  = real(rad2bt(fc,a.avrc));
dbtm  = real(rad2bt(fd,a.avrd));
incd  = find(ismember(fc, fd));
biasm = dbtm - cbtm(incd);

figure(1);clf;plot(fa,abtm,'b',fc,cbtm,'r',fd,dbtm,'g');grid;axis([640 2700 200 290]);
figure(2);clf;h1=subplot(2,1,1);plot(fd,biasm,'m');grid;axis([640 2400 -1 1]);
   h2=subplot(2,1,2);plot(fd(ci),biasm(ci),'m.');grid;axis([640 2400 -0.5 0.5]);
   linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp=get(h1,'position');
   set(h1,'position',[pp(1) pp(2)-pp(4)*0.1 pp(3) pp(4)*1.1])
   pp=get(h2,'position'); set(h2,'position',[pp(1) pp(2) pp(3) pp(4)*1.1]);

abtm  = real(rad2bt(fa,a.gavrm));
cbtm  = real(rad2bt(fc,a.gcvrm));
dbtm  = real(rad2bt(fd,a.gdvrm));
incd  = find(ismember(fc, fd));
biasm = dbtm - cbtm(incd);

btm   = 0.5*(dbtm + cbtm(incd));
mdr   = 1E-3*(1./drdbt(fd,btm') );
drse  = sqrt((a.gddrs.^2 + a.gcdrs(incd).^2))/sqrt(sum(a.nSam));
dbse  = mdr.*drse';

figure(3);clf;plot(fa(nig),abtm(nig),'b',fc,cbtm,'r',fd,dbtm,'g');axis([640 2700 210 270]);
  grid;xlabel('wavenumber cm^{-1}');ylabel('BT K');
  title('Airs (b) CrIS (r) A2C (g) SNO 2012-15 All scenes')
figure(3);clf;plot(fd,biasm,'m');axis([640 2650 -0.6 0.6]);grid;


%}
