function [s a] = load_sno_airs_cris_jpl_mat(sdate1, sdate2, cchns)
%
% function load_sno_airs_cris_jpl_mat() loads up radiances for a selected number
%   of channels, specified by CrIs channel number, from the JPL SNO mat files 
%   and for the specified year and months. Unlike the sister function 'read_sno...'
%   this function cacluates statistics during load and subsets of CrIS FOV.
%
% Synopsis: load_sno_airs_cris_jpl_mat('date1','date2',[chan1...chan10]);
%           date1: first month as string: 'YYYY/MM/DD'
%           date2: last month as string: 'YYYY/MM/DD'
%           [chan1, ...]: numeric list of CrIS channels to load (max 10).
%           eg [403 499 737 884 905 998 1021 1297]
%
% Output:  Two structures of arrays. 
%             s: the SNO single fields.
%             a: whole spectrum averages and first moment.
%
% Notes: If the CrIS channel is associated with a bad AIRS channel or has been
%        modified by AIRS L1C (clean and fill) then the next neares good channel
%        is substituted.
%
% Dependencies: i) AIRS good channel list; ii) nominal AIRS and CrIS frequency grids.
%    iii) fixed path and file name syntax for SNO files.
%
% Notes: i) No QA is applied. ii) time separation of SNO pairs from file is positive-only
%    so is recomputed here.
%
 
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

xx=load('/asl/data/airs/airs_freq.mat'); fa=xx.freq; clear xx;
xx=load('cris_freq_nogrd.mat'); f_cris=xx.vchan; clear xx;  % 1305 chns (no guard chans)
xx=load('cris_freq_2grd.mat');  fc = xx.vchan; clear xx;    % 1317 chns (12 guard chans)
fd = f_cris([1:1037 1158:1305]);                            % prior knowledge from Howards decon routine

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
s.Wavs  = cWavs; 
s.sWavs = sWavs;

% ************* load up SNO data ********************
dp     = '/asl/s1/chepplew/projects/sno/airs_cris/v10_0_0/standard/'; % standard/';
snoLst = dir(strcat(dp,'sno_airs_cris_*.mat'));
fprintf(1,'Found %d total SNO files\n',numel(snoLst));

% subset range by date as requested:
dstart = datenum([nyr1 nmn1 ndy1]);
dlast  = datenum([nyr2 nmn2 ndy2]);
for i=1:numel(snoLst)
  junk = snoLst(i).name(15:22);
  thisdat = datenum( [str2num(junk(1:4)) str2num(junk(5:6)) str2num(junk(7:8))] );
  if(thisdat <= dstart) ifn1 = i; end
  if(thisdat <= dlast) ifn2 = i; end
end
fprintf(1,'Processing SNO files from: %s to %s\n',snoLst(ifn1).name, snoLst(ifn2).name);

% Load requested SNO data
s.td   = [];   s.arad = [;]; s.crad = [;]; s.drad = [;]; s.ctim = [];  s.atim = []; 
s.arlat = []; s.arlon = [];  s.dsn  = []; s.crlat = []; s.crlon = []; s.csolz = [];  
s.alnfr = []; s.clnfr = [];  s.cifv  = [];
a.nSam = [];  a.avrd = [;]; a.avra = [;]; a.avrc = [;]; a.sdra = [;];  a.sdrc = [;]; 
a.sdrd = [;];
a.avrx = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
a.cvrx = a.avrx; a.dvrx = a.avrx; a.avrs = a.avrx; a.cvrs = a.avrx;  a.dvrs = a.avrx; 
a.s_nSam = a.avrx; SNames = fieldnames(a.avrx);

for ifnum = ifn1:ifn2
  vars = whos('-file',strcat(dp,snoLst(ifnum).name));
  if( ismember('raDecon', {vars.name}) )              % AIRS->CrIS is present
    if(snoLst(ifnum).bytes > 1.0E4)
      g = load(strcat(dp,snoLst(ifnum).name));
      if(size(g.alat,1) < size(g.alat,2)) fprintf(1,'Unexpected row/column swapped\n'); end
      s.arad  = [s.arad, g.ra(achn,:)];                % [arad, [ra(achn,:); avaw]]; etc
      s.crad  = [s.crad, g.rc(cchn,:)];                % 1317 chns (12 guard chans)
      s.drad  = [s.drad, g.raDecon(dchn,:)];           %
      s.atim  = [s.atim,g.atime'];
      s.ctim  = [s.ctim,g.ctime'];
      s.arlat = [s.arlat,g.alat'];        s.arlon = [s.arlon,g.alon'];
      s.crlat = [s.crlat,g.clat'];        s.crlon = [s.crlon,g.clon'];
      s.cifv  = [s.cifv,g.cifov'];
      s.td    = [s.td,g.tdiff'];                        % use [x,xd'] for JPL SNO
      s.dsn   = [s.dsn,g.dist'];
      s.csolz = [s.csolz,g.csolzen'];
      s.alnfr = [s.alnfr, g.alandfrac'];  s.clnfr = [s.clnfr,g.clandfrac']; 

      a.nSam  = [a.nSam,single(size(g.ra,2))];
      a.avra  = [a.avra,nanmean(g.ra,2)];       a.sdra = [a.sdra,nanstd(g.ra,1,2)];  
      a.avrc  = [a.avrc,nanmean(g.rc,2)];       a.sdrc = [a.sdrc,nanstd(g.rc,1,2)];
      a.avrd  = [a.avrd,nanmean(g.raDecon,2)];  a.sdrd = [a.sdrd,nanstd(g.raDecon,1,2)];

      for k=1:9 incf{k} = find(g.cifov == k); 
        a.avrx.(SNames{k}) = [a.avrx.(SNames{k}), nanmean(g.ra(:,incf{k}),2)]; 
        a.cvrx.(SNames{k}) = [a.cvrx.(SNames{k}), nanmean(g.rc(:,incf{k}),2)]; 
        a.dvrx.(SNames{k}) = [a.dvrx.(SNames{k}), nanmean(g.raDecon(:,incf{k}),2)];
        a.avrs.(SNames{k}) = [a.avrs.(SNames{k}), nanstd(g.ra(:,incf{k}),1,2)]; 
        a.cvrs.(SNames{k}) = [a.cvrs.(SNames{k}), nanstd(g.rc(:,incf{k}),1,2)]; 
        a.dvrs.(SNames{k}) = [a.dvrs.(SNames{k}), nanstd(g.raDecon(:,incf{k}),1,2)]; 
	a.s_nSam.(SNames{k}) = [a.s_nSam.(SNames{k}), size(g.ra(:,incf{k}),2)]; 
      end

    else
      fprintf('%d insufficient number samples\n',ifnum)
    end
  else
    fprintf('skip %s ',snoLst(ifnum).name(15:22) );
  end
  fprintf('.');
end    
fprintf('\n');
s.td2 = s.atim - s.ctim;              % tdiff from JPL are exclusively real positive -WRONG

% Compute averages and standard deviations

azero=zeros(2378,1); czero=zeros(1317,1); dzero=zeros(1185,1);
ratpm=0; rctpm=0; rdtpm=0; ratps=0; rctps=0; rdtps=0; ratxs=0; rctxs=0; rdtxs=0;
ratxm = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
rctxm = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
rdtxm = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
ratxs = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
rctxs = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
rdtxs = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
for k=1:9 
  ratxm.(SNames{k}) = azero;   ratxs.(SNames{k}) = azero;
  rctxm.(SNames{k}) = czero;   rctxs.(SNames{k}) = czero;
  rdtxm.(SNames{k}) = dzero;   rdtxs.(SNames{k}) = dzero;
end
for i = 1:numel(a.nSam)
  ratpm = ratpm + a.avra(:,i).*a.nSam(i);    rctpm = rctpm + a.avrc(:,i).*a.nSam(i);
  rdtpm = rdtpm + a.avrd(:,i).*a.nSam(i);
  ratps = ratps + ( a.sdra(:,i).*a.sdra(:,i) + a.avra(:,i).*a.avra(:,i) )*a.nSam(i);
  rctps = rctps + ( a.sdrc(:,i).*a.sdrc(:,i) + a.avrc(:,i).*a.avrc(:,i) )*a.nSam(i);
  rdtps = rdtps + ( a.sdrd(:,i).*a.sdrd(:,i) + a.avrd(:,i).*a.avrd(:,i) )*a.nSam(i);
  for k = 1:9;
    ratxm.(SNames{k}) = ratxm.(SNames{k}) + a.avrx.(SNames{k})(:,i).*a.s_nSam.(SNames{k})(i);    
    rctxm.(SNames{k}) = rctxm.(SNames{k}) + a.cvrx.(SNames{k})(:,i).*a.s_nSam.(SNames{k})(i);
    rdtxm.(SNames{k}) = rdtxm.(SNames{k}) + a.dvrx.(SNames{k})(:,i).*a.s_nSam.(SNames{k})(i);
    ratxs.(SNames{k}) = ratxs.(SNames{k}) + ( a.avrs.(SNames{k})(:,i).^2  + ...
                         a.avrx.(SNames{k})(:,i).^2 ).*a.s_nSam.(SNames{k})(i);
    rctxs.(SNames{k}) = rctxs.(SNames{k}) + ( a.cvrs.(SNames{k})(:,i).^2  + ...
                         a.cvrx.(SNames{k})(:,i).^2 ).*a.s_nSam.(SNames{k})(i);
    rdtxs.(SNames{k}) = rdtxs.(SNames{k}) + ( a.dvrs.(SNames{k})(:,i).^2  + ...
                         a.dvrx.(SNames{k})(:,i).^2 ).*a.s_nSam.(SNames{k})(i);
  end
end
%
a.ragxm = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
a.rcgxm = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
a.rdgxm = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
%
a.gavrm = ratpm/sum(a.nSam);  a.gcvrm = rctpm/sum(a.nSam);  a.gdvrm = rdtpm/sum(a.nSam);
for k = 1:9
  a.ragxm.(SNames{k}) = ratxm.(SNames{k})/sum( a.s_nSam.(SNames{k})(:) );
  a.rcgxm.(SNames{k}) = rctxm.(SNames{k})/sum( a.s_nSam.(SNames{k})(:) );
  a.rdgxm.(SNames{k}) = rdtxm.(SNames{k})/sum( a.s_nSam.(SNames{k})(:) );
end
%
a.ragxs = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
a.rcgxs = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
a.rdgxs = struct('A',[],'B',[],'C',[],'D',[],'E',[],'F',[],'G',[],'H',[],'I',[]);
%
a.gadrs = real(sqrt( ratps/sum(a.nSam) - a.gavrm.*a.gavrm ));
a.gcdrs = real(sqrt( rctps/sum(a.nSam) - a.gcvrm.*a.gcvrm ));
a.gddrs = real(sqrt( rdtps/sum(a.nSam) - a.gdvrm.*a.gdvrm ));
a.garse = a.gadrs/sqrt(sum(a.nSam));   a.gcrse = a.gcdrs/sqrt(sum(a.nSam));
a.gdrse = a.gddrs/sqrt(sum(a.nSam));
for k = 1:9
  a.ragxs.(SNames{k}) = real(sqrt( ratxs.(SNames{k})/sum(a.s_nSam.(SNames{k})) - ...
                      a.ragxm.(SNames{k}).^2 ));
  a.ragxe.(SNames{k}) = a.ragxs.(SNames{k})./sqrt(sum(a.s_nSam.(SNames{k})));
  a.rcgxs.(SNames{k}) = real(sqrt( rctxs.(SNames{k})/sum(a.s_nSam.(SNames{k})) - ...
                      a.rcgxm.(SNames{k}).^2 ));
  a.rcgxe.(SNames{k}) = a.rcgxs.(SNames{k})./sqrt(sum(a.s_nSam.(SNames{k})));
  a.rdgxs.(SNames{k}) = real(sqrt( rdtxs.(SNames{k})/sum(a.s_nSam.(SNames{k})) - ...
                      a.rdgxm.(SNames{k}).^2 ));
  a.rdgxe.(SNames{k}) = a.rdgxs.(SNames{k})./sqrt(sum(a.s_nSam.(SNames{k})));
end


%{
% plot options
addpath /asl/matlab2012/aslutil                  % simplemap.m
junk = ones(1,size(s.td2,2)); junk = s.td2*86400; junk = s.dsn;
figure(1);clf;simplemap(s.arlat, s.arlon, junk); title('AIRS CRIS SNO map 2013 w/delay (s)');
  % aslprint('AirsCrisLrC_SNO_2013_locMap_wDelay.png');
figure(2);clf;semilogy(fa,a.avra(:,1),'b',fc,a.avrc(:,1),'m',fd,a.avrd(:,1),'c');grid;
figure(3);clf;hist(s.arlat,180),xlim([-90 90]);
   xlabel('latitude bin');ylabel('population');title('AIRS CRIS stnd SNO All data Distribution');
   % aslprint('AirsCris_AllSno_lat_hist.png')

abtm  = real(rad2bt(fa,a.gavrm));
cbtm  = real(rad2bt(fc,a.gcvrm));
dbtm  = real(rad2bt(fd,a.gdvrm));
incd  = find(ismember(fc, fd));
biasm = dbtm - cbtm(incd);

btm   = 0.5*(dbtm + cbtm(incd));
mdr   = 1E-3*(1./drdbt(fd,btm') );
drse  = sqrt((gdrsd.^2 + gcrsd(incd).^2))/sqrt(sum(sd.nSam));
dbse  = mdr.*drse;



figure(3);clf;hist(abt(4,:),100); xlabel('B.T. (K)');ylabel('population');
   title('AIRS All Sample SNO 1414 wn population'); % set(gca, 'YScale', 'log')
   % aslprint('AirsCris_AllSno_1414wn_hist.png')
figure(3);clf;plot(fa,abt,'b',fc,cbt,'c',fd,dbt,'g');
figure(3);clf;plot(fd,biasm,'m');axis([640 2650 -1 1]);grid;
%}

end
