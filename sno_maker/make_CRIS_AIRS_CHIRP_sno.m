function [] = make_CRIS_AIRS_CHIRP_sno(rdate, cris_res, cris, opt1, vers)
%
% function make_CRIS_AIRS_CHIRP_sno(rdate, cris_res, cris)
% 
% to produce files of SNO using CRIS L1C from CCAST processing at UMBC.ASL 
%   and CHIRP generated from AIRS L1C files.
% 
% INPUTS: rdate:    date to processes: string: 'YYYY/MM/DD'
%         cris_res: two strings: [{'high','mid','low'}{'high','mid','low'}]
%                   1st: CrIS source SDR spectral resolution.
%                   2nd: CrIS SNO spectral resolution (limits apply)
%         cris:     numeric: [1 or 2] CrIS mission NPP (1) or JPSS-1 (2)
%         opt1:     structure of fields to control airs2cris conversion,
%                   if entered with no fields no conversion takes place.
%         vers:     Based on CHIRP version & reference to add to the SNO mat file.
%                   maximum allowed length 6 characters. 
%
% 1. Sources:
%    CHIRP: /asl/hpcnfs1/chirp/airs_L1c/
%    CRIS:  /asl/cris/ccast/sdr45_{npp,j01}_HR/YYYY/JJJ/SDR_d20130828_t0006589.mat
% 
% Notes:  HM: L1a_err, in the SDR files is a 30 x nscan array, one value for each FOR:
%         1 = bad, 0 = OK.   
%
% version 2: retain FOR -by- x-track -by- a-track array structures
%            include neighbours.

cd /home/chepplew/projects/sno/makeSNO;
addpath /home/chepplew/projects/sno/makeSNO/src/       % read_netcdf_lls.m
addpath /asl/matlib/time                               % tai2dnum, airs2dnum
addpath /home/chepplew/myLib/matlib                    % read_airs_l1c
addpath /asl/rtp_prod/airs/utils                       % f_default_1lb.mat

addpath /home/chepplew/gitLib/airs_deconv/test         % k2cris.m
addpath /asl/packages/ccast/source                     % inst_params.m
addpath /asl/matlib/h4tools                            % h4sdread.m
addpath /asl/packages/airs_decon/source                % airs2cris.m

warning('off');
run cris_neighbor_LUT                                % load CrIS adjacent FOV look-up tables

% Check all input parameters are entered.
if(nargin ~= 5) error('Please enter all 5 input arguments'); return; end

% Check date string and record the following day
whos rdate; disp(rdate); fprintf('\n');
try 
   D = datenum(rdate,'yyyy/mm/dd');
   Dp1 = D+1;
catch
   error('Incorrect Date Format')
   return
end

% Check CrIS Resolution
cris_res = upper(cris_res);
all_res  = {'HIGH','MID','LOW'};
if(~ismember(cris_res,all_res))
  error('Invalid CrIS resolution. Choices are: high mid low');
  return
end
  
% Check CrIS mission (NPP=1, J-1=2)
if(~ismember(cris,[1,2])) error('CrIS can only take values 1 or 2'); return; end
if(cris == 1); CX = '';  NX = 'npp'; end
if(cris == 2); CX = '2'; NX = 'j01'; end

% check version reference
vers = lower(vers);
if(length(vers) > 8) error('Version reference too long, 6 chars max'); return; end
allvers = {'v01a',''};
if(~ismember(vers,allvers)) error('Not valid CHIRP data version');return; end
 
% slelect source and destination directories.   
switch cris_res{1}
  case 'HIGH'
    CR = 'HR';
  case 'MID'
    CR = 'MR';
  case 'LOW'
    CR = 'LR';
end
switch cris_res{2}
  case 'HIGH' 
    savDir = ['/asl/xfs3/sno/chirp_cris' CX '/ASL/HR/'];
    if(strcmp(vers,'noaa'))
       savDir = ['/asl/xfs3/sno/chirp_cris' CX '/NOAA_pon/HR/'];
    end
    npLW = 717;  npMW = 869;   npSW = 637;
  case 'MID' 
    savDir = ['/asl/xfs3/sno/chirp_cris' CX '/ASL/MR/'];
    npLW = 717;  npMW = 865;   npSW = 633;
  case 'LOW'
    savDir = ['/asl/xfs3/sno/chirp_cris' CX '/ASL/LR/'];
    npLW = 717;  npMW = 437;   npSW = 163;
end

% Check conversion options, names and values
%opt1  = struct;
opt1
names = fieldnames(opt1);
if(isempty(names)) clear opt1; 
elseif(~all(isfield(opt1,{'hapod','inst_res','user_res','nguard','bfile','scorr','dvb'})))
  disp('Not all options for airs2cris are present, plz chceck');
end

% Criteria for separation and delay
maxDtim   = 0.006944;             % day. 600/86400 secs =  Mattime is decimal day (20 mins=0.0139)
maxDphi   = 0.0719;               % deg. 0.0719 = 8.0 km. Earth radius: 6371 km. 1 deg = 111 km.

% Get day to process & convert to required formats
[nYr1 nMn1 nDy1] = datevec(D);
[nYr2 nMn2 nDy2] = datevec(Dp1);
cYr     = rdate(1:4);     cMn = rdate(6:7);     cDy = rdate(9:10);
day0    = sprintf('%4d/%02d/%02d',nYr1 - 1,12,31);
jday    = datenum(rdate) - datenum(day0);  clear junk;           % needed for CRIS directory
jday2   = fix(Dp1 - datenum(day0));

% ************    Get CRIS SDR granule files  *****************
if(isempty(CX) )
  %crDir  = ['/asl/data/cris/ccast/sdr60_npp_' CR '/'];               % pre-2018
  %crDir  = ['/asl/data/cris/ccast/sdr45_npp_' CR '/'];               % from 2018
  crDir  = ['/asl/cris/ccast/sdr45_' NX '_HR/'];                     % from Dec 2018
  crDir1 = sprintf('%s%s/%03d/',crDir, cYr,jday);
  crLst  = dir(strcat(crDir1, 'CrIS_SDR_npp_s45_d', cYr, cMn, cDy,'_t*', 'v20a.mat'));    % from 2018
  disp([strcat(crDir1, 'SDR_d', cYr, cMn, cDy,'_t*.mat')]);
end
if(CX == '2')
  crDir  = ['/asl/cris/ccast/sdr45_j01_' CR '/'];
  crDir1 = sprintf('%s%s/%03d/',crDir, cYr,jday);
  disp(crDir1)
  crLst  = dir(strcat(crDir1, 'CrIS_SDR_j01_s45_d', cYr, cMn, cDy,'_t*.mat'));
  if(strcmp(vers,'noaa'))
    crDir  =  '/home/chepplew/data/cris/noaa/sdr_j01_hr/';
    crDir1 = sprintf('%s%s/%03d/',crDir, cYr,jday);
    crLst  = dir(strcat(crDir1, 'CrIS_SDR_j01_s45_d', cYr, cMn, cDy,'_tXXX_noaa_PLon_*.mat'));
  end
end
disp(['Found ' num2str(numel(crLst)) ' CrIS L1c Files']);
if(numel(crLst) < 1) error('Insufficient CrIS CCAST data'); return; end

disp([crLst(1).folder '  ' crLst(1).name])
disp([' cris mission: ' num2str(cris) '. CHIRP vers: ' vers])


% Get first granule of the following day and append to full list:
crDir2 = sprintf('%s%s/%03d/',crDir, cYr,jday2);
crLst2 = dir([crDir2 'SDR_d', cYr, sprintf('%02d', nMn2), sprintf('%02d',nDy2),'_t*.mat']);
for i=1:numel(crLst) cFnams{i} = strcat(crLst(i).folder, '/', crLst(i).name); end

% add one granule from next day (duplicate last one of the day
cFnams{numel(crLst)+1} = strcat(crLst(end).folder,'/',crLst(end).name);
 
clat = [];   clon = []; ctim = [];  cfov = [];  cxtrak = []; catrak = []; cgran = [];
csozn = []; csazn = [];
for fn = 1:numel(crLst)
  load(cFnams{fn});
% Check granule and generate FOV, a-track, x-track indexes - only works with fixed array size
  if(ndims(geo.Latitude) ~= 3) fprintf('ERROR: wrong ndims of geo.Lat\n'); end
  [sz1 sz2 sz3] = size(geo.Latitude);
  disp([num2str(fn) ': ' num2str(sz1) '  ' num2str(sz2) '  ' num2str(sz3)])
  if(sz1 ~= 9 | mod(sz2,30) ~= 0) fprintf('ERROR: granule size is wrong\n'); end
% granule size is okay:- proceed;
  count = 0; clear tmpv tmpx tmpa;
  for m=1:sz3
    for i=1:sz2 
      for j=1:sz1
        tmpv(j, i, m) = j;          % <- nominal (9 x 30 x 60)
        tmpx(j, i, m) = i;
        tmpa(j, i, m) = m;
      end
    end
  end
  tmpx(:,L1a_err) = NaN;
  tmpz(:,L1a_err) = NaN;
  tmpv(:,L1a_err) = NaN;
  cxtrak   = cat(3,cxtrak, tmpx);
  catrak   = cat(3,catrak, tmpa);
  cfov     = cat(3,cfov,   tmpv);
  
% Get granule ID, should return one value per aTrack (60), but deal with xxxxxxxxxxxxx
  clear granID junk B C
  k = 1;
  for i=1:length(geo.MidTime)  % length(geo.Granule_ID)
    %if(strcmp(geo.Granule_ID(i,1:3),'NPP'))
      %junk(k) = str2num(geo.Granule_ID(i,4:end));
      %junk(k) = str2num(cFnams{fn}(98:100));
      junk(k) = str2num(cell2mat(regexp(cFnams{fn}, '(?<=_d)\d+(?=_t)','match')));
      %junk(k) = str2num(cell2mat(regexp(cFnams{fn}, '(?<=_g)\d+(?=_v20a)','match')));
      k = k + 1;
    %end
  end
% Expand to full size - take first valid value.
  B = repmat(junk(1),sz1,sz3,sz2); 
  C = permute(B,[1 3 2]);
  C(:,L1a_err) = NaN; 
  cgran = cat(3,cgran, C);

  errFlg = reshape(L1a_err,[],1);                     % (30 x nscan) 1: bad, 0: good.
   junk  = geo.FORTime;   junk(L1a_err) = NaN;
  tmpt   = iet2dnum(junk);    clear junk;
  % Need to expand ctim to match every FOV so that sub-setting works
  junk   = repmat(tmpt,[1 1 9]);
  tmpt   = permute(junk,[3 1 2]);    clear junk;
%  ctim   = [ctim; reshape(junk,[],1)]; clear junk;           % <- (16200 x 1) all FOVs
  ctim   = cat(3,ctim, tmpt);
    junk = geo.Latitude;              junk(:,L1a_err) = NaN;
  clat   = cat(3,clat, junk);
    junk = geo.Longitude;             junk(:,L1a_err) = NaN;
  clon   = cat(3,clon, junk);
    junk = geo.SolarZenithAngle;      junk(:,L1a_err) = NaN;
  csozn  = cat(3,csozn, junk);
    junk = geo.SatelliteZenithAngle;  junk(:,L1a_err) = NaN;
  csazn  = cat(3,csazn, junk);
  
%  fprintf('.');
end
fprintf('\n');
whos clat clon ctim cxtrak catrak cfov cgran csozn csazn

% Subset to center tracks 15,16 in preparation for SNO processing.
  cCnLat  = clat(:,15:16,:);
  cCnLon  = clon(:,15:16,:);
  cCnTim  = ctim(:,15:16,:);
  cCnFov  = cfov(:,15:16,:);
  cCnAtr  = catrak(:,15:16,:);
  cCnXtr  = cxtrak(:,15:16,:);
  cCnGid  = cgran(:,15:16,:);
  cCnSozn = csozn(:,15:16,:);
  cCnSazn = csazn(:,15:16,:);
  ncCn    = numel(cCnLat);
  disp(['Found ' num2str(ncCn) ' CrIS center track FOVs'])

% *******************   Get CHIRP granule files  *******************
fprintf('Loading CHIRP geo\n'); 
xDir = strcat('/asl/hpcnfs1/chirp/airs_L1c/',cYr,'/',sprintf('%03d',jday),'/');
xLst = dir(sprintf('%sCHIRP_AIRS-L1C_d%4d%02d%02d*.nc',xDir,nYr1,nMn1,nDy1));
disp(['Found ' num2str(numel(xLst)) ' CHIRP Files']);

if(numel(xLst) < 1) error('Insufficient CHIRP data'); return; end

alat  = []; alon = []; atim = []; aatrak = []; axtrak = []; asolzn = []; agran = [];
alnfr = []; aAscFlg = [];
for fn = 1:numel(xLst);
  chrp = read_netcdf_lls(strcat(xLst(fn).folder,'/',xLst(fn).name));
  alat     = [alat;   chrp.lat];
  alon     = [alon;   chrp.lon]; 
  atim     = [atim;   airs2dnum(chrp.obs_time_tai93)];               % convert to matlab time.
  aatrak   = [aatrak; chrp.atrack_ind];
  axtrak   = [axtrak; chrp.xtrack_ind];
  asolzn   = [asolzn; chrp.sol_zen];
%  agran    = [agran;  l1c.gindex];
  alnfr    = [alnfr;  chrp.land_frac];
  aAscFlg  = [aAscFlg; chrp.asc_flag];
  fprintf('.');
end

fprintf('\n');
whos alat alon atim aatrak axtrak asolzn agran alnfr;

icnt = find(axtrak == 43 | axtrak == 44 | axtrak == 45 | axtrak == 46 | ...
            axtrak == 47 | axtrak == 48);
aCnTim  = atim(icnt)';            % NB ': to match CrIS center arrays
aCnLat  = alat(icnt)';
aCnLon  = alon(icnt)';
aCnatr  = aatrak(icnt)';
aCnxtr  = axtrak(icnt)';
%aCnGrn  = agran(icnt)';
aCnSozn = asolzn(icnt)';
aCnLnfr = alnfr(icnt)';
aCnAsc  = aAscFlg(icnt)';
disp(['Found ' num2str(numel(aCnTim)) ' CHIRP center track FOVs'])

% Get ascending flag (1=asc, 0=desc)
%lat45   = alat(:,45);  clear junk;
%for i=1:length(lat45)-1 junk(i+1) = lat45(i+1)-lat45(i); end
%  junk(length(lat45)) = junk(i);
%aAsc = junk;
%aAsc(aAsc<=0) = 0;
%aAsc(aAsc>0)  = 1;
%clear junk;
% Expand to match near-nadir arrays:
%aAsc = repmat(aAsc,6,1);


%  ****************** SECTION 3: Get indexes of SNOs  ********************
%            pre-compute position vectors - saves a heap of time *********
%            NaNs are conserved
fprintf('computing position vectors\n');
P1 = zeros(numel(aCnLat),3);                % CHIRP position
P2 = zeros(numel(cCnLat),3);                % CRIS position
for ii = 1:numel(aCnLat)                                   % Number of Swath.
    P1(ii,:) = [ cos(aCnLat(ii)*pi/180.0) * cos(aCnLon(ii)*pi/180.0), ...
         cos(aCnLat(ii)*pi/180.0)*sin(aCnLon(ii)*pi/180.0), sin(aCnLat(ii)*pi/180.0) ];
end
for ii = 1:numel(cCnLat)
    P2(ii,:) = [ cos(cCnLat(ii)*pi/180.0) * cos(cCnLon(ii)*pi/180.0), ...
         cos(cCnLat(ii)*pi/180.0)*sin(cCnLon(ii)*pi/180.0), sin(cCnLat(ii)*pi/180.0) ];
end    
whos P1 P2;

% Faster Algorithm 

dta1  = 0.0225/86400;      % <- AIRS time between adjacent samples within Centre group
dta2  = 2.555/86400;       % <- AIRS time between each centre FOV group
dtc1  = 0.200/86400;       % <- CRIS time between adjacent FORs within center group
dtc2  = 7.800/86400;       % <- CRIS time between each center FOR group
wndta = fix(6*maxDtim/dta2)+1;    % <- no. AIRS samples in time window set by maxDtim.

arsa  = find(aCnTim > cCnTim(1),1) + wndta;    % sample to start AIRS
aren  = arsa + 2*wndta;

%{
smPos = [;];  smTd = []; smPhi = 75.0; 
for jj = 5:18:length(cCnLat)
  for ii = arsa:6:length(aCnLat)
    cnPhi = real(acos( sum(P1(ii,:).*P2(jj,:)) )*180/pi);
    if(cnPhi < smPhi) 
      smPhi   = [smPhi,cnPhi]; 
      smPos   = [smPos;[ii,jj]];                    % pos(ii: AIRS index. jj: CRIS)
      smTd    = [smTd,(aCnTim(ii) - cCnTim(jj))];
    end
  end
  if(~mod(jj,365)) fprintf('.'); end
end
%} 
% save('/asl/s1/chepplew/projects/sno/airs_cris/HR/20130827_snapshot.mat');

%        Compute separations and save indexes when criteria are met
fprintf('Computing separations\n');
dist = 0.0; m = 0; k = 1;
pos = []; tdiff = [];
tic
for jj = 1:1:numel(cCnLat)
  %arst = max(1,find(aCnTim > cCnTim(jj),1) - wndta);
  %arsp = min(arst + 2*wndta, length(aCnLat));
  for ii = 1:1:numel(aCnLat) % arst:1:arsp
    Dtim = abs(aCnTim(ii) - cCnTim(jj));
    if Dtim <= maxDtim       
      m = m+1;
      % dist = distance(arCnLat(ja),arCnLon(ja),iaCnLat(ji),iaCnLon(ji)); % way too slow
      Phi = real(acos( sum(P1(ii,:).*P2(jj,:)) )*180/pi);
      if Phi <= maxDphi                % phi = 0.18 deg (0.0031 rad) => 20 km.
        dist(k) = Phi;                % dist;
        pos     = [pos;[ii,jj]];      % pos(ii: AIRS index. jj: CRIS)
        tdiff   = [tdiff,(aCnTim(ii) - cCnTim(jj))];
        k = k+1;
      end
    end
  end
  if(~mod(jj,1000)) fprintf('.'); end
end
toc
dist = real(dist);                         % can get complex phi.
fprintf('\n');

% save(['/home/chepplew/data/sno/airs_cris/' strRes '/' strYr '/' strYr strMn strDy '_postTest.mat']);

fprintf(1,'Number matches found %d\n',size(pos,1));

% This gives us duplicate hits - so need to select only unique pairs.
%  this is done by sequentially testing for uniqueness on each sensor.
upos = [];
if(size(pos,1) > 10 )
  un1 = [];  un2 = [];  upos = [];
  [x,ib,ix] = unique(pos(:,1)); un1  = [x,pos(ib,2)];
  [x,ib,ix] = unique(un1(:,2)); upos = [un1(ib,1),x]; clear un1;
  fprintf('Number unique SNOs %d %d\n',size(upos));
  nSNO = size(upos,1);
end
if(size(upos,1) > 2 )
% Record time and space separations (pre reloading of files)
  sno = struct;
  sno.aTime     = aCnTim(upos(:,1));    sno.cTime   = cCnTim(upos(:,2));
  sno.aLat      = aCnLat(upos(:,1));    sno.cLat    = cCnLat(upos(:,2));
  sno.aLon      = aCnLon(upos(:,1));    sno.cLon    = cCnLon(upos(:,2));
  sno.aAtrack   = aCnatr(upos(:,1));    sno.cAtrack = cCnAtr(upos(:,2));
  sno.aXtrack   = aCnxtr(upos(:,1));    sno.cXtrack = cCnXtr(upos(:,2));
%  sno.aGran     = aCnGrn(upos(:,1));    sno.cGran   = cCnGid(upos(:,2));
  sno.aAsc      = aCnAsc(upos(:,1));     sno.cFov    = cCnFov(upos(:,2));
  sno.aSolzen   = aCnSozn(upos(:,1));   sno.cSolzen = cCnSozn(upos(:,2));
  sno.alnfrac   = aCnLnfr(upos(:,1));   sno.clnfrac = [];
  sno.tdiff     = sno.aTime - sno.cTime;    
%
  uP1    = P1(upos(:,1),:);      uP2 = P2(upos(:,2),:);    % P1: AIRS, P2: CRIS.
  for i = 1:size(upos,1)
     uPhi(i)   = real( acos(sum(uP1(i,:).*uP2(i,:))) * 180/pi );
  end
  sno.dist = uPhi';

% ---------------------------------------------------------------------
% **** Section 4 Reload CHIRP and CRIS files and save only required Obs.
% ---------------------------------------------------------------------
  fprintf('loading CHIRP obs\n');
  xCntr = [43:48];
  ps = 1; pe = 0; aCount = 0;  k = 1;
  tarad  = [;]; taAtrk = [];  taXtrk = [];  tasozn = []; tasazn  = []; 
  talnfr = [];   tatim = [];  talat  = [];  talon  = []; taPrc   = [];
  taReas = [];   taRaQc = [];  taSyQc = [];
  nbr_ra  = zeros(nSNO,2,2645);
  for fn = 1:numel(xLst); % fn = 1;
    chrp     = read_netcdf_lls(strcat(xDir,xLst(fn).name));
    atim     = airs2dnum(chrp.obs_time_tai93);     % conv to matlab time.
    aatrak   = chrp.atrack_ind;
    axtrak   = chrp.xtrack_ind;
    tCnInd   = find(axtrak == 43 | axtrak == 44 | axtrak == 45 | ...
                    axtrak == 46 | axtrak == 47 | axtrak == 48); 
    pe       = ps + numel(tCnInd) - 1;
    %  tests which data points we want to keep match to appropriate FOR
    fprintf('%d: %d %d ,',fn,ps,pe);
    junk   = ismember(upos(:,1),[ps:pe]);       aSams = upos(junk,1);    clear junk;
    if ( length(aSams >= 1)  )
      fprintf('%d \t',length(aSams));
      [nx ny] = size(chrp.rad);
      % subset onto center track b4 using upos index.
      junk    = chrp.rad(:,[tCnInd]);
      % aCnRad  = permute(junk,[3,2,1]);            clear junk;
      tarad   = [tarad,  junk(:,aSams-ps+1)];
      junk    = chrp.rad_qc(tCnInd);
      % aCnPrc  = permute(junk,[3,2,1]);
      taRaQc   = [taRaQc; junk(aSams-ps+1)];    clear junk;
      %junk    = chrp.syn_qc(tCnInd);
      % aCnReas = permute(junk,[3,2,1]);
      %taSyQc  = [taSyQc; junk(aSams-ps+1)];    clear junk;
      taSyQc  = [taSyQc; chrp.syn_qc];
      
      %for j = 1:length(aSams)
      %  xxt = aSams(j)-ps+1;
      %  nbr_ra(k,1,:) = chrp.rad(:,aCnatr(xxt),aCnxtr(xxt)-1);
      %  nbr_ra(k,2,:) = chrp.rad(:,aCnatr(xxt),aCnxtr(xxt)+1);
      %  k = k+1;
      %end
    end
    ps = pe+1;
    aCount = aCount + length(aSams);  fprintf('%d \n',aCount);
%   fprintf('.');
  end
  fprintf('\n');
  % nbr_ra = permute(nbr_ra, [3,2,1]);              % match tarad.
  % whos tarad nbr_ra taReas taPrc;

% ************  Re-load CRIS SNO Obs data ************

  fprintf('Loading CIRS SNO Obs\n');
  ps = 1; pe = 0; cCount = 0;   k = 1;       tcrLW = [];   tcrMW = [];  tcrSW = [];
  tcatrk = [];  tcxtrk = [];    tcFov = [];
  tcTim  = [];   tcLat = [];    tcLon = [];  tcSozn = []; tcSazn = []; tcLnfr = [];
  snVara = [;];  cnObs = [;]; nbr_rLW = zeros(nSNO,4,npLW);
  for fn = 1:numel(crLst)
    clear g2 rLW2 rLWp1;
    load(cFnams{fn});
% Check granule and generate FOV, a-track, x-track indexes - only works with fixed array size
    if(ndims(geo.Latitude) ~= 3) fprintf('ERROR: wrong ndims of geo.Lat\n'); end
    [sz1 sz2 sz3] = size(geo.Latitude);

    rLW(:,:,L1a_err) = NaN;  
    cradLW = rLW(:,:,15:16,:);          % <- (717 x 9 x 2 x 60) per gran or fewer
    rMW(:,:,L1a_err) = NaN;
    cradMW = rMW(:,:,15:16,:);          % <- (869 x 9 x 2 x 60) per gran
    rSW(:,:,L1a_err) = NaN;
    cradSW = rSW(:,:,15:16,:);          % <- (637 x 9 x 2 x 60) per gran
 % grab the next granule and append the first a-track to current granule to create rLWp1 etc. 
      g2        = load(cFnams{fn+1});
      rLW2      = g2.rLW;
      rLW2(:,:,g2.L1a_err) = NaN;
      rLW2_atr1 = rLW2(:,:,:,1);
      rLWp1     = cat(4,rLW, rLW2_atr1);

    
    ncObs  = sz1 * 2 * sz3;              % 9x2x60 = 1080
    pe     = ps + ncObs -1;
    junk   = ismember(upos(:,2),[ps:pe]);
    cSams  = upos(junk,2);                clear junk;
    fprintf('%d: %d %d ,',fn,ps,pe);

    if ( length(cSams) >= 1 )
      fprintf(1,'\t %d ',length(cSams));
      tcrLW   = [tcrLW,  cradLW(:,cSams-ps+1)];
      tcrMW   = [tcrMW,  cradMW(:,cSams-ps+1)];
      tcrSW   = [tcrSW,  cradSW(:,cSams-ps+1)];
      for j = 1:length(cSams)
        clear trLW;
        % fold at granule boundary edges if needed
        map   = cFOR( cCnXtr(cSams(j)) ).d{cCnFov(cSams(j))};
        trLW(1,:) = rLWp1(:,map(1,3),map(1,2),cCnAtr(cSams(j))+map(1,1));
        trLW(2,:) = rLWp1(:,map(2,3),map(2,2),cCnAtr(cSams(j))+map(2,1));
        trLW(3,:) = rLWp1(:,map(3,3),map(3,2),cCnAtr(cSams(j))+map(3,1));
        if(cCnAtr(cSams(j)) > 1 | map(4,1) >= 0 )
         trLW(4,:) = rLWp1(:,map(4,3),map(4,2),cCnAtr(cSams(j))+map(4,1));
        else
         trLW(4,:) = rLWp1(:,map(4,3),map(4,2),cCnAtr(cSams(j)) );
        end
      % record time stamps for each neighbour for validation
        % TBD
        nbr_rLW(k,:,:) = trLW;
        k = k+1;
      end         % end for cSams
    end           % end: if cSams
    ps = pe + 1;
    cCount = cCount + length(cSams); fprintf('%d \n',cCount);
  end             % for fn loop
  nbr_rLW = permute(nbr_rLW,[3,2,1]);
  
  rx    = [tcrLW; tcrMW; tcrSW];
  fx    = [vLW(:); vMW(:); vSW(:)];

  % whos tcr* nbr_rLW rc fc

% ----------------------------------------------------------------
%          Convert CrIS hires to CrIS mid (CHIRP) 
% ----------------------------------------------------------------
if( all(strcmp(cris_res, {'HIGH','MID'})) )
  disp(['Converting CrIS SDR from FSR to MSR'])

  wlaser = 773.1301;  % nominal wlaser
  bstr   = {'LW' 'MW' 'SW'};
  for i = 2:3
    [instX, userX]   = inst_params(bstr{i}, wlaser, opt1);
    iix              = [find(fx == userX.v1) find(fx == userX.v2)];
    X(i).v           = fx(iix(1):iix(2));
    X(i).r           = rx(iix(1):iix(2),:);

    % interpolate to the new user grid
    [Y(i).r, Y(i).v] = finterp(X(i).r, X(i).v, userX.dv);
  end 
  rc  = [tcrLW; Y(2).r; Y(3).r];
  fc  = [vLW;   Y(2).v; Y(3).v];
elseif( strcmp(cris_res{1}, cris_res{2}) )
  disp(['Retaining original CrIS SDR spectral resolution'])
  rc = rx;
  fc = fx;
  clear rx fx;
end

%btx = real(rad2bt(fc, rc));

%{
% ----------------------------------------------------------------
%  Convert AIRS to CrIS 
% ----------------------------------------------------------------
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

load('/home/chepplew/projects/cris/cris_freq_2grd.mat'); fcris = vchan;   % fcris [1317x1]
% get AIRS frequency grid 
load('/home/chepplew/projects/airs/airs_f.mat');  % f [2378] fairs [2645] AIRS frequency grids
      f2645 = fairs;                     % the filled 2645 AIRS frequency grid
f_chirp = chrp.wnum;

kn = 1000;

if(exist('opt1','var')) 
  switch opt1.user_res
    case 'lowres', opt1.cfile = 'corr_lowres.mat';
    case 'midres', opt1.cfile = 'corr_midres.mat';
    case 'hires',  opt1.cfile = 'corr_hires.mat';
  end  
  disp('Converting AIRS to CrIS');
  [ns sz] = size(tarad)
  if(ns ~= 1483 & sz == 1483) xra = tarad'; sz=ns; else xra = tarad; end
  ra2c = [];
  for j = 1:kn:sz
    sx = j:min(j+kn-1,sz);
    [rrx fa2c] = airs2cris(xra(:,sx),f2645,sfile,opt1);
    ra2c       = [ra2c, rrx];
    fprintf(1,'.');
  end

disp('Completed conversion airs2cris');
end
%}
%%%%%%%%%%%%%%%%%%%% Saving data %%%%%%%%%%%%%%%%%%%%%
  clear atime alat alon asolzen asatzen ra;
  clear ctime clat clon cfov csolzen csatzen tdiff dist;
  r_chirp    = tarad;  f_chirp = chrp.wnum;  %ichans = l1c.chanID; 
  radQc  = taRaQc;
  synQc  = taSyQc;
  savFN   = strcat(savDir, cYr, ...
            '/sno_chirp_cris_asl_',cYr,cMn,cDy,'_',vers,'.mat');
  savVars = {'sno','r_chirp','f_chirp','rc','fc','upos',...
             'maxDtim','maxDphi','nbr_rLW','nbr_ra',...
	     'rdate','cris_res','xDir','xLst','cFnams','radQc',...
	     'synQc','opt1'};
  fprintf(1,'Saving data to file: %s\n',savFN);
  save(savFN, savVars{:});


end

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Plotting Section   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /asl/matlib/aslutil
addpath /asl/matlib/plotutils
addpath /asl/packages/airs_decon/source
% distribution maps

figure(1);clf;simplemap(sno.cLat, sno.cLon, sno.dist);
figure(1);clf;simplemap(sno.cLat, sno.cLon, sno.tdiff*24*60);

% separation plots
figure(2);clf;plot(sno.cLat - sno.aLat,'.');ylabel('Latitude difference (deg)')
figure(2);clf;plot(sno.cLon - sno.aLon,'.');ylabel('Longitude difference (deg)');
figure(2);scatter(sno.dist, sno.tdiff*24*60,[]);xlabel('g.c. angular dist');
   ylabel('time diff (mins)');

% BT plots
% CrIS
[sz1 sz2] = size(rc);
cbt     = real(rad2bt(fc, rc));
cbt_ham = single(hamm_app(double(cbt)));
cbtm    = nanmean(cbt_ham,2); 
cbtsd   = nanstd(cbt_ham,0,2);
figure(3);clf;subplot(2,1,1);plot(fc, cbtm,'-'); grid on;xlim([640 2550]);
   subplot(2,1,2);plot(fc,cbtsd,'-'); grid on; xlim([640 2550]);

figure(3);clf;plot(cbt_ham(403,:),'.');

% AIRS
abt   = real(rad2bt(f, ra));
abtm  = nanmean(abt,2);
abtsd = nanstd(abt,0,2);
figure(3);clf;subplot(2,1,1);plot(f,abtm,'-');grid on;axis([640 1200 190 270]);
   subplot(2,1,2);plot(f,abtsd,'-');grid on;axis([640 1200 0 20]);
figure(3);clf;plot(abt(759,:),'.');
 
% together
nSam = size(ra,2);   
figure(3);clf;plot([1:nSam], abt(759,:),'.', [1:nSam], cbt_ham(403,:),'o');
figure(3);clf;plot([1:nSam], abt(759,:) - cbt_ham(403,:),'.');
dbtbins = [-15:0.2:15];  dbtcens = [-14.9:0.2:14.9];
bias_pdf = histcounts(abt(759,:) - cbt_ham(403,:), dbtbins);
figure(3);plot(dbtcens, bias_pdf,'.-');grid on;

figure(3);plot(fc,cbtm,'-',fa,abtm,'-');axis([640 1100 200 270]);grid on;
dbt900 = abt(759,:) - cbt_ham(403,:);
figure(1);clf;simplemap(sno.cLat, sno.cLon, dbt900');

% Cris uniformity
new_rLW  = permute(nbr_LW,[3 2 1]);         % 1st dim must be channel
nbr_btLW = rad2bt(vLW, new_rLW);            % [717 x 4 x nSNO]
dxm_btLW = max(nbr_btLW,[],2) - min(nbr_btLW,[],2);
mn_btLW  = nanmean(nbr_btLW,2);

figure(2);clf;simplemap(pre.cSnoLat, pre.cSnoLon, squeeze(mn_btLW(403,1,:)) );
  title('2016.01.01 AC SNO CrIS 900wn 4-neighbor mean BT')
figure(1);clf;simplemap(pre.cSnoLat, pre.cSnoLon, squeeze(cbt_ham(403,:))' );
  title('2016.01.01 AC SNO CrIS 900wn BT');
figure(3);clf;simplemap(pre.cSnoLat, pre.cSnoLon, squeeze(dxm_btLW(403,1,:)) );
  title('2016.01.01 AC SNO CrIS 900wn 4-neighbor max-min BT');

figure(1);clf;scatter(dbt900, squeeze(dxm_btLW(403,1,:)), [], cbt_ham(403,:) ); grid on;

%}
