function [] = make_IASI_CRIS_SNO_frmL1(rdate,cris_res,iasi,cris,opt1,vers);
%
% function [] = make_IASI_CRIS_SNO_frmL1(date,cris_res,iasi,cris,opts,vers)
%
% produce files of SNO using CRIS from CCAST processing at UMBC.ASL and IASI.
%
% INPUT PARAMETERS:
%   1.    rdate:    string 'YYYY/MM/DD' to process
%   2.    cris_res: two strings {'high','mid','low'}{ditto}
%                   1st: the CCAST SDR source data.
%                   2nd: SNO converted CrIS resolution.
%   3.    iasi:     first, 2nd or 3rd (MetOp-A or B, C) [1,2,3]
%   4.    cris:     first or second (NPP, JPSS-1, JPSS-2)  [1,2,3]
%   5.    opt1:     structure of options to control iasi2cris conversion
%                   if entered with no fields, no conversion takes place.
%   6.    vers:     for special tests or revision of CrIS data.
%                   current options: {'', 'noaa_poff','noaa_pon'}
%
%1. Sources & Dependencies:
%   /asl/data/cris/ccast/sdr60_hr/YYYY/JJJ/SDR_d20130828_t0006589.mat
%   /asl/data/IASI/L1C/YYYY/MM/DD/IASI_xxx_1C_M01_20131020005656Z_20131020005959Z
%   
%2. Destination & Outputs:
%   /asl/s1/chepplew/projects/sno/iasi_cris/HR/ or /LR/
%
% Notes:  HM: L1a_err, in the SDR files is a 30 x nscan array, one value for each FOR:
%         1 = bad, 0 = OK.   
%
% 2-Jan-2018 CLH, added iasi and cris mission numbers
% Jul 2020. CLH: New load functions.

%cd /home/chepplew/projects/sno/makeSNO/;

addpath /asl/matlib/time                              % tai2dnum
addpath /home/chepplew/gitLib/asl_sno/sno_maker
addpath /home/chepplew/gitLib/asl_sno/utilities
addpath /asl/matlib/aslutil                           % simplemap, findfiles
%addpath /asl/rtp_prod/iasi/readers                   % readl1c_epsflip_all
%addpath /home/chepplew/myLib/matlib/aslutil          % utc2tai2000 used by readl1c_epsflip_all
%addpath /home/chepplew/myLib/matlib                  % tai2utc1958.m
%
addpath /asl/packages/iasi_decon                      % iasi2cris
addpath /asl/packages/ccast/source                    % inst_params
addpath /asl/packages/airs_decon/source               % hamm_app

% Criteria for separation and delay
maxDtim   = 0.014;                  % day. 600/86400 secs =  Mattime is decimal day (20 mins=0.0139)
maxDphi   = 0.18;                   % deg. 0.07 = 7.8 km Earth radius: 6371 km. 1 deg = 111 km.

% ------------------------------------------------------
%   check and validate input parameters
% ------------------------------------------------------

% Check all four input parameters are entered.
if(nargin ~= 6) error('Please enter all 6 input arguments'); return; end

% Check date string
whos rdate; disp(rdate); fprintf('\n');
try 
   D = datenum(rdate,'YYYY/MM/DD');
catch
   error('Incorrect Date Format')
   return
end

% Check CrIS CCAST source and SNO converted resolution
cris_res = upper(cris_res);
all_res  = {'HIGH','MID','LOW'};
if(~ismember(cris_res,all_res))
  error('Invalid CrIS resolution. Choices are: high or low');
  return
end

% check version reference
vers = lower(vers);
if(length(vers) > 9) error('Version reference too long, 9 chars max'); return; end
allvers = {'','noaa_pon','noaa_poff','v20a','v20d','nasa'};
if(~ismember(vers,allvers)) error('Not valid CrIS L1C version');return; end

% Check IASI mission
if(~ismember(iasi,[1,2])) error('IASI can take only values 1 or 2'); return; end
if(iasi == 1); IX = ''; end
if(iasi == 2); IX = '2'; end

% Check CrIS mission (NPP=1, J-1=2)
if(~ismember(cris,[1,2])) error('CrIS can only take values 1 or 2'); return; end
if(cris == 1); CX = ''; end
if(cris == 2); CX = '2'; end

% Check conversion options, names and values
%opts  = struct;
opt1
names = fieldnames(opt1);
if(isempty(names)) clear opt1; end 
if(~all(ismember(names,{'hapod','user_res','nguard'})))
  error('Not all options for iasi2cris are present, plz chceck'); return;
end

% --------------------------------------------------------
%    Initialize and setup
% --------------------------------------------------------

%  convert date to required formats
strYr   = rdate(1:4);   strMn = rdate(6:7);    strDy = rdate(9:10);
numYr   = str2num(strYr);  numMn = str2num(strMn);   numDy = str2num(strDy);
junk    = sprintf('%4d/%02d/%02d',numYr-1,12,31);
jday    = datenum(rdate)-datenum(junk);  clear junk;           % needed for CRIS directory

% Create output directory and channel numbers per band (TBC)   
% slelect source {1} and destination directories {2}.
switch cris_res{1}
  case 'HIGH'
    CR = 'HR';
    npLW = 717;  npMW = 869;   npSW = 637;
  case 'MID'
    CR = 'MR';
    npLW = 717;  npMW = 865;   npSW = 633;
  case 'LOW'
    CR = 'LR';
    npLW = 717;  npMW = 437;   npSW = 163;
end
switch cris_res{2}
  case 'HIGH'
    savDir = ['/asl/s1/chepplew/data/sno/iasi' IX '_cris' CX '/ASL/HR/' strYr '/'];
    if(strcmp(vers,'noaa_pon'))
       savDir = ['/asl/xfs3/sno/iasi' IX '_cris' CX '/NOAA_pon/HR/' strYr '/'];
    end
    if(strcmp(vers,'noaa_poff'))
       savDir = ['/asl/xfs3/sno/iasi' IX '_cris' CX '/NOAA_poff/HR/' strYr '/'];
    end
  case 'MID'
    savDir = ['/asl/s1/chepplew/data/sno/iasi' IX '_cris' CX '/ASL/MR/' strYr '/'];
  case 'LOW'
    savDir = ['/asl/s1/chepplew/data/sno/iasi' IX '_cris' CX '/ASL/LR/' strYr '/'];
end

disp(['CrIS: ' num2str(cris) ' IASI: ' num2str(iasi) ' Res: ' cris_res ' Vers: ' vers]);

% CCAST granules for NPP are in different locations & names after 2017
%if (numYr <= 2017) ccast = 'sdr60_'; pattern='SDR_*d'; end
%if (numYr >= 2018) ccast = 'sdr45_'; pattern='CrIS_SDR_*d'; end
ccast = 'sdr45_'; pattern='CrIS_SDR_*d';

% ************ Get CRIS granule files list *****************
if(strcmp(cris_res{1},'HIGH') & cris == 1) 
  crDir = ['/asl/cris/ccast/' ccast 'npp_HR/'];
  crDir = sprintf('%s%s/%03d/',crDir, strYr,jday);
  crLst = dir(strcat(crDir, pattern, strYr, strMn, strDy,'_t*.mat'));
end
if(strcmp(cris_res{1},'HIGH') & cris == 2) 
  crDir = '/asl/cris/ccast/sdr45_j01_HR/';
  crDir = sprintf('%s%s/%03d/',crDir, strYr,jday);
  crLst = dir(strcat(crDir, pattern, strYr, strMn, strDy,'_t*.mat'));
  if(strcmp(vers,'noaa_pon'))
    crDir  =  '/home/chepplew/data/cris/noaa/sdr_j01_hr/';
    crDir  = sprintf('%s%s/%03d/',crDir, strYr,jday);
    crLst  = dir(strcat(crDir, 'CrIS_SDR_j01_s45_d', strYr, strMn, strDy,'_tXXX_noaa_PLon_*.mat'));
  end
  if(strcmp(vers,'noaa_poff'))
    crDir  =  '/home/chepplew/data/cris/noaa/sdr_j01_hr/';
    crDir  = sprintf('%s%s/%03d/',crDir, strYr,jday);
    crLst  = dir(strcat(crDir, 'CrIS_SDR_j01_s45_d', strYr, strMn, strDy,'_tXXX_noaa_PLoff_*.mat'));
  end
end
if(strcmp(cris_res{1},'LOW') & cris == 1)  
  crDir = ['/asl/data/cris/ccast/' ccast 'npp_LR/'];      
  crDir = sprintf('%s%s/%03d/',crDir, strYr,jday);
  crLst = dir(strcat(crDir, pattern, strYr, strMn, strDy,'_t*.mat'));
end
if(strcmp(cris_res{1},'LOW') & cris == 2) 
  crDir = '/asl/data/cris/ccast/sdr45_j01_LR/';  
  crDir = sprintf('%s%s/%03d/',crDir, strYr,jday);
  crLst = dir(strcat(crDir, pattern, strYr, strMn, strDy,'_t*.mat'));
end
if(strcmp(vers,'nasa') & cris == 1)
  crDir = '/umbc/isilon/rs/strow/asl/cris/nasa1b/';
  crDir = sprintf('%s%s/%03d/', crDir, strYr, jday);
  crLst = dir(strcat(crDir, 'SNDR.SNPP.CRIS.*.nc'));
end
if(strcmp(vers,'nasa') & cris == 2)
  crDir = '/umbc/isilon/rs/strow/asl/n20/cris/nasa1b/';
  crDir = sprintf('%s%s/%03d/', crDir, strYr, jday);
  crLst = dir(strcat(crDir, 'SNDR.J1.CRIS.*.nc'));
end

if(numel(crLst) < 1) error('Insufficient CrIS L1C data'); return; end

% Collect all granule filenames:
for i=1:numel(crLst)
  cFnams{i} = strcat(crLst(i).folder, '/', crLst(i).name);
end

disp(['Found ' num2str(numel(crLst)) ' CrIS CCAST L1c Files in ' crDir]);

% ****************  Get IASI granule files list *******************
if(iasi == 1) 
  iaDir = strcat('/asl/data/IASI/L1C/',rdate,'/');
  iaLst = dir(strcat(iaDir,'IASI_xxx_1C_M02*'));
end
if(iasi == 2) 
  iaDir = strcat('/asl/data/IASI/L1C/',rdate,'/');
  iaLst = dir(strcat(iaDir,'IASI_xxx_1C_M01*'));
end
disp(['Found ' num2str(numel(iaLst)) ' IASI L1C files in ' iaDir]);

if(numel(iaLst) < 1) error('Insufficient IASI L1C data'); return; end

% ------------ Load CrIS geo data --------------
fprintf('loading CrIS geo data\n');

cgeo = load_cris_geo_for_sno(cFnams,vers);

% Subset to center tracks 15,16 in preparation for SNO processing.
cCn.Lat  = cgeo.lat(:,15:16,:);
cCn.Lon  = cgeo.lon(:,15:16,:);
cCn.Tim  = cgeo.tim(:,15:16,:);
cCn.Fov  = cgeo.fov(:,15:16,:);
cCn.Atr  = cgeo.atrak(:,15:16,:);
cCn.Xtr  = cgeo.xtrak(:,15:16,:);
cCn.Gid  = cgeo.gran(:,15:16,:);
cCn.Sozn = cgeo.sozn(:,15:16,:);
cCn.Sazn = cgeo.sazn(:,15:16,:);
ncCn    = numel(cCn.Lat);
disp(['Found ' num2str(ncCn) ' CrIS center track FOVs'])

whos cCn*

% ------------ Load IASI geo data --------------
fprintf('loading IASI geo data\n');

igeo = load_iasi_geo_for_sno(iaLst);

% Subset center track FOVs
icntr = find(igeo.xtrak(:,1) == 15 | igeo.xtrak(:,1) == 16);
disp(['Found ' num2str(numel(icntr)) ' IASI center track FOVs'])
iCn.Tim   = igeo.tim(icntr,:);                           % <- [N x 1] arrays
iCn.Lat   = igeo.lat(icntr,:);
iCn.Lon   = igeo.lon(icntr,:);
iCn.Fov   = igeo.fov(icntr,:);
iCn.Qual  = igeo.qual(icntr,:);
% iCn.Scln  = igeo.scln(icntr,:);
iCn.Solz  = igeo.solzen(icntr,:);
iCn.Satz  = igeo.satzen(icntr,:);
whos iCn*


% --------------------------------------------------------------
%        pre-compute position vectors - saves a heap of time
% --------------------------------------------------------------
fprintf('computing position vectors\n');
P1 = zeros(numel(cCn.Lat),3);                              % CRIS
P2 = zeros(numel(iCn.Lat),3);                              % IASI
for ic = 1:numel(cCn.Lat)               % Number of centres.
    P1(ic,:) = [cos(cCn.Lat(ic)*pi/180.0) * cos(cCn.Lon(ic)*pi/180.0), ...
         cos(cCn.Lat(ic)*pi/180.0)*sin(cCn.Lon(ic)*pi/180.0), sin(cCn.Lat(ic)*pi/180.0) ];
end
for ii = 1:numel(iCn.Lat)
    P2(ii,:) = [cos(iCn.Lat(ii)*pi/180.0) * cos(iCn.Lon(ii)*pi/180.0), ...
         cos(iCn.Lat(ii)*pi/180.0)*sin(iCn.Lon(ii)*pi/180.0), sin(iCn.Lat(ii)*pi/180.0) ];
end    
%
% save(strcat(savDir,'sno_Iasi_hrCris_ln145_',strYr, strMn, strDy,'.mat'));

%{
% test windows for nearest samples (not used)
dti1  = 2.5347e-06;        % <- (day) IASI time between adjacent FORs within Centre group
dti2  = 9.2593e-05;        % <- IASI time between each centre FOV group
dtc1  = 0.200/86400;       % <- CRIS time between adjacent FORs within center group
dtc2  = 7.800/86400;       % <- CRIS time between each center FOR group
wndta = fix(8*maxDtim/dti2)+1;    % <- no. IASI samples in time window set by maxDtim.

irsa  = max(1,find(iCnTim > cCnTim(1),1) - wndta);    % sample to start IASI
iren  = irsa + 2*wndta;

fprintf('Finding nearest approach\n');
smPos = [;];  smTd = []; smPhi = 75.0; 
for jj = 5:18:length(cCnLat)
  for ii = irsa:8:length(iCnLat)
    cnPhi = real(acos( sum(P1(jj,:).*P2(ii,:)) )*180/pi);
    if(cnPhi < smPhi) 
      smPhi   = [smPhi,cnPhi]; 
      smPos   = [smPos;[jj,ii]];                    % pos(ii: AIRS index. jj: CRIS)
      smTd    = [smTd,(iCnTim(ii) - cCnTim(jj))];
    end
  end
  if(~mod(jj,365)) fprintf('.'); end
end
%}
%{ 
jsz = size(smPhi,2);
figure(1);clf(1);plot((1:jsz),smPhi,'k.',(1:jsz-1),smTd*24*60,'g.');
  legend('smPhi','smTd mins');
figure(1);scatter(smTd*24*60, smPhi(1:end-1));xlabel('delay mins');ylabel('distance deg');
  grid on;
% save('/asl/s1/chepplew/projects/sno/airs_cris/HR/20130827_snapshot.mat');
%}

%        Compute separations and save indexes when criteria are met
fprintf(1,'Computing separations\n');
dist = 0.0; m = 0; k = 1;
pos  = [;]; tdiff = [];
tic
for i = 1:1:numel(cCn.Lat)
  for j = 1:1:numel(iCn.Tim)
    Dtim = abs(iCn.Tim(j) - cCn.Tim(i));
    if Dtim <= maxDtim      
      m = m+1;
      Phi = real(acos( sum(P1(i,:).*P2(j,:)) )*180/pi);    % P1: CrIS, P2: IASI.
      if Phi <= maxDphi                  % phi = 0.18 deg (0.0031 rad) => 20 km.
        dist(k) = Phi;                   % dist;
        pos     = [pos;[i,j]];           % pos(i: CrIS index. j: IASI)
        tdiff   = [tdiff,(iCn.Tim(j) - cCn.Tim(i))];
        k = k+1;
      end
    end
  end
  if(~mod(i,1000)) fprintf(1,'.'); end
end
toc
dist = real(dist);                         % can get complex phi.
fprintf(1,'Found %i matchups\n',k);

% save('/asl/s1/chepplew/projects/sno/iasi_cris/HR/20130827_t007d15km.mat');

% This gives us duplicate hits - so need to select only unique pairs.
%  this is done by sequentially testing for uniquness on each sensor.
if(size(pos,1) > 2 )
  un1 = []; un2 = []; upos = [];
  [x,ib,ix] = unique(pos(:,1)); un1  = [x,pos(ib,2)];
  [x,ib,ix] = unique(un1(:,2)); upos = [un1(ib,1),x]; clear un1;
  fprintf('Number unique SNOs %d %d\n',size(upos));

% Record time and space separations (pre reloading of files)
  sno=struct;
  sno.itim     = iCn.Tim(upos(:,2));    sno.ctim    = cCn.Tim(upos(:,1));
  sno.ilat     = iCn.Lat(upos(:,2));    sno.clat    = cCn.Lat(upos(:,1));
  sno.ilon     = iCn.Lon(upos(:,2));    sno.clon    = cCn.Lon(upos(:,1));
  sno.ifov     = iCn.Fov(upos(:,2));    sno.cfov    = cCn.Fov(upos(:,1));
  sno.isolz    = iCn.Solz(upos(:,2));   sno.csolz   = cCn.Sozn(upos(:,1));
  sno.isatz    = iCn.Satz(upos(:,2));   sno.csatz   = cCn.Sazn(upos(:,1));
  sno.iqual    = iCn.Qual(upos(:,2));
  sno.tdiff    = sno.itim - sno.ctim;    
%
  uP1    = P1(upos(:,1),:);      uP2 = P2(upos(:,2),:);    % P1: CRIS. P2: IASI
  for i = 1:size(upos,1)
     uPhi(i)   = real( acos(sum(uP1(i,:).*uP2(i,:))) * 180/pi );
  end
  sno.dist = uPhi;
  sno.upos = upos;

% ------------------------------------------------------------------------
% *********** SECTION 3 - Get Sensor Obs for the SNOs ***********
% ************  Load up CRIS SNO Obs data ************

  fprintf(1,'Loading CIRS SNO Obs\n');
  
  [fcx rcx] = load_cris_rad_for_ic_sno(cFnams, sno, cCn, vers);

% ****************  Re Load IASI SNO Obs data ************
  fprintf('Loading IASI SNO Obs\n');

  [rix] = load_iasi_rad_for_ic_sno(iaLst, sno, iCn); 

% ----------------------------------------------------------------
%  Convert IASI to CrIS 
% ----------------------------------------------------------------
  load('/asl/data/iremis/danz/iasi_f.mat');                      % fiasi [8461x1]

  if(exist('opt1','var')) 
    disp('Converting IASI to CrIS');
    [ns nz] = size(rix);
    if(ns ~= 8461 & nz == 8461) xri = rix'; else xri = rix; end
    [ri2c fi2c] = iasi2cris(xri,fiasi,opt1);

  end
  whos *i2c
% ----------------------------------------------------------------
%  Convert CrIS hires to CrIS mid or low res if needed
% ----------------------------------------------------------------

  if( all(strcmp(cris_res, {'HIGH','MID'})) | all(strcmp(cris_res, {'HIGH','LOW'})) )
    disp(['Converting CrIS SDR from FSR to MSR or NSR'])
  
%frq.LW = fcx(1:717);
%frq.MW = fcx(718:1586);
%frq.SW = fcx(1587:2223);
%rad.LW = rcx(1:717,:);
%rad.MW = rcx(718:1586,:);
%rad.SW = rcx(1587:2223,:);
  frq.LW  = fcx.LW;
  frq.MW  = fcx.MW;
  frq.SW  = fcx.SW;
  rad.LW  = rcx.LW;
  rad.MW  = rcx.MW;
  rad.SW  = rcx.SW;

  [fc rc] = translate_cris(frq, rad, 'midres', 1);


%{
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
%}

  elseif( strcmp(cris_res{1}, cris_res{2}) )
    disp(['Retaining original CrIS SDR spectral resolution'])
    rc2c = rcx;
    fc2c = fcx;
    clear rx fx;

  end


% ----------------------------------------------------------------
%             Save Data
% ----------------------------------------------------------------
  ri      = rix'; 
  savVars = {'strYr','strMn','strDy','maxDtim','maxDphi','sno',...
             'ri','rc','fc','ri2c','fi2c','fiasi',...
             'opt1','iaLst','crLst','iasi','cris'};
  savFN   = strcat(savDir, 'sno_iasi_cris_asl_',strYr,strMn,strDy,'_',vers,'.mat');
  disp(['Saving file: ' savFN]);
  save(savFN,savVars{:});


end   % ********** END if(size(pos,1) > 1)

%{
load('~strow/Matlab/Iasi/iasi_f.mat');     % fiasi [8461 x 1]
load('/home/chepplew/myLib/data/fcris_hires_2grd.mat');
np = size(rc,2);
bc = real(rad2bt(fc,rc));
bi = real(rad2bt(fiasi,ri'));
 ichn = find(fiasi > 900,1);     % 1022
 cchn = find(fc > 900,1);        % 404
 figure(3);clf(3);plot( [1:np],bc(cchn,:),'b.',[1:np],bi(ichn,:),'g.');
 figure(3);clf;scatter(dist,tdiff*24*60,'.');
 nanmean(bc(cchn,:),2) - nanmean(bi(ichn,:),2)  % 236.9370 - 236.9432
%}
