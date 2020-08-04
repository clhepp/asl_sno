function [] = make_AIRS_CRIS_sno_frmL1c(rdate, cris_res, cris, opt1, vers)
%
% function make_AIRS_CRIS_sno_frmL1c(rdate, cris_res, cris)
% 
% to produce files of SNO using CRIS L1C from CCAST processing at UMBC.ASL 
%   and AIRS L1C files.
% 
% INPUTS: rdate:    date to processes: string: 'YYYY/MM/DD'
%         cris_res: two strings: [{'high','mid','low'}{'high','mid','low'}]
%                   1st: CrIS source SDR spectral resolution.
%                   2nd: CrIS SNO spectral resolution (limits apply)
%         cris:     numeric: [1 or 2] CrIS mission NPP (1) or JPSS-1 (2)
%         opt1:     structure of fields to control airs2cris conversion,
%                   if entered with no fields no conversion takes place.
%         vers:     short string for a version reference to add to the SNO mat file.
%                   maximum allowed length 6 characters. 
% ------         -----------          ----------          ----------
% version 2: retain FOR -by- x-track -by- a-track array structures
%            include neighbours.
% 9.Jun.2020 CLH; new CrIS translation function. New load functions.
% 1 Aug 2020 CLH: new source directories 
%
% 1. Possible Sources:
%    /asl/cris/ccast/sdr{45,60}_{npp,j01}_{LR,MR,HR}/YYYY/JJJ/SDR_dYYYYMMdd*.mat
%       or /CrIS_SDR_npp_s45_dYYYYMMdd*.mat
%    /asl/airs/L1C_{v612,v672}/YYYY/JJJ/AIRS.2013.08.28.239.L1c.AIRS_Rad.v5.0.22.0.G13343142628.hdf
%    from Dec 2018 use: /asl/cris/ccast/sdr45_{npp,j01}_HR/2018/JJJ/SDR_dYYYYMMDD_t*.mat
%    
% 2. Destination:
%    /home/chepplew/data/sno/airs_cris{'',2}/ASL/{HR,MR,,LR}/YYYY/   
%   
% 
% Notes:  HM: L1a_err, in the SDR files is a 30 x nscan array, one value for each FOR:
%         1 = bad, 0 = OK.   
%

cd /home/chepplew/gitLib/asl_sno/sno_maker;
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
if(cris == 1); CX = '';  end
if(cris == 2); CX = '2'; end

% check version reference
vers = lower(vers);
if(length(vers) > 8) error('Version reference too long, 6 chars max'); return; end
allvers = {'v20a','v20d','a2v4_ref','a2v4_m20','a2v4_p20','a2_test1','','noaa','nasa'};
if(~ismember(vers,allvers)) error('Not valid CCAST data version');return; end
 
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
    savDir = ['/home/chepplew/data/sno/airs_cris' CX '/ASL/HR/'];
    if(strcmp(vers,'noaa'))
       savDir = ['/home/chepplew/data/sno/airs_cris' CX '/NOAA_pon/HR/'];
    end
    npLW = 717;  npMW = 869;   npSW = 637;
  case 'MID' 
    savDir = ['/home/chepplew/data/sno/airs_cris' CX '/ASL/MR/'];
    npLW = 717;  npMW = 865;   npSW = 633;
  case 'LOW'
    savDir = ['/home/chepplew/data/sno/airs_cris' CX '/ASL/LR/'];
    npLW = 717;  npMW = 437;   npSW = 163;
end

% Check conversion options, names and values
%opts  = struct;
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

% ************    Get source of CRIS SDR granules  *****************
if(isempty(CX) & ~strcmp(vers,'nasa') )
  %crDir  = ['/asl/data/cris/ccast/sdr60_npp_' CR '/'];               % pre-2018
  %crDir  = ['/asl/data/cris/ccast/sdr45_npp_' CR '/'];               % from 2018
  crDir  = ['/asl/cris/ccast/sdr45_npp_' CR '/'];                     % from Dec 2018
  crDir1 = sprintf('%s%s/%03d/',crDir, cYr,jday);
  %crLst  = dir(strcat(crDir1, 'SDR_d', cYr, cMn, cDy, '_t*.mat'));              % pre 2018
  crLst1 = dir(strcat(crDir1, 'CrIS_SDR_npp_s45_d', cYr, cMn, cDy,'_t*', vers, '.mat'));    % from 2018
  crDir2 = sprintf('%s%s/%03d/',crDir, cYr,jday2);
  crLst2 = dir([crDir2 'CrIS_SDR_npp_s45_d', cYr, sprintf('%02d', nMn2), sprintf('%02d',nDy2),'_t*.mat']);
  disp([strcat(crDir1, 'CrIS_SDR_npp_s45_d', cYr, cMn, cDy,'_t*.mat')]);
end
if(CX == '2'  & ~strcmp(vers,'nasa'))
  crDir  = ['/asl/cris/ccast/sdr45_j01_' CR '/'];
  crDir1 = sprintf('%s%s/%03d/',crDir, cYr,jday);
  disp(crDir1)
  crLst1 = dir(strcat(crDir1, 'CrIS_SDR_j01_s45_d', cYr, cMn, cDy,'_t*.mat'));
  crDir2 = sprintf('%s%s/%03d/',crDir, cYr,jday2);
  crLst2 = dir([crDir2 'CrIS_SDR_j01_s45_d', cYr, sprintf('%02d', nMn2), sprintf('%02d',nDy2),'_t*.mat']);
  if(strcmp(vers,'noaa'))
    crDir  =  '/home/chepplew/data/cris/noaa/sdr_j01_hr/';
    crDir1 = sprintf('%s%s/%03d/',crDir, cYr,jday);
    crLst  = dir(strcat(crDir1, 'CrIS_SDR_j01_s45_d', cYr, cMn, cDy,'_tXXX_noaa_PLon_*.mat'));
  end
end
if(strcmp(vers,'nasa') & CX == 1)
  crhome = '/asl/isilon/cris/nasa_l1b/npp/';
  crDir1 = sprintf('%s%s/%03d/',crhome, cYr,jday);
  crDir2 = sprintf('%s%s/%03d/',crhome, cYr,jday+1);
  crLst1 = dir(strcat(crDir1, 'SNDR.SNPP.CRIS.', cYr, cMn, sprintf('%02d',nDy1),'T*.nc')); 
  crLst2 = dir(strcat(crDir2, 'SNDR.SNPP.CRIS.', cYr, cMn, sprintf('%02d',nDy2),'T*.nc')); 
end  

if(numel(crLst1) < 1) error('Insufficient CrIS CCAST data'); return; end

disp(['Found ' num2str(numel(crLst1)) ' CrIS L1c Files']);
disp([crLst1(1).folder '  ' crLst1(1).name])
disp(['cris_res: ' cris_res ' cris mission: ' num2str(cris) ' vers: ' vers])

% Collect all granule filenames:
for i=1:numel(crLst1) 
  cFnams{i} = strcat(crLst1(i).folder, '/', crLst1(i).name);
end
% Get first granule of the following day and append to full list:
if(numel(crLst2))
  cFnams{numel(crLst1)+1} = strcat(crLst2(1).folder,'/',crLst2(1).name);
end


% ------------ Load CrIS geo data --------------

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

% --------------  Load AIRS geo data  ----------------------
fprintf('Loading AIRS geo\n'); 
arDir = strcat('/umbc/hpcnfs1/strow/asl/airs/L1C/',cYr,'/',sprintf('%03d',jday),'/');
arLst = dir(sprintf('%sAIRS.%4d.%02d.%02d*.hdf',arDir,nYr1,nMn1,nDy1));
disp(['Found ' num2str(numel(arLst)) ' AIRICRAD Files']);

if(numel(arLst) < 1) error('Insufficient AIRS L1C data'); return; end

alat  = [];  alon = []; atim = []; aatrak = []; axtrak = []; asolzn = []; agran = [];
alnfr = []; ascnt = [];
for fn = 1:numel(arLst);
  l1c = read_airs_l1c(strcat(arLst(fn).folder,'/',arLst(fn).name));
  alat     = [alat;   l1c.Latitude];
  alon     = [alon;   l1c.Longitude]; 
  atim     = [atim;   airs2dnum(l1c.Time)];               % convert to matlab time.
  aatrak   = [aatrak; l1c.atrack];
  axtrak   = [axtrak; l1c.xtrack];
  asolzn   = [asolzn; l1c.solzen];
  agran    = [agran;  l1c.gindex];
  alnfr    = [alnfr;  l1c.landFrac];
  ascnt    = [ascnt;  repmat(l1c.scan_node_type',1,90)];   % scan node type is once per scan line
  fprintf('.');
end

fprintf('\n');
whos alat alon atim aatrak axtrak asolzn agran alnfr ascnt;

aCntr = [43 44 45 46 47 48];
aCn.Tim  = atim(:,aCntr)';            % NB ': to match CrIS center arrays
aCn.Lat  = alat(:,aCntr)';
aCn.Lon  = alon(:,aCntr)';
aCn.atr  = aatrak(:,aCntr)';
aCn.xtr  = axtrak(:,aCntr)';
aCn.Grn  = agran(:,aCntr)';
aCn.Sozn = asolzn(:,aCntr)';
aCn.Lnfr = alnfr(:,aCntr)';
aCn.Scnt = ascnt(:,aCntr)';
disp(['Found ' num2str(numel(aCn.Tim)) ' AIRS center track FOVs'])

% Get ascending flag (1=asc, 0=desc)
lat45   = alat(:,45);  clear junk;
for i=1:length(lat45)-1 junk(i+1) = lat45(i+1)-lat45(i); end
  junk(length(lat45)) = junk(i);
aAsc = junk;
aAsc(aAsc<=0) = 0;
aAsc(aAsc>0)  = 1;
clear junk;
% Expand to match near-nadir arrays:
aAsc = repmat(aAsc,6,1);


%  ****************** SECTION 3: Get indexes of SNOs  ********************
%            pre-compute position vectors - saves a heap of time *********
%            NaNs are conserved
fprintf('computing position vectors\n');
P1 = zeros(numel(aCn.Lat),3);                              % AIRS
P2 = zeros(numel(cCn.Lat),3);                              % CrIS

for ii = 1:numel(aCn.Lat)                                   % Number of Swath.
    P1(ii,:) = [ cos(aCn.Lat(ii)*pi/180.0) * cos(aCn.Lon(ii)*pi/180.0), ...
         cos(aCn.Lat(ii)*pi/180.0)*sin(aCn.Lon(ii)*pi/180.0), sin(aCn.Lat(ii)*pi/180.0) ];
end
for ii = 1:numel(cCn.Lat)
    P2(ii,:) = [ cos(cCn.Lat(ii)*pi/180.0) * cos(cCn.Lon(ii)*pi/180.0), ...
         cos(cCn.Lat(ii)*pi/180.0)*sin(cCn.Lon(ii)*pi/180.0), sin(cCn.Lat(ii)*pi/180.0) ];
end    
whos P1 P2;

%        Compute separations and save indexes when criteria are met
fprintf('Computing separations\n');
dist = 0.0; m = 0; k = 1;
pos = []; tdiff = [];
tic
for jj = 1:1:numel(cCn.Lat)
  for ii = 1:1:numel(aCn.Lat) % arst:1:arsp
    Dtim = abs(aCn.Tim(ii) - cCn.Tim(jj));
    if Dtim <= maxDtim       
      m = m+1;
      % dist = distance(arCnLat(ja),arCnLon(ja),iaCnLat(ji),iaCnLon(ji)); % way too slow
      Phi = real(acos( sum(P1(ii,:).*P2(jj,:)) )*180/pi);
      if Phi <= maxDphi                % phi = 0.18 deg (0.0031 rad) => 20 km.
        dist(k) = Phi;                % dist;
        pos     = [pos;[ii,jj]];      % pos(ii: AIRS index. jj: CRIS)
        tdiff   = [tdiff,(aCn.Tim(ii) - cCn.Tim(jj))];
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
if(size(pos,1) > 1 )
  un1 = [];  un2 = [];  upos = [];
  [x,ib,ix] = unique(pos(:,1)); un1  = [x,pos(ib,2)];
  [x,ib,ix] = unique(un1(:,2)); upos = [un1(ib,1),x]; clear un1;
  fprintf('Number unique SNOs %d %d\n',size(upos));
  nSNO = size(upos,1);
% Record time and space separations (pre reloading of files)
  sno = struct;
  sno.aTime     = aCn.Tim(upos(:,1));    sno.cTime   = cCn.Tim(upos(:,2));
  sno.aLat      = aCn.Lat(upos(:,1));    sno.cLat    = cCn.Lat(upos(:,2));
  sno.aLon      = aCn.Lon(upos(:,1));    sno.cLon    = cCn.Lon(upos(:,2));
  sno.aAtrack   = aCn.atr(upos(:,1));    sno.cAtrack = cCn.Atr(upos(:,2));
  sno.aXtrack   = aCn.xtr(upos(:,1));    sno.cXtrack = cCn.Xtr(upos(:,2));
  sno.aGran     = aCn.Grn(upos(:,1));    sno.cGran   = cCn.Gid(upos(:,2));
  sno.aAsc      = aAsc(upos(:,1));       sno.cFov    = cCn.Fov(upos(:,2));
  sno.aSolzen   = aCn.Sozn(upos(:,1));   sno.cSolzen = cCn.Sozn(upos(:,2));
  sno.alnfrac   = aCn.Lnfr(upos(:,1));   sno.clnfrac = [];
  sno.aScnt     = aCn.Scnt(upos(:,1));
  sno.tdiff     = sno.aTime - sno.cTime;    
%
  uP1    = P1(upos(:,1),:);      uP2 = P2(upos(:,2),:);    % P1: AIRS, P2: CRIS.
  for i = 1:size(upos,1)
     uPhi(i)   = real( acos(sum(uP1(i,:).*uP2(i,:))) * 180/pi );
  end
  sno.dist = uPhi';

% ---------------------------------------------------------------------
% **** Section 4 Reload AIRS and CRIS files and save only required Obs.
% ---------------------------------------------------------------------

  [fax, rax] = load_airs_rad_for_sno(arLst, upos, aCn);
 
  [fcx, rcx] = load_cris_rad_for_sno(cFnams, upos, npLW, cCn, vers);


% ----------------------------------------------------------------
%  Convert CrIS hires to CrIS mid or low res if needed 
% ----------------------------------------------------------------
  if( all(strcmp(cris_res, {'HIGH','MID'})) || all(strcmp(cris_res, {'HIGH','LOW'})) )
    disp(['Converting CrIS SDR from FSR to MSR or NSR']);
    resx = opt1.user_res;
    rcy.LW = rcx.LW;
    rcy.MW = rcx.MW;
    rcy.SW = rcx.SW;
    
    [fc rc] = translate_cris(fcx, rcy, resx, 1);
  elseif( strcmp(cris_res{1}, cris_res{2}) )
    disp(['Retaining original CrIS SDR spectral resolution'])
    rc = [rcx.LW; rcx.MW; rcx.SW];
    fc = [fcx.LW; fcx.MW; fcx.SW];
    clear rcx fcx;
  end


% ----------------------------------------------------------------
%  Convert AIRS to CrIS 
% ----------------------------------------------------------------
  sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
  % new L1C ref to Sep2010.
  sfile = '/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf';

  load('/home/chepplew/projects/cris/cris_freq_2grd.mat'); fcris = vchan;   % fcris [1317x1]
  % get AIRS frequency grid 
  load('/home/chepplew/projects/airs/airs_f.mat');  % f [2378] fairs [2645] AIRS frequency grids
      f2645 = fairs;                     % the filled 2645 AIRS frequency grid
  kn = 1000;

  opt1
  if(exist('opt1','var')) 
    switch opt1.user_res
      case 'lowres', opt1.cfile = 'corr_lowres.mat';
      case 'midres', opt1.cfile = 'corr_midres.mat';
      case 'hires',  opt1.cfile = 'corr_hires.mat';
    end  
    disp('Converting AIRS to CrIS');
    [ns sz] = size(rax.rad);
    if(ns ~= 2645 & sz == 2645) xra = rax.rad'; sz=ns; else xra = rax.rad; end
    ra2c = [];
    for j = 1:kn:sz
      sx = j:min(j+kn-1,sz);
      [rrx fa2c] = airs2cris(xra(:,sx),f2645,sfile,opt1);
      ra2c       = [ra2c, rrx];
      fprintf(1,'.');
    end

    disp('Completed conversion airs2cris');
  end
  opt1

%%%%%%%%%%%%%%%%%%%% Saving data %%%%%%%%%%%%%%%%%%%%%
  clear atime alat alon asolzen asatzen ra;
  clear ctime clat clon cfov csolzen csatzen tdiff dist;
  nbr_rLW        = rcx.nbr_rLW;
  nbr_ra         = [];
  ra             = xra; 
  l1cProc        = rax.Prc;
  l1cSynthReason = rax.Reas;
  airs.f      = cell2mat(l1c.freq);  
  airs.idchns = cell2mat(l1c.chanID);
  savFN   = strcat(savDir, cYr, ...
            '/sno_airs_cris_asl_wngbr_',cYr,cMn,cDy,'_frmL1c_',vers,'.mat');
  savVars = {'sno','ra','airs','rc','fc','ra2c','fa2c','upos','maxDtim','maxDphi',...
             'nbr_rLW','nbr_ra','rdate','cris_res','arDir','arLst','cFnams','l1cProc',...
	     'l1cSynthReason','opt1','sfile'};
  fprintf(1,'Saving data to file: %s\n',savFN);
  save(savFN, savVars{:});


end    % end size(pos) > 1

