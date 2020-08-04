function [] = make_AIRS_CRIS_sno_frmL1b(req_date, cris_res)
%
% function [] = make_AIRS_CRIS_sno_frmL1b(req_date, cris_res)
% 
% to produce files of SNO using CRIS hi-res from CCAST processing at UMBC.ASL
%            and AIRIBRAD granules
% 
% 1. Sources & Inputs:
%    /asl/data/cris/ccast/sdr60_hr/YYYY/JJJ/SDR_d20130828_t0006589.mat
%    /asl/data/airs/AIRIBRAD/YYYY/JJJ/AIRS.2013.08.28.239.L1B.AIRS_Rad.v5.0.22.0.G13343142628.hdf
%    
% 2. Destination & Outputs:
%    /asl/s1/chepplew/projects/sno/airs_cris/HR/   
% 
% Notes:  HM: L1a_err, in the SDR files is a 30 x nscan array, one value for each FOR:
%         1 = bad, 0 = OK.   
%
% version 2: retain FOR -by- x-track -by- a-track array structures

cd /home/chepplew/projects/sno/makeSNO;
addpath /asl/matlib/time                             % tai2dnum
addpath /asl/rtp_prod/airs/readers                   % readl1b_all
addpath /asl/rtp_prod/airs/utils                     % f_default_1lb.mat
warning('off');
run cris_neighbor_LUT                                % load CrIS adjacent FOV look-up tables

% Check date string and record the following day
whos req_date; disp(req_date); fprintf('\n');
try 
   D = datenum(req_date,'yyyy/mm/dd');
   Dp1 = D+1;
catch
   error('Incorrect Date Format')
   return
end

% Check CrIS Resolution
cris_res = upper(cris_res);
all_res  = {'HIGH','LOW'};
if(~ismember(cris_res,all_res))
  error('Invalid CrIS resolution. Choices are: high or low');
  return
end
  
if(strcmp(cris_res,'HIGH')) 
  savDir = '/asl/s1/chepplew/data/sno/airs_cris/HR/';
  strRes = 'HR';
  npLW = 717;  npMW = 869;   npSW = 637;
end
if(strcmp(cris_res,'LOW'))
  savDir = '/asl/s1/chepplew/data/sno/airs_cris/ASL/LR/';
  strRes = 'LR';
  npLW = 717;  npMW = 437;   npSW = 163;
end


% Criteria for separation and delay
maxDtim   = 0.006944;             % day. 600/86400 secs =  Mattime is decimal day (20 mins=0.0139)
maxDphi   = 0.071900;             % deg. 0.07 = 7.8 km Earth radius: 6371 km. 1 deg = 111 km.

% Get day to process & convert to required formats
[nYr1 nMn1 nDy1] = datevec(D);
[nYr2 nMn2 nDy2] = datevec(Dp1);
cYr     = req_date(1:4);     cMn = req_date(6:7);     cDy = req_date(9:10);
day0    = sprintf('%4d/%02d/%02d',nYr1 - 1,12,31);
jday    = datenum(req_date) - datenum(day0);  clear junk;           % needed for CRIS directory
jday2   = fix(Dp1 - datenum(day0));

% ************    Get CRIS SDR granule files  *****************
if(strcmp(cris_res,'HIGH')) crDir = '/asl/data/cris/ccast/sdr60_hr/';  end
if(strcmp(cris_res,'LOW'))  crDir = '/asl/data/cris/ccast/sdr60/';  end
crDir1 = sprintf('%s%s/%03d/',crDir, cYr,jday);
crLst  = dir(strcat(crDir1, 'SDR_d', cYr, cMn, cDy,'_t*.mat'));
disp(crLst(1).folder)
disp(['Found ' num2str(numel(crLst)) ' CrIS CCAST L1c Files']);
% Get first granule of the following day and append to full list:
crDir2 = sprintf('%s%s/%03d/',crDir, cYr,jday2);
crLst2 = dir([crDir2 'SDR_d', cYr, sprintf('%02d', nMn2), sprintf('%02d',nDy2),'_t*.mat']);
for i=1:numel(crLst) cFnams{i} = strcat(crLst(i).folder, '/', crLst(i).name); end
cFnams{numel(crLst)+1} = strcat(crLst2(1).folder,'/',crLst2(1).name);
 
clat = [];   clon = []; ctim = [];  cfov = [];  cxtrak = []; catrak = []; cgran = [];
csozn = []; csazn = [];
for fn = 1:numel(crLst);
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
  
% Get granule ID, should return one value per aTrack (60)
  clear granID junk
  for i=1:length(geo.Granule_ID)
    if(strcmp(geo.Granule_ID(i,1:3),'NPP'))
      junk(i) = str2num(geo.Granule_ID(i,4:end));
    end
  end
% Expand to full size 
  B = repmat(junk,9,1,sz2); 
  C = permute(B,[1 3 2]);
  C(:,L1a_err) = NaN; 
  cgran = cat(3,cgran, C);

  errFlg = reshape(L1a_err,[],1);                     % (30 x nscan) 1: bad, 0: good.
%   junk  = reshape(geo.FORTime(~L1a_err),[],1);       % <- (30 x 60) per granule. nu-secs since 1958
%  clat   = [clat; reshape(geo.Latitude(:,~L1a_err),[],1)]; 
%  clon   = [clon; reshape(geo.Longitude(:,~L1a_err),[],1)];         % <- (16200 x 1) per granule
%  cfov   = [cfov;   reshape(tmpv(:,~L1a_err),[],1)];   clear tmpv;
%  cxtrak = [cxtrak; reshape(tmpx(:,~L1a_err),[],1)];   clear tmpx;
%  catrak = [catrak; reshape(tmpa(:,~L1a_err),[],1)];   clear tmpa;

   junk  = geo.FORTime;   junk(L1a_err) = NaN;
  tmpt   = iet2dnum(junk);    clear junk;
  % Need to expand ctim to match every FOV so that sub-setting works
  junk   = repmat(tmpt,[1 1 9]);
  tmpt   = permute(junk,[3 1 2]);    clear junk;
%  ctim   = [ctim; reshape(junk,[],1)]; clear junk;           % <- (16200 x 1) times for every FOV
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

% *******************   Get AIRS daily files  *******************
fprintf('Loading AIRS geo\n'); 
arDir = strcat('/asl/data/airs/AIRIBRAD/',cYr,'/',sprintf('%03d',jday),'/');
arLst = dir(sprintf('%sAIRS.%4d.%02d.%02d*.hdf',arDir,nYr1,nMn1,nDy1));
disp(['Found ' num2str(numel(arLst)) ' AIRIBRAD Files']);

alat  = []; alon = []; atim = []; aatrak = []; axtrak = []; asolzn = []; agran = [];
alnfr = [];
for fn = 1:numel(arLst);
 [eqXtai,f,gdata] = readl1b_all(strcat(arDir,arLst(fn).name));
 alat     = [alat; single(gdata.rlat)'];
 alon     = [alon; single(gdata.rlon)'];
 atim     = [atim; airs2dnum(gdata.rtime)'];               % convert to matlab time.
 aatrak   = [aatrak; single(gdata.atrack)'];
 axtrak   = [axtrak; single(gdata.xtrack)'];
 asolzn   = [asolzn; single(gdata.solzen)'];
 agran    = [agran; single(gdata.findex)'];
 alnfr    = [alnfr; single(gdata.landfrac)'];
 fprintf('.');
end

fprintf('\n');
whos alat alon atim aatrak axtrak asolzn agran alnfr;

aCntr = [43 44 45 46 47 48];
cntr = find(axtrak == 43 | axtrak == 44 | axtrak == 45 | axtrak == 46 | ...
            axtrak == 47 | axtrak == 48);
disp(['Found ' num2str(numel(cntr)) ' AIRS center track FOVs'])

aCnTim  = atim(cntr);
aCnLat  = alat(cntr);
aCnLon  = alon(cntr);
aCnatr  = aatrak(cntr);
aCnxtr  = axtrak(cntr);
aCnGrn  = agran(cntr);
aCnSozn = asolzn(cntr);
aCnLnfr = alnfr(cntr);

% Get ascending flag (1=asc, 0=desc)
cntr45  = find(axtrak == 45);
lat45   = alat(cntr45);  clear junk;
for i=1:length(lat45)-1 junk(i+1) = lat45(i+1)-lat45(i); end
  junk(length(lat45)) = junk(i);
aAsc = junk;
aAsc(aAsc<=0) = 0;
aAsc(aAsc>0)  = 1;
clear junk;
% Expand to match near-nadir arrays:
aAsc = repmat(aAsc,6,1);
aAsc = aAsc(:);


%  ****************** SECTION 3: Get indexes of SNOs  ********************
%            pre-compute position vectors - saves a heap of time *********
%            NaNs are conserved
fprintf('computing position vectors\n');
P1 = [;];                % AIRS position
P2 = [;];                % CRIS position
for ii = 1:numel(aCnLat)                                   % Number of Swath.
    P1 = [P1;[ cos(aCnLat(ii)*pi/180.0) * cos(aCnLon(ii)*pi/180.0), ...
         cos(aCnLat(ii)*pi/180.0)*sin(aCnLon(ii)*pi/180.0), sin(aCnLat(ii)*pi/180.0) ]];
end
for ii = 1:numel(cCnLat)
    P2 = [P2;[ cos(cCnLat(ii)*pi/180.0) * cos(cCnLon(ii)*pi/180.0), ...
         cos(cCnLat(ii)*pi/180.0)*sin(cCnLon(ii)*pi/180.0), sin(cCnLat(ii)*pi/180.0) ]];
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
pos = [;]; tdiff = [];
tic
for jj = 1:1:numel(cCnLat)
  arst = max(1,find(aCnTim > cCnTim(jj),1) - wndta);
  arsp = min(arst + 2*wndta, length(aCnLat));
  for ii = arst:1:arsp
    Dtim = abs(aCnTim(ii) - cCnTim(jj));
    if Dtim <= maxDtim       
      m = m+1;
      % dist = distance(arCnLat(ja),arCnLon(ja),iaCnLat(ji),iaCnLon(ji)); % way too slow
      Phi = real(acos( sum(P1(ii,:).*P2(jj,:)) )*180/pi);
      if Phi <= maxDphi               % phi = 0.18 deg (0.0031 rad) => 20 km.
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
  sno.aGran     = aCnGrn(upos(:,1));    sno.cGran   = cCnGid(upos(:,2));
  sno.aAsc      = aAsc(upos(:,1));      sno.cFov    = cCnFov(upos(:,2));
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
% **** Section 4 Reload AIRS and CRIS files and save only required Obs.
% ---------------------------------------------------------------------
  fprintf('loading AIRS obs\n');
  ps = 1; pe = 0; aCount = 0;  k = 1;
  tarad   = [;];  taAtrk = []; taXtrk = [];  tasozn = [];
  tasazn  = [];   talnfr = [];  tatim = [];  talat  = []; talon  = []; 
  nbr_ra  = zeros(nSNO,2,2378);
  for fn = 1:numel(arLst); % fn = 1;
    [eqXtai,f,gdata] = readl1b_all(strcat(arDir,arLst(fn).name));
    atim     = airs2dnum(gdata.rtime);                   % convert to matlab time.
    tCnInd   = find(gdata.xtrack == 43 | gdata.xtrack == 44 | gdata.xtrack == 45 | ...
                    gdata.xtrack == 46 | gdata.xtrack == 47 | gdata.xtrack == 48);
    pe       = ps + numel(gdata.rlat(tCnInd)) - 1;
    ext_indx = find(ismember(gdata.xtrack, [42:49]));
    %  tests which data points we want to keep match to appropriate FOR
    fprintf('%d: %d %d ,',fn,ps,pe);
    junk   = ismember(upos(:,1),[ps:pe]);       aSams = upos(junk,1);    clear junk;
    if ( length(aSams >= 1)  )
      fprintf('%d \t',length(aSams));
      % subset onto center track b4 using pos index.
      taXtrk = [taXtrk, gdata.xtrack(tCnInd(aSams-ps+1))];
      tarad  = [tarad, gdata.robs1(:,tCnInd(aSams-ps+1))];
      
%{      
      indx = find(gdata.xtrack == 43 | gdata.xtrack == 44 | gdata.xtrack == 45 | ...
                  gdata.xtrack == 46 | gdata.xtrack == 47 | gdata.xtrack == 48);
      junk = gdata.robs1(:,indx);  tarad  = [tarad, junk(:,aSams-ps+1) ]; clear junk;
      junk = gdata.xtrack(indx);   taXtrk = [taXtrk, junk(aSams-ps+1)];   clear junk;
      junk = gdata.atrack(indx);   taAtrk = [taAtrk, junk(aSams-ps+1)];   clear junk;
      junk = gdata.rlat(indx);     talat  = [talat, junk(aSams-ps+1)];    clear junk;
      junk = gdata.rlon(indx);     talon  = [talon, junk(aSams-ps+1)];    clear junk;
      junk = gdata.solzen(indx);   tasozn = [tasozn, junk(aSams-ps+1)];   clear junk;
      junk = atim(indx);           tatim  = [tatim, junk(aSams-ps+1)];    clear junk;
      junk = gdata.landfrac(indx); talnfr = [talnfr,junk(aSams-ps+1)];    clear junk;
      junk = gdata.satzen(indx);   tasazn = [tasazn,junk(aSams-ps+1)];    clear junk;
      %aCnra    = gdata.robs1(:,indx);
%}
      junk     = gdata.xtrack(tCnInd);
      for j = 1:length(aSams)
        xxt = junk(aSams(j)-ps+1);
        if(aSams(j)-ps > ps)    nbr_ra(k,1,:) = gdata.robs1(:,tCnInd(aSams(j)-ps)); end
	if(aSams(j)-ps <= ps)   nbr_ra(k,1,:) = gdata.robs1(:,tCnInd(aSams(j)-ps+1) - 1); end
	if(aSams(j)+ps+2 <= pe) nbr_ra(k,2,:) = gdata.robs1(:,tCnInd(aSams(j)-ps+2)); end
	if(aSams(j)+ps+2 > pe)  nbr_ra(k,2,:) = gdata.robs1(:,tCnInd(aSams(j)-ps+1) + 1); end
        k = k+1;
      end
    end
    ps = pe+1;
    aCount = aCount + length(aSams);  fprintf('%d \n',aCount);
%   fprintf('.');
  end
  fprintf('\n');
  % whos talat talon tatim talnfr taAtrk taXtrk tarad;

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

    if ( length(cSams >= 1) )
      fprintf(1,'\t %d ',length(cSams));
      nrLW = [];

      tcLat   = [tcLat;  cCnLat(cSams-ps+1)];
      tcLon   = [tcLon;  cCnLon(cSams-ps+1)];
      tcTim   = [tcTim;  cCnTim(cSams-ps+1)];
      tcxtrk  = [tcxtrk; cCnXtr(cSams-ps+1)];
      tcatrk  = [tcatrk; cCnAtr(cSams-ps+1)];
      % Get landfrac
      tcSozn  = [tcSozn; cCnSozn(cSams-ps+1)];
      tcFov   = [tcFov;  cCnFov(cSams-ps+1)];
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
	%elseif(cCnAtr(cSams(j)) == 1 & map(4,1) == -1 )
        else
	  trLW(4,:) = rLWp1(:,map(4,3),map(4,2),cCnAtr(cSams(j)) );
	end
      % record time stamps for each neighbour for validation
        % TBD
        %nrLW   = cat(1, nrLW, trLW);
	nbr_rLW(k,:,:) = trLW;
	k = k+1;
      end         % end for cSams
    end           % end: if cSams
    ps = pe + 1;
    cCount = cCount + length(cSams); fprintf('%d \n',cCount);
  end             % for fn loop
  % whos tcLat tcLon tcxtrk tcatrk tcTim tcFov tcrLW tcrMW tcrSW nbr_rLW

%%%%%%%%%%%%%%%%%%%% Saving data %%%%%%%%%%%%%%%%%%%%%
  clear atime alat alon asolzen asatzen ra;
  clear ctime clat clon cfov csolzen csatzen rc tdiff dist;
  ra    = tarad;  fa = f; 
  rc    = [tcrLW; tcrMW; tcrSW]; fc = [vLW; vMW; vSW];
  savFN   = strcat(savDir, cYr, '/sno_airs_cris_asl_w4ngbr_',cYr,cMn,cDy,'_asc.mat')
  savVars = {'sno','ra','fa','rc','fc','upos','maxDtim','maxDphi','nbr_rLW',...
             'nbr_ra','req_date','cris_res','arDir','arLst','cFnams'};
  fprintf('Saving data to file: %s\n',savFN);
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
