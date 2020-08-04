function [] = make_airs_iasi_sno(rdate,airs,iasi,opts);
%
%   function: make_airs_iasi_sno(A,B)
%   find SNO coincidences between AIRS and IASI on given date, save data.
%   
%   inputs: rdate: the date string, format YYYY/MM/DD
%           airs: the AIRS data level string. Options: {'l1b','l1c'}
%           iasi: Mission no. Options: [1,2] for IASI-1 or IASI-2.
%           opts: structure of fields to control iasi2cris
%
% close all; clear all;

cd      /home/chepplew/projects/sno/makeSNO/;
addpath /home/chepplew/myLib/matlib                  % read_airs_l1c
addpath /asl/matlib/time                             % tai2dnum
addpath /home/chepplew/gitLib/asl_sno/sno_maker
addpath /home/chepplew/gitLib/asl_sno/utilities    
addpath /asl/matlib/aslutil                          % findfiles

addpath /asl/packages/iasi_decon
addpath /asl/matlib/h4tools                          % h4sdread

warning('off');

load('/home/chepplew/projects/airs/airs_f.mat');
load('/home/chepplew/projects/iasi/f_iasi.mat');
 ach = find(fairs > 900,1);
 ich = find(f_iasi > 900,1);

% Check date string and record the following day
whos rdate; disp(rdate); fprintf('\n');
try 
   D = datenum(rdate,'yyyy/mm/dd');
   Dp1 = D+1;
catch
   error('Incorrect Date Format')
   return
end

% Get day to process & convert to required formats
[nYr1 nMn1 nDy1] = datevec(D);
[nYr2 nMn2 nDy2] = datevec(Dp1);
cYr     = rdate(1:4);     cMn = rdate(6:7);     cDy = rdate(9:10);
day0    = sprintf('%4d/%02d/%02d',nYr1 - 1,12,31);
jday    = datenum(rdate) - datenum(day0);  clear junk;           % needed for CRIS directory
jday2   = fix(Dp1 - datenum(day0));

fprintf('Processing %s\n',rdate);
datStr = strcat(rdate(1:4),rdate(6:7),rdate(9:10));
jday = sprintf('%03d',datenum(rdate) - datenum(strcat(cYr, '/01/01')) + 1);

% Check AIRS data source
airs = upper(airs)
if(~ismember(airs,{'L1B','L1C'})) error('Invalid AIRS data type'); return; end

% CHeck IASI mission number
if(~ismember(iasi,[1,2,3])); error('Invalid IASI mission number'); end
if(iasi == 1) IX = 'M02'; IR = '';  end
if(iasi == 2) IX = 'M01'; IR = '2'; end
if(iasi == 3) IX = 'M03'; IR = '3'; end

% Check conversion options, names and values
%opts  = struct;
opts.sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
opts
names = fieldnames(opts);
if(isempty(names)) clear opts; 
elseif(~all(isfield(opts,{'sfile'})))
  disp('Not all options for iasi2airs are present, plz chceck');
end

% Criteria for separation and delay
maxDtim   = 0.0139; % 0.00694;                % day. 600/86400 secs =  Mattime is decimal day (20 mins=0.0139)
maxDphi   = 0.18;   % 0.07;                   % deg. 0.07 = 7.8 km Earth radius: 6371 km. 1 deg = 111 km.

% ----------------------------------------------
% Get IASI data
% ----------------------------------------------
fprintf('loading IASI data\n');
switch iasi
  case 1 
    iadir = strcat('/asl/data/IASI/L1C/',cYr,'/',cMn,'/',cDy,'/');
  case 2 
    iadir = strcat('/asl/data/IASI2/L1C/',cYr,'/',cMn,'/',cDy,'/');
  case 3
    iadir = strcat('/asl/data/IASI3/L1C/',cYr,'/',jday,'/');
end

iaLst = dir(strcat(iadir,'IASI_xxx_1C_',IX,'*'));
disp(['Found ' num2str(numel(iaLst)) ' IASI granules in dir: ' iadir ' (IASI ' IX ')']);

if(numel(iaLst) < 1) error('No IASI granules found'); return; end

ia = load_iasi_geo_for_sno(iaLst);

%{
iafov  = [];  ialat = []; ialon = []; iatim = []; iasolzen = []; iaatrak = [];
iaAtrk = []; iaXtrk = []; iqual = [];
x_iObs = [];
for fn = 1: numel(iaLst);
 s = readl1c_epsflip_all(strcat(iadir,iaLst(fn).name)); 
 iafov    = [iafov;    reshape(s.IASI_FOV,[],1)];              % [N x 4] arrays
 ialat    = [ialat;    reshape(s.Latitude,[],1)];
 ialon    = [ialon;    reshape(s.Longitude,[],1)];
    junk  = s.Time2000;                                        % 2000 epoch
 iatim    = [iatim;    reshape(iasi2mattime(junk),[],1)];      % convert to matlab time
 iasolzen = [iasolzen; reshape(s.Solar_Zenith,[],1)];
 iaatrak  = [iaatrak;  reshape(s.Scan_Line,[],1)];             % [690 x 4] values 1:23 (4 times)/granule
 iqual    = [iqual;    reshape(s.GQisFlagQual,[],1)];
    %junk  = permute(s.IASI_Radiances,[2,1,3]);
 x_iObs   = [x_iObs; reshape(s.IASI_Radiances(:,:,ich),[],1)];
 fprintf(1,'.');
end
fprintf(1,'\n');
  
% create across-track array. Check along track has 30 FORs values of 1:23
 iaxtrak = zeros(size(iaatrak,1)/4,4);
 for j = 1:fix(size(iaatrak,1)/120)              % retain row x col structure
  is = (j-1)*30 +1;  ie = j*30;
  iaxtrak(is:ie,[1 2 3 4]) = [1:30; 1:30; 1:30; 1:30]';    
 end
iaxtrak  = reshape(iaxtrak,[],1);

%  whos iafov ialat ialon iatim iasolzen iaatrak iaxtrak iqual x_iObs

% get center-track geo-location for finding SNOs
icntr    = find(ia.xtrak == 15 | ia.xtrak == 16);
iCnLat   = ia.lat(icntr);
iCnLon   = ia.lon(icntr);
iCnTim   = ia.tim(icntr);
iCnFov   = ia.fov(icntr);
iCnSzn   = ia.solzen(icntr);
iCnQual  = ia.qual(icntr);
%%iCnObs = ia.xobs(icntr);
disp(['Kept ' num2str(length(iCnTim)) ' IASI center track FORs']);
  whos iCnLat iCnLon iCnTim x_iCnObs
%}
%{
% screen for geo errors (can get invalid latitudes)
 igeo_err = find(iCnLat < -90);

 iCnLat(igeo_err) = NaN;
 iCnLon(igeo_err) = NaN;
 iCnTim(igeo_err) = NaN;
 arr = ones(size(iCnLat));  arr(igeo_err) = NaN;
%} 
% -------------------------------------------------
% Get AIRS data
% -------------------------------------------------
% NB. AIRS v6.1.0.0 only between 2015-Apr to end 2015-Jul. v6.1.2.0 otherwise.
fprintf('loading AIRS data\n');
%ardir = strcat('/asl/data/airs/AIRIBRAD/',strDate,'/');
%arLst = dir(strcat(ardir,'AIRS.',yr,'.',mn,'.',dy,'.*.L1B.AIRS_Rad*.hdf'));
ardir = strcat('/asl/data/airs/',airs,'/',cYr,'/',jday,'/');
arLst = dir(strcat(ardir,'AIRS.',cYr,'.',cMn,'.',cDy,'.*.',airs,'.AIRS_Rad.v*.hdf'));
disp(['Found ' num2str(numel(arLst)) ' AIRS granules'])

if(numel(arLst) < 1) error('No AIRS granules found'); return; end

ar = load_airs_geo_for_ai_sno(arLst,airs);

%{
ar.lat   = []; ar.lon = [];  ar.tim = []; ar.atrak = []; ar.xtrak = []; 
ar.solzn = []; ar.lanfr = []; ar.xobs = []; ar.proc = []; ar.reas = [];

if L1B
for fn = 1:numel(arLst);
  [eqXtai,f,gdata] = readl1b_all(strcat(ardir,arLst(fn).name));
  ar.lat     = [ar.lat,   gdata.rlat];
  ar.lon     = [ar.lon,   gdata.rlon];
    junk     = gdata.rtime;                           % AIRS epoch
  ar.tim     = [ar.tim,   tai2mattime(junk)];             % convert to matlab time.
  ar.atrak   = [ar.atrak,  gdata.atrack];
  ar.xtrak   = [ar.xtrak,  gdata.xtrack];
  ar.solzn   = [ar.solzn, gdata.solzen];
  fprintf('.');
end
fprintf('\n');
end
if L1C
for fn = 1:numel(arLst);
  l1c = read_airs_l1c(strcat(arLst(fn).folder,'/',arLst(fn).name));
  ar.lat     = [ar.lat;   l1c.Latitude];
  ar.lon     = [ar.lon;   l1c.Longitude];
    junk     = l1c.Time;                             % AIRS epoch
  ar.tim     = [ar.tim;   airs2dnum(junk)];             % convert to matlab time.
  ar.atrak   = [ar.atrak;  l1c.atrack];
  ar.xtrak   = [ar.xtrak;  l1c.xtrack];
  ar.solzn   = [ar.solzn; l1c.solzen];
  ar.lanfr   = [ar.lanfr; l1c.landFrac];
   %ar.xobs  = [ar.xobs; l1c.radiances(:,:,ach)];
  %ar.proc   = [ar.proc;  l1c.L1cProc];
  %ar.reas   = [ar.reas;  l1c.L1cSynthReason];
 fprintf('.');
end
fprintf('\n');
end
ar

% Subset center track FOVs
if L1B
aCnInd = find(ar.xtrak == 43 | ar.xtrak == 44 | ar.xtrak == 45 | ar.xtrak == 46 | ...
              ar.xtrak == 47 | ar.xtrak == 48);
aCnTim = artim(aCnInd);
aCnLat = arlat(aCnInd);
aCnLon = arlon(aCnInd);
aCnSzn = arsolzn(aCnInd);
end

if L1C
aCntr    = [43 44 45 46 47 48];
aCnTim   = reshape(ar.tim(:,aCntr),[],1); 
aCnLat   = reshape(ar.lat(:,aCntr),[],1);
aCnLon   = reshape(ar.lon(:,aCntr),[],1);
aCnatr   = ar.atrak(:,aCntr);
aCnxtr   = ar.xtrak(:,aCntr);
aCnSozn  = ar.solzn(:,aCntr);
aCnLnfr  = ar.lanfr(:,aCntr);
%aCnObs  = reshape(ar.xobs(:,aCntr),[],1);
%aCnProc  = reshape(permute(ar.proc(:,aCntr,:),[3,1,2]),2645,[],1);
%aCnReas  = reshape(permute(ar.reas(:,aCntr,:),[3,1,2]),2645,[],1); 
disp(['Found ' num2str(numel(aCnTim)) ' AIRS center track FOVs'])
end
%}


% ---------- pre-compute position vectors - saves a heap of time
fprintf('computing position vectors\n');
P1 = zeros(numel(ar.cnlat),3);                              % AIRS
P2 = zeros(numel(ia.cnlat),3);                              % IASI
for i = 1:numel(ar.cnlat)              % Number of centres.
    P1(i,:) = [cos(ar.cnlat(i)*pi/180.0) * cos(ar.cnlon(i)*pi/180.0), ...
         cos(ar.cnlat(i)*pi/180.0)*sin(ar.cnlon(i)*pi/180.0), sin(ar.cnlat(i)*pi/180.0) ];
end
for j = 1:numel(ia.cnlat)
    P2(j,:) = [cos(ia.cnlat(j)*pi/180.0) * cos(ia.cnlon(j)*pi/180.0), ...
         cos(ia.cnlat(j)*pi/180.0)*sin(ia.cnlon(j)*pi/180.0), sin(ia.cnlat(j)*pi/180.0) ];
end    
  whos P1 P2



%{
% find the time and sample numbers of closest approach
minPhi = 1000.0;
for ji = 1:10:length(iaCnLat)
  for ja = 1:100:length(arCnLat)
     Phi = acos( sum(P1(:,:,ja+1).*P2(:,:,ji+1)) )*180/pi;
     if (Phi < minPhi) minPhi = Phi; minPos = [ji,ja]; end
  end
  [mintDff, arTIndx] = min( abs(arCnTim - iaCnTim(ji)) );
  if(rem(ji,1000) == 0 ) fprintf('.'); end
end
[x,iaTIndx] = min( abs(arCnTim(arTIndx) - iaCnTim) );
%}

fprintf(1,'Computing matchups\n')
dist = 0.0; m = 0; k = 1;
pos = [];   tdiff = [];
tic
for jj = 1:1:numel(ia.cnlat)
   for ii = 1:1:numel(ar.cnlat)
      Dtim = abs(ia.cntim(jj) - ar.cntim(ii));
      if Dtim <= maxDtim
         m = m+1;
         Phi = acos( sum(P1(ii,:).*P2(jj,:)) )*180/pi;
	 if Phi <= maxDphi                % phi = 0.18 deg (0.0031 rad) => 20 km.
	    dist(k) = Phi;                % dist;
	    pos     = [pos;[ii,jj]];      % pos(ii: AIRS index. jj: IASI index)
	    tdiff   = [tdiff,(ar.cntim(ii)-ia.cntim(jj))];
	    k = k+1;
	 end
      end
   end
   if(rem(jj,1001) == 0 ) fprintf('.'); end
end
toc
fprintf('\n');
fprintf('There are %d matchups\n',size(pos,1));

%{
% Search for close separations given time range of obs. (10 mins = 0.006944 hrs) 
fprintf('searching for matches\n');
AIpos = [;];
for ja = 1:length(arCnLat)
  ParTim = arCnTim(ja);
  junk = find(iaCnTim <= ParTim + 0.0069 & iaCnTim > ParTim - 0.0069);
%  find their spatial separation
  for k = 1:numel(junk)
    phi(k) = acos( sum(P1(:,:,ja+1).*P2(:,:,junk(k)+1)) )*180/pi;
    if(abs(phi(k)) <= 0.18) AIpos = [AIpos,[ja;junk(k)]]; end
  end
%figure(2); plot(phi); pause(2);clf(2);
  if(rem(ja,500) == 0) fprintf('.'); end
end
fprintf('\n');
%}

% This gives us duplicate hits - so need to select only unique pairs.
sno = struct;
if(size(pos,1) > 1 )
  upos = [];
  un1 = []; un2 = [];
  [x,ib,ix] = unique(pos(:,1)); un1  = [x,pos(ib,2)];
  [x,ib,ix] = unique(un1(:,2)); upos = [un1(ib,1),x]; clear un1;
%  [x,ib,ix] = unique(pos(:,1));
%  upos = [x,pos(ib,2)];               % keep (:,1)=AIRS. (:,2)=IASI.
  fprintf('There are %d unique SNOs\n',size(upos,1));

% Record time and space separations
  sno.atim    = aCnTim(upos(:,1));   sno.itim = iCnTim(upos(:,2));
  sno.tdiff   = sno.atim - sno.itim;
  sno.alat    = aCnLat(upos(:,1));    sno.ilat    = iCnLat(upos(:,2));
  sno.alon    = aCnLon(upos(:,1));    sno.ilon    = iCnLon(upos(:,2));
  sno.atim    = aCnTim(upos(:,1));    sno.itim    = iCnTim(upos(:,2));
  sno.asolzen = aCnSozn(upos(:,1));   sno.isolzen = iCnSzn(upos(:,2));
  sno.alandfr = aCnLnfr(upos(:,1));     
                                      sno.iqual   = iCnQual(upos(:,2));
                                      sno.ifov    = iCnFov(upos(:,2));
%
  uP1    = P1(upos(:,1),:); uP2 = P2(upos(:,2),:);
  for i=1:size(upos,1) 
    sno.dist(i)  = real( acos(sum(uP1(i,:).*uP2(i,:))) * 180/pi ); 
  end
%
%  l1cProc        = aCnProc(:,upos(:,1));
%  l1cSynthReason = aCnReas(:,upos(:,1));
%{
% Save SNO geolocation here:
  savVars = {'arLst','iaLst','maxDphi','maxDtim','sno','pos','upos',...
             'P1','P2','x_aCnObs','x_iCnObs','rdate','airs','iasi'};
  savFn = '/home/chepplew/data/sno/airs_iasi/ASL/20180123_sno_geo.mat';
  save(savFn, savVars{:});
%}    
% -----------------------------------------------------
% Reload AIRS granules and save only required Obs.
% ------------------------------------------------------
  fprintf('loading AIRS obs\n');
  ps = 1; pe = 0; aCount = 0; 
  aObs  = zeros(size(upos,1),2645); 
  aProc = uint8(zeros(size(upos,1),2645)); 
  aReas = uint8(zeros(size(upos,1),2645));
  if L1B
    for fn = 1:numel(arLst); % fn = 1;
      [eqXtai,f,gdata] = readl1b_all(strcat(ardir,arLst(fn).name));
      tXtrk = gdata.xtrack;                      % select centre track FOVs
      tCnInd = find(tXtrk == 43 | tXtrk == 44 | tXtrk == 45 | tXtrk == 46 | ...
               tXtrk == 47 | tXtrk == 48);
      pe = ps + numel(gdata.rlat(tCnInd)) - 1;

    %  tests which data points we want to keep match to appropriate FOR
      aInter = ismember(upos(:,1),[ps:pe]);                      % returns 1 at positions where values match
      if (length(find(aInter)) >= 1) 
         jnkObs  = gdata.robs1(:,tCnInd);
         arSnObs = [arSnObs, jnkObs(:,upos(aInter,1)-ps+1)]; clear jnkObs;
      end
    % Update counters:
      ps = pe+1;
      cumInter = cumInter + length(find(aInter));
      fprintf('.');
    end
  fprintf('\n');
  end
  %%%%
  if L1C
    for fn = 1:numel(arLst);
      l1c  = read_airs_l1c(strcat(arLst(fn).folder,'/',arLst(fn).name));
      for j = 1:numel(sno.alat)
        xv = ismember(l1c.Latitude,  sno.alat(j));
        yv = ismember(l1c.Longitude, sno.alon(j));
        zv = ismember(airs2dnum(l1c.Time), sno.atim(j));
        awant = intersect(find(xv),intersect(find(yv),find(zv)));
        if(awant)
	  aCount = aCount + 1; 
	  disp(['fn ' num2str(fn) '   aCount ' num2str(aCount)]) 
          %disp(['fn  numel(awant): ' num2str(fn) ': ' num2str(numel(awant))]);
          [I,J] = ind2sub(size(l1c.Latitude),awant(1));
          aObs(aCount,:)  = squeeze(l1c.radiances(I,J,:))';
	  aProc(aCount,:) = squeeze(l1c.L1cProc(I,J,:))';
	  aReas(aCount,:) = squeeze(l1c.L1cSynthReason(I,J,:))';
        end
      end
      fprintf('.')
    end
  end

%{
      pe   = ps + numel(l1c.Latitude(:,aCntr)) - 1;
      junk = ismember(upos(:,1),[ps:pe]);       aSams = upos(junk,1);  clear junk;
      fprintf('%d\t %d %d\t %d ',fn,ps,pe,length(aSams));
      if (length(aSams) >= 1) 
        cnObs   = reshape(l1c.radiances(:,aCntr,:),[],2645);
        aObs    = [aObs; cnObs(aSams-ps+1,:)]; clear cnObs;
      end
      % Update counters:
      ps = pe+1;
      aCount = aCount + length(aSams);  
      fprintf(1,'%d\n', aCount);
    end
%}

% ----------------------------------------------
%  Load IASI SNO Obs
% ----------------------------------------------
  fprintf('\n loading IASI obs\n');
  ps = 1; pe = 0; iSams = []; iObs = zeros(size(upos,1),8461); iCount = 0;

  for fn = 1:numel(iaLst);
    s = readl1c_epsflip_all(strcat(iadir,iaLst(fn).name));
    for j = 1:numel(sno.ilat)
     xv = ismember(s.Latitude,  sno.ilat(j));
     yv = ismember(s.Longitude, sno.ilon(j));
     zv = ismember(iasi2mattime(s.Time2000),  sno.itim(j));
     iwant = intersect(find(xv),intersect(find(yv),find(zv)));
       if(iwant)
         iCount = iCount + 1;
	 disp(['fn ' num2str(fn) '   iCount ' num2str(iCount) ' j ' num2str(j)]) 
         %disp(['fn  numel(iwant): ' num2str(fn) ': ' num2str(numel(iwant))]);
         [I,J] = ind2sub(size(s.Latitude),iwant(1));
         iObs(iCount,:)  = squeeze(s.IASI_Radiances(I,J,:))';
       end
     end
     fprintf('.')
   end

%{
    [nx ny ~]  = size(s.IASI_Radiances);
    [mx my]    = size(s.Latitude);
    if(nx ~= mx | ny ~= my) error('IASI radiance and latitude array mis-match'); end
    t_ilat     = reshape(s.Latitude,[],1);
    t_ilon     = reshape(s.Longitude,[],1);
    t_itim     = iasi2mattime(reshape(s.Time2000,[],1));
    t_iatrak   = reshape(s.Scan_Line,[],1);
    t_ixtrak   = zeros(size(t_iatrak,1)/4,4);
    for j = 1:fix(size(t_iatrak,1)/120)              % retain row x col structure
     is = (j-1)*30 +1;  ie = j*30;
     t_ixtrak(is:ie,[1 2 3 4]) = [1:30; 1:30; 1:30; 1:30]';    
    end
    t_ixtrak  = reshape(t_ixtrak,[],1); 
% subset on center track (FORs 15 and 16) to update counters
    tCnInd   = find(t_ixtrak == 15 | t_ixtrak == 16);
    mx = size(s.Scan_Line,1)/30;
    xCnInd = []; for i=1:mx xCnInd(i,:)   = [15*i,16*i]; end; xCnInd = sort(xCnInd(:));
    tCnLat   = t_ialat(tCnInd);
    num_iCn  = size(s.Latitude,1)/15 * 4;
    pe = ps + length(tCnInd) - 1;
    junk = ismember(upos(:,2),[ps:pe]);   iSams = upos(junk,2);  clear junk;
    fprintf('%d\t %d %d %d\t',fn,ps,pe,length(iSams));
    if (length(iSams) >= 1)
      t_ilat   = t_ilat(tCnInd);
      t_ilon   = t_ilon(tCnInd);
      t_itim   = t_itim(tCnInd);
      x_ilat   = [x_ilat; t_ilat(iSams-ps+1)];
      x_ilon   = [x_ilon; t_ilon(iSams-ps+1)];
      x_itim   = [x_itim; t_itim(iSams-ps+1)];
      [nx ny ns] = size(s.IASI_Radiances);
      %iCnObs = s.IASI_Radiances(tCnInd,:,:);
      junk   = permute(s.IASI_Radiances(xCnInd,:,:),[1 2 3]);
      junk2  = reshape(junk(:,:,ich),[],1);
      %t_irad = junk2(tCnInd,:);
      iObs   = [iObs; junk2(iSams-ps+1,:)]; clear junk junk2 t_irad;
    end
  % Update counters:
    ps = pe+1;
    iCount = iCount + length(iSams);
    fprintf(1,'%d\n', iCount);
  end
%}  
  whos iObs aObs;
  
%{
 figure(1);clf;plot(sno.alat, sno.alon,'+', sno.ilat,sno.ilon,'o')
 figure(2);clf;scatter(sno.dist,sno.tdiff*24*60,[])
 %[nx nz] = size(aCnObs)
 figure(3);clf;plot([1:77],aCnObs(pos(:,1)),'.-', [1:77],iCnObs(pos(:,2)),'o-')


figure(4);clf;plot([1:17], aObs(:,ach),'.-', [1:17], iObs(:,ich),'o-')
figure(3);clf;plot([1:77], x_aCnObs(pos(:,1)),'.-', [1:77],x_iCnObs(pos(:,2)),'o-')
figure(3);clf;plot([1:17], x_aCnObs(upos(:,1)),'.-', [1:17],x_iCnObs(upos(:,2)),'o-')

%}

% ------------------------------------------------
% Optional convert iasi2airs
% ------------------------------------------------
% get IASI user-grid parameters
iasi_p = iasi_params;
ifrq1  = iasi_p.freq;

kn = 1000;

if(exist('opts','var')) 
  disp('Converting IASI to AIRS');
  [ns sz] = size(iObs);
  if(ns ~= 8461 & sz == 8461) xri = iObs'; sz=ns; else xri = iObs; end
  ri2a = [];
  for j = 1:kn:sz
    sx = j:min(j+kn-1,sz);
    [rrx fi2a] = iasi2airs(xri(:,sx),ifrq1,opts.sfile,fairs);
    ri2a       = [ri2a, rrx];
    fprintf(1,'.');
  end

disp('Completed conversion iasi2airs');
end

% ---------------------------------------------
% Save SNO data
% ---------------------------------------------
  ra = aObs';
  ri = iObs';
  l1cProc        = aProc';
  l1cSynthReason = aReas';
  savDir = ['/home/chepplew/data/sno/airs_iasi' IR '/ASL/' cYr '/'];
  savFil = strcat('sno_airs_iasi_',datStr,'_frm',airs,'.mat');
  fprintf('saving %s\n',savFil);
  savVars = {'rdate','arLst','iaLst','upos','sno','ra','ri',...
             'airs','iasi','opts','ri2a','fi2a','l1cProc','l1cSynthReason'};
  save(strcat(savDir,savFil),savVars{:}); 
	      
end    %% <<<<<< if upos >=1  >>>>>>>>	   

%{
% Get spectrum for SNOs.
% IASI requires both FOV and time (four FOVs per time stamp)
 adtim = datetime(sno.atim,'ConvertFrom','datenum');
 idtim = datetime(sno.itim,'ConvertFrom','datenum');







%}
