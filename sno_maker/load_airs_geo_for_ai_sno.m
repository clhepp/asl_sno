function ar = load_airs_geo_for_ai_sno(airslist,airsvers)

%
%
%
%
%

addpath /home/chepplew/gitLib/asl_sno/sno_maker
addpath /home/chepplew/gitLib/asl_sno/utilities
addpath /home/chepplew/myLib/matlib                  % read_airs_l1c
addpath /asl/matlib/time

ar.lat   = []; ar.lon = [];  ar.tim = []; ar.atrak = []; ar.xtrak = [];
ar.solzn = []; ar.lanfr = []; ar.xobs = []; ar.proc = []; ar.reas = [];

switch airsvers
  case 'L1B'
  for fn = 1:numel(airslist);
    [eqXtai,f,gdata] = readl1b_all([airslist(fn).folder '/' airslist(fn).name]);
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

  case 'L1C'
  for fn = 1:numel(airslist);
    l1c = read_airs_l1c([airslist(fn).folder '/' airslist(fn).name]);
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

end

% Subset center track FOVs
switch airsvers
  case 'L1B'
    aCnInd = find(ar.xtrak == 43 | ar.xtrak == 44 | ar.xtrak == 45 | ar.xtrak == 46 | ...
              ar.xtrak == 47 | ar.xtrak == 48);
    ar.cntim  = ar.tim(aCnInd);
    ar.cnlat  = ar.lat(aCnInd);
    ar.cnlon  = ar.lon(aCnInd);
    ar.cnsozn = ar.solzn(aCnInd);

  case 'L1C'
    aCntr    = [43 44 45 46 47 48];
    ar.cntim   = reshape(ar.tim(:,aCntr),[],1);
    ar.cnlat   = reshape(ar.lat(:,aCntr),[],1);
    ar.cnlon   = reshape(ar.lon(:,aCntr),[],1);
    ar.cnatr   = ar.atrak(:,aCntr);
    ar.cnxtr   = ar.xtrak(:,aCntr);
    ar.cnsozn  = ar.solzn(:,aCntr);
    ar.cnlnfr  = ar.lanfr(:,aCntr);
    %ar.cnobs   = reshape(ar.xobs(:,aCntr),[],1);
    %ar.cnproc  = reshape(permute(ar.proc(:,aCntr,:),[3,1,2]),2645,[],1);
    %ar.cnreas  = reshape(permute(ar.reas(:,aCntr,:),[3,1,2]),2645,[],1);
   disp(['Found ' num2str(numel(ar.cntim)) ' AIRS center track FOVs'])

end


