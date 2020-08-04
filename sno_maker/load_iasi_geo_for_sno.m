function ia = load_iasi_geo_for_sno(filelist)

%
%
%
%
%
%
%

addpath /home/chepplew/gitLib/asl_sno/sno_maker
addpath /home/chepplew/gitLib/asl_sno/utilities


ia.fov  = [];  ia.lat = []; ia.lon = []; ia.tim = []; ia.solzen = []; 
ia.atrak = []; ia.atrk = []; ia.xtrk = []; ia.qual = [];
ia.xobs = []; ia.satzen = [];

for fn = 1: numel(filelist);
 s = readl1c_epsflip_all([filelist(fn).folder '/' filelist(fn).name]);
 ia.fov    = [ia.fov;    reshape(s.IASI_FOV,[],1)];              % [N x 4] arrays
 ia.lat    = [ia.lat;    reshape(s.Latitude,[],1)];
 ia.lon    = [ia.lon;    reshape(s.Longitude,[],1)];
    junk  = s.Time2000;                                        % 2000 epoch
 ia.tim    = [ia.tim;    reshape(iasi2mattime(junk),[],1)];      % convert to matlab time
 ia.solzen = [ia.solzen; reshape(s.Solar_Zenith,[],1)];
 ia.satzen = [ia.satzen; reshape(s.Satellite_Zenith,[],1)];
 ia.atrak  = [ia.atrak;  reshape(s.Scan_Line,[],1)];             % [690 x 4] values 1:23 (4 times)/granule
 ia.qual   = [ia.qual;   reshape(s.GQisFlagQual,[],1)];
    %junk  = permute(s.IASI_Radiances,[2,1,3]);
 %ia.xobs   = [ia.xobs; reshape(s.IASI_Radiances(:,:,ich),[],1)];
 fprintf(1,'.');
end
fprintf(1,'\n');

% create across-track array. Check along track has 30 FORs values of 1:23
 ia.xtrak = zeros(size(ia.atrak,1)/4,4);
 for j = 1:fix(size(ia.atrak,1)/120)              % retain row x col structure
  is = (j-1)*30 +1;  ie = j*30;
  ia.xtrak(is:ie,[1 2 3 4]) = [1:30; 1:30; 1:30; 1:30]';
 end
ia.xtrak  = reshape(ia.xtrak,[],1);

% get center-track geo-location for finding SNOs
icntr    = find(ia.xtrak == 15 | ia.xtrak == 16);
ia.cnlat   = ia.lat(icntr);
ia.cnlon   = ia.lon(icntr);
ia.cntim   = ia.tim(icntr);
ia.cnfov   = ia.fov(icntr);
ia.cnsolzn = ia.solzen(icntr);
ia.cnsatzn = ia.satzen(icntr);
ia.cnqual  = ia.qual(icntr);
%%iCnObs = ia.xobs(icntr);
disp(['Kept ' num2str(length(ia.cntim)) ' IASI center track FORs']);

%{
% screen for geo errors (can get invalid latitudes)
 igeo_err = find(iCnLat < -90);

 iCnLat(igeo_err) = NaN;
 iCnLon(igeo_err) = NaN;
 iCnTim(igeo_err) = NaN;
 arr = ones(size(iCnLat));  arr(igeo_err) = NaN;
%}

