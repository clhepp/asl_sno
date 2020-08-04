function hsk = get_orbit_phase_sno(sdate)
%
% Calculate orbit phase for AIRS frequency Drift Correction
%    from existing AIRS:CRIS SNO files.
% collect basic housekeeping: latitude, longitude and time.
%
% INPUT: [string] Date for processing, format 'yyyy/MM/dd'
%
% OUTPUT: structure hsk with fields:
%        asc   [N x 1] single.   AIRS ascending flag (1: asc, 0: desc)
%        phi   [N x 1] single.   AIRS orbit phase (0:180) per LLS definition.
%        alat  [N x 1] single.   AIRS Observation Latitude
%        alon  [N x 1] single.   AIRS Observation Longitude
%        atime [N x 1] datetime. AIRS Observation time.
%
% C Hepplewhite. 21-Nov-2017
%
% 
% Check date string and record the following day
whos sdate; disp(sdate); fprintf('\n');
try 
   D = datetime(sdate,'format','yyyy/MM/dd');
   Dp1 = D+1;
catch
   error('Incorrect Date Format')
   return
end

cyear = num2str(year(D));
cmon  = sprintf('%02i', month(D));
cday  = sprintf('%02i', day(D));

snoD  = ['/home/chepplew/data/sno/airs_cris/ASL/LR/' cyear '/'];
snoFn = ['sno_airs_cris_asl_wngbr_' cyear cmon cday '_frmL1c.mat'];
d     = dir(strcat(snoD,snoFn));
 
%disp(['found ' num2str(length(d)) ' SNO files'])


ifn=1;
g = load([snoD d(ifn).name],'sno');

% Get orbit phase (0: desc eqX, 45: S.P., 90: asc eqX , 135: N.P., 180: desc eqX.)
% with interpolation
nobs = size(g.sno.aLat,1);
disp(['Number of SNO pairs: ', num2str(nobs)]);
phi  = [];
for i = 1:nobs
  if(g.sno.aAsc(i) == 1) 
    %phi(i)   = fix(0.5*(g.sno.aLat(i) + 180)); end
    phi(i,:) = [floor(0.5*(g.sno.aLat(i) + 180)), 0.5*(g.sno.aLat(i) + 180), ...
                   ceil(0.5*(g.sno.aLat(i) + 180))]; end
  if(g.sno.aAsc(i) == 0) 
    %phi(i)   = fix(0.5*mod(360 - g.sno.aLat(i), 360)); end
    phi(i,:) = [floor(0.5*mod(360 - g.sno.aLat(i), 360)), ...
          0.5*mod(360 - g.sno.aLat(i), 360), ceil(0.5*(mod(360 - g.sno.aLat(i), 360)))]; 
  end
end

% reduce the sample by distancd and delay (to match original set: Dt=0.006944 Ds=0.0719)
ii = find(g.sno.dist  <= 0.071900 & g.sno.tdiff <= 0.006944);
ii=':';

% Get Lat, Lon, time
hsk.asc  = g.sno.aAsc(ii);
hsk.phi  = phi(ii,:);
hsk.alat = g.sno.aLat(ii);
hsk.alon = g.sno.aLon(ii);
hsk.atime = datetime(g.sno.aTime(ii),'convertFrom','datenum');


