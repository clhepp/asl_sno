function hsk = run_get_sno_housekeeping(nyear)
%
% function run_get_sno_housekeeping(nyear)
%
% a wrapper for get_orbit_phase_sno.m
%
% Used for AIRS frequency drift correction work to accompany stand-alone SNO
%    radiance files in /home/chepplew/data/sno/airs_cris/radiance/
%      
% Process a year of ASL AIRS:CRIS low-res SNO mat files.
% 
%
cd /home/chepplew/gitLib/asl_sno/run

addpath /home/chepplew/gitLib/asl_sno/source

srcd  = '/home/chepplew/data/sno/airs_cris/ASL/LR/2013/';
dlist = dir([srcd 'sno_airs_cris_asl_wngbr_*_frmL1c.mat']);
disp(['Found ' num2str(length(dlist)) ' SNO files'])

savd  = '/home/chepplew/data/sno/airs_cris/radiances/';

hsk = [];
for ifn = 1:2:length(dlist)
  junk = regexp(dlist(ifn).name,'[0-9]','match');                    % format yyyy/MM/dd
  sdate = [cell2mat(junk(1:4)) '/' cell2mat(junk(5:6)) '/' cell2mat(junk(7:8))];
  hsk   = [hsk; get_orbit_phase_sno(sdate)];
  fprintf(1,'.')
end
fprintf(1,'\n');
  
snum = 0; alat = []; alon = []; atime = []; phi = []; asc = [];
for i = 1:length(dlist)  
   snum  = snum + length(hsk(i).alat);
   alat  = [alat; hsk(i).alat];
   alon  = [alon; hsk(i).alon];
   atime = [atime; hsk(i).atime];
   phi   = [phi; hsk(i).phi];
   asc   = [asc; hsk(i).asc];
end
disp(['total no. SNO pairs included in housekeeping: ' num2str(snum)])

% apply QC screening/subsetting from load_sno_airs_cris_asl_mat()
alat_ig  = alat(s.ig);
alon_ig  = alon(s.ig);
atime_ig = atime(s.ig);
phi_ig   = phi(s.ig)';
asc_ig   = asc(s.ig);

whos *_ig

% save these data
svars = {'alat_ig','alon_ig','atime_ig','phi_ig','asc_ig'};
savfn = '/home/chepplew/data/sno/airs_cris/radiances/2013even_asl_ac_sno_hskeeping_mw.mat';
save(savfn, svars{:});
