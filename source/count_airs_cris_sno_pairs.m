function count_airs_cris_sno_pairs()


cyear = '2013';

snoD1  = ['/home/chepplew/data/sno/airs_cris/ASL/LR/' cyear '/'];
snoFn1 = ['sno_airs_cris_asl_wngbr_' cyear '*_frmL1c.mat'];
d1     = dir(strcat(snoD1,snoFn1));
disp(['Found ' num2str(length(d1)) ' type 1 SNO files'])

snoD2  = ['/home/chepplew/data/sno/airs_cris/ASL/LR/' cyear '/'];
snoFn2 = ['sno_airs_cris_asl_wngbr_' cyear '*_frmL1c.mat'];
d2     = dir(strcat(snoD2,snoFn2));
disp(['Found ' num2str(length(d2)) ' type 1 SNO files'])


disp(['Found ' num2str(length(d1)) ' set 1 and ' num2str(length(d2)) ' set 2 files']);

num1 = [];  dist1 = []; tdiff1 = [];
for ifn1=1:length(d1)
  g = load(strcat(snoD1, d1(ifn1).name),'sno');
  num1   = [num1; length(g.sno.aLat)];
  dist1  = [dist1; g.sno.dist];
  tdiff1 = [tdiff1; g.sno.tdiff];
  fprintf(1,'.');
end

num2 = []; dist2 = []; tdiff2 = [];
for ifn2=1:length(d2)
  g = load(strcat(snoD2, d2(ifn2).name),'sno');
  num2   = [num2; length(g.sno.aLat)];
  dist2  = [dist2; g.sno.dist];
  tdiff2 = [tdiff2; g.sno.tdiff];
  fprintf(1,'.');
end

fprintf(1,'%8d SNO pairs from type-1 files\n', sum(num1))
fprintf(1,'%8d SNO pairs from type-2 files\n', sum(num2))
