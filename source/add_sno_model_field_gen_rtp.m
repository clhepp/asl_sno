% function add_sno_model_field_gen_rtp.m
%
% Take existing SNO data file, add model fields and create RTP suitable
% for forward model calc.
%
%

addpath /home/chepplew/gitLib/asl_sno/source       % load sno
addpath /asl/matlib/time                           % dnum2tai
addpath /home/chepplew/gitLib/rtp_prod2/grib       % fill_{era, ecmwf, merra}
addpath /home/chepplew/gitLib/rtp_prod2/util       % set_attr
addpath /home/chepplew/gitLib/rtp_prod2/emis       % rtp_add_emis

% hardwire configuration settings
cfg.model = 'ecmwf';

% load of a chunk of SNO
sno = load_sno_airs_cris_asl_mat('2018/01/01','2018/01/03',[1:1305],'low',1,'v20a');

% Use AIRS Observations - cross ref channels for AIRS SARTA
load('/home/chepplew/myLib/data/airs_f2378.mat');    % afchan [2378]
[sf si] = sort(x.afchan);
[xi xj] = seq_match(sf, sno.fa);

nchan = size(xj,1);
vchan = sno.fa(xj);
ichan = xi;
nobs  = size(sno.aLat,1);

% create RTP structure fields:
head = struct;
head.pfields = 5;  % robs1, no calcs. 5 = (1) prof + (4) obs
head.ptype   = 0;
head.ngas    = 2;
head.glist   = int32([1 3]');
head.gunit   = int32([21 21]');

% Assign header variables
head.instid = 800;           % AIRS
head.pltfid = -9999;
head.nchan  = nchan;
head.ichan  = ichan;
head.vchan  = vchan;
head.vcmax  = max(head.vchan);
head.vcmin  = min(head.vchan);


% Assign header attribute strings
hattr={ {'header' 'pltfid' 'Aqua'}, ...
        {'header' 'instid' 'AIRS'}, ...
        {'header' 'githash' ' '}, ...
        {'header' 'rundate' ' '}, ...
        {'header' 'klayers_exec' ' '}, ...
        {'header' 'sarta_exec' ' '} };

% write prof structure
prof = struct;
prof.robs1    = sno.ra(xj,:);
prof.rlon     = sno.aLon';
prof.rlat     = sno.aLat';
prof.rtime    = dnum2tai(sno.aTime');
prof.landfrac = sno.alnfr';              % no CrIS landfrac
prof.solazi   = zeros(1,nobs);
prof.solzen   = ones(1,nobs).*150.0;
prof.satazi   = zeros(1,nobs);           % nadir
prof.satzen   = zeros(1,nobs);           % nadir
prof.salti    = zeros(1,nobs);
prof.scanang  = zeros(1,nobs);
prof.zobs     = ones(1,nobs).*705000.0;
prof.upwell   = ones(1,nobs);

% write pattr structure
pattr = struct;
pattr = { ...
   {'profiles', 'iudef(1,:)', 'Dust flag:[1=true,0=false,-1=land,-2=cloud,-3=bad data] {dustflag}'}, ...
   {'profiles', 'iudef(2,:)', 'Dust_score'}, ...
   {'profiles', 'iudef(3,:)', 'L1cProc'}, ...
   {'profiles', 'iudef(4,:)', 'scan_node_type'}, ...
   {'profiles', 'iudef(5,:)', 'L1cSynthReason'}, ...
   {'profiles', 'iudef(6,:)', 'SceneInhomogeneous'}, ...
   {'profiles', 'udef(6,:)',  'NeN'} };


% Add in model data
fprintf(1, '>>> Add model: %s...', cfg.model)
switch cfg.model
  case 'ecmwf'
    [prof,head,pattr]  = fill_ecmwf(prof,head,pattr);
  case 'era'
    [prof,head,pattr]  = fill_era(prof,head,pattr);
  case 'merra'
    [prof,head,pattr]  = fill_merra(prof,head,pattr);
end

% adding surface emissivity and reflectivity
fprintf(1, '>>> Running rtp_add_emis...');
try
    [prof,pattr] = rtp_add_emis(prof,pattr);
catch
    fprintf(1, '>>> ERROR: rtp_add_emis failure');
    %return;
end

% Save RTP file (only short filenames allowed)
savdr   = '/home/chepplew/data/sno/rtp/';
savfn   = 'sno_airs_cris_lr_2018d001e003.ip.rtp';
disp(['Saving RTP data to: ' [savdr savfn] ]);
rtpwrite([savdr 'temp.rtp'], head, hattr, prof, pattr);
syscmd = ['mv ' [savdr 'temp.rtp'] ' ' [savdr savfn]];
[status,cmdout] = system(syscmd);

fprintf(1, 'Done\n');


% convert to layers file for sarta:
addpath /home/sbuczko1/git/matlib/clouds/sarta

klayers_exec = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';
sarta_exec   = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';
code0        = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1        = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';

run_sarta.klayers_code    = klayers_exec;
run_sarta.sartaclear_code = sarta_exec;
run_sarta.sartacloud_code = code0;
run_sarta.cloud           =+1;
run_sarta.clear           =+1;
run_sarta.cumsum          =9999;   % was -1

% Run KLAYERS
%  rtpwrite(outP, hxn,ha,pxn,pa);
rtpin  = [savdr savfn];
rtpout = [savdr 'temp.op.rtp'];
errF   ='/home/chepplew/logs/klayers/klayers.err';
  klayerser = ['!' klayers_exec ' fin=' rtpin ' fout=' rtpout ' >& ' errF];
  eval(klayerser);
[hd2,ha2,pd2,pa2] = rtpread(rtpout);

pd2.clwc = prof.clwc;
pd2.ciwc = prof.ciwc'
pd2.cc   = prof.cc;
pd2.tcc  = prof.tcc;

[prof0, oslabs] = driver_sarta_cloud_rtp(hd2,ha2,pd2,pa2,run_sarta);
[p2]            = driver_sarta_cloud_rtp(head,hattr,prof,pattr,run_sarta);

rtpsar  = [savdr 'sarta_out.rtp'];
command = [sarta_exec ' fin=' rtpout ' fout=' rtpsar ' > ' ...
                  '/home/chepplew/logs/sarta/sar_stdout.txt'];
unix(command)


% profile attribute changes for airicrad
pa = set_attr('profiles', 'robs1', infile);
pa = set_attr(pa, 'rtime', 'TAI:1958');

% Compare Obs with calcs - doo PDFs since samples may not align simply with each pair.

sbt = rad2bt(head.vchan, prof0.rcld);


btclr  = rad2bt(head.vchan, p2.rclr);
btcld  = rad2bt(head.vchan, p2.rcld);
bto    = rad2bt(head.vchan, p2.robs1);
btom   = nanmean(bto,2);
btcldm = nanmean(btcld,2);
btclrm = nanmean(btclr,2);
whos bt*

figure(1);clf;plot(head.vchan, btom,'-', head.vchan, btclrm,'-')
   hold on; plot(head.vchan,btcldm,'-')


