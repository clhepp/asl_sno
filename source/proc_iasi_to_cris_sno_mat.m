function proc_iasi_to_cris_sno_mat()

% Takes IASI spectral observations from {IASI CrIS} or {AIRS IASI} SNO files obtained from JPL
%   converts then to the CrIS grid and saves them back to the original SNO
%   file.
%
% set CRIS and AIRS to either 0,1. to use the paired sensor with IASI.

% C Hepplewhite
%
% version 1: Nov 2015.

cd /home/chepplew/gitLib/asl_sno/run

addpath /home/chepplew/gitLib/asl_sno/source
addpath /home/chepplew/gitLib/asl_sno/data
addpath /asl/packages/iasi_decon
addpath /asl/packages/ccast/source                             % inst_params
addpath /asl/packages/airs_decon/source                        % hamm_app

load('/asl/data/iremis/danz/iasi_f.mat');                      % fiasi [8461x1]
xx=load('cris_freq_2grd.mat');  fcris = xx.vchan; clear xx;    % 1317 chns (12 guard chans)

if(CRIS)
  dp    = '/asl/s1/chepplew/projects/sno/iasi_cris/JPL/';
  fpatt = '''sno_iasi_cris*.mat''';               % 3x' required to include ' in the variable
end
if(AIRS)
  dp    = '/asl/s1/chepplew/projects/sno/airs_iasi/JPL/';
  fpatt = '''sno_airs_iasi_2011*.mat''';               % 3''' required to include ' in the variable
end

clear x cc ccs;
unix(['cd ' dp '; find . -noleaf -maxdepth 1 -type f -name ' fpatt ' -printf ''%P\n'' > /tmp/fn.txt;']);
fh = fopen('/tmp/fn.txt');
x  = fgetl(fh);
i  = 1;
while ischar(x)
   cc{i} = x;
   i = i + 1;
   x = fgetl(fh);
end
fclose(fh);
cc  = unique(cellstr(cc));
ccs = sort(cc);
%fullfile(dp,ccs{i})  Note that i = 1:length(ccs)
fprintf(1,'Found %d total SNO files\n',numel(ccs));

opt1.hapod   = 0;  % 1;                            % no hamming gives better result
opt1.resmode = 'lowres';
opt1.nguard  = 2;                                  % 2 guard channels per band edge.

for ifn = 1:numel(ccs)
  g = load(strcat(dp,ccs{ifn}));
  if(numel(g.ilat) > 100)
    %%if(~isfield(g, 'i2rc'))
      fprintf(1,'processing: %s\n',ccs{ifn}); 
      [xrad xfrq] = iasi2cris(g.ri,fiasi,opt1);

      i2rc  = single(xrad);  fi2c = single(xfrq);
      matfn = matfile( strcat(dp,ccs{ifn}),'Writable',true );
      matfn.i2rc = i2rc;
      matfn.fi2c = fi2c;
      matfn = matfile( strcat(dp,ccs{ifn}),'Writable',false );
    %%else
    %%  fprintf(1,'skip %d\n',ifn)
    %%end 
  else
    fprintf(1,'skip %d\n',ifn)
  end
%  fprintf(1,'.');
  clear matfn xrad xfrq i2rc fi2c g;
end

%clear g;

%{ 
% Sanity check - pick one file and re-load g.
clear cbt ibt dbt cbm ibm dbm
junk = single(hamm_app(double(g.rc)));
  cbt = real(rad2bt(fcris,junk));
junk = single(hamm_app(double(g.ri)));
  ibt = real(rad2bt(fiasi,junk));
junk = single(hamm_app(double(g.i2rc)));
  dbt = real(rad2bt(fcris,junk));
cbm = nanmean(cbt,2);
ibm = nanmean(ibt,2);
dbm = nanmean(dbt,2);
whos cbt ibt dbt cbm ibm dbm 

figure(1);clf;plot(fcris,cbm,'b',fiasi,ibm,'g',fcris,dbm,'r');grid on;xlim([640 1100]);
figure(1);clf;plot(fcris,cbm - dbm, 'm');grid on;axis([640 1100 -0.5 0.5])
%}
