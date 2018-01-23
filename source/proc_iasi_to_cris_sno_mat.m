function proc_iasi_to_cris_sno_mat(CRIS, cris_res, iasi)

% Takes IASI spectral observations from {IASI CrIS} or {AIRS IASI} SNO files
%   converts them to the CrIS grid and saves them back to the original SNO
%   file.
%
% INPUTS: 
%    CRIS: Logical x 2. CrIS mission numbers NPP or J1
%          e.g. [1 0] NPP. [0 1] J1. if [0 0] then the SNO pair is AIRS:IASI
%    cris_res: Specral resolution for CrIS {'high','low'}
%    IASI: Numeric x 1. mission number 1 or 2, MetOpA or MetOpB resp.

% C Hepplewhite
%
% version 1: Nov 2015.
% 2-Jan-2018 CLH. Added logical paramaters cor CrIS and AIRS.
%

cd /home/chepplew/gitLib/asl_sno/run

addpath /home/chepplew/gitLib/asl_sno/source
addpath /home/chepplew/gitLib/asl_sno/data
addpath /asl/packages/iasi_decon
addpath /asl/packages/ccast/source                             % inst_params
addpath /asl/packages/airs_decon/source                        % hamm_app

load('/asl/data/iremis/danz/iasi_f.mat');                      % fiasi [8461x1]
xx=load('cris_freq_2grd.mat');  fcris = xx.vchan; clear xx;    % 1317 chns (12 guard chans)

% Hard code the year to process
cyear = '2018';

% Check Input Parameter for CrIS mission
if(length(CRIS) ~= 2) error('Invalid CrIS Mission Options'); return;
end
if(~all(ismember(CRIS,[0, 1]))) error('Invalid CrIS Mission selection'); return;
end

% Check input parameter for IASI mission
if(~ismember(iasi, [1,2])) error('Invalid IASI mission number'); return end
if(iasi == 1) IX = ''; end
if(iasi == 2) IX = '2'; end

% Check input Paramaeter for CrIS Resolution
cris_res = upper(cris_res);
if(~ismember(cris_res,{'HIGH','LOW'})) error('Invalid CrIS Resolution'); return;
end
if(strcmp(cris_res,'HIGH')) RES='HR'; end
if(strcmp(cris_res,'LOW'))  RES='LR'; end

% switch sources depending on CrIS Mission (or AIRS:IASI SNO)
if any(and(CRIS,[0,0]))
  dp    = '/asl/s1/chepplew/data/sno/airs_iasi/JPL/';
  fpatt = '''sno_airs_iasi_2011*.mat''';               % 3''' required to include ' in the variable
  snoLst = dir([dp 'sno_airs_iasi_' cyear '*.mat'];
end
if any(and(CRIS,[1,0]))
  CX    = '';
  dp    = ['/asl/s1/chepplew/data/sno/iasi' IX '_cris' CX '/ASL/' RES '/' cyear '/'];
  snoLst = dir([dp 'sno_iasi_cris_asl_' cyear '*.mat']);
  fpatt = '''sno_iasi_cris*.mat''';               % 3x' required to include ' in the variable
end
if any(and(CRIS,[0,1]))
  CX    = '2';
  dp    = ['/asl/s1/chepplew/data/sno/iasi' IX '_cris' CX '/ASL/' RES '/' cyear '/'];
  snoLst = dir([dp 'sno_iasi_cris_asl_' cyear '*.mat']);
  fpatt = '''sno_iasi_cris*.mat''';               % 3x' required to include ' in the variable
end
  
fprintf(1,'Found %d total SNO files\n',numel(snoLst));

opt1.hapod   = 1;             % 0: no hamming gives better result for lowres
opt1.resmode = 'hires3';      % 'lowres','hires3',...
opt1.nguard  = 2;             % 2 guard channels per band edge.

for ifn = 1:numel(snoLst)
  g = load(strcat(snoLst(ifn).folder, '/', snoLst(ifn).name));
  if(numel(g.ilat) > 1)
    if(~isfield(g, 'ri2c'))
      fprintf(1,'processing: %s\n',strcat(dp,snoLst(ifn).name)); 
      [ns nz] = size(g.ri);
      if(ns ~= 8461 & nz == 8461) xri = g.ri'; else xri = g.ri; end
      [ri2c fi2c] = iasi2cris(xri,fiasi,opt1);

      ri2c  = single(ri2c);
      matfn = matfile( strcat(dp,snoLst(ifn).name),'Writable',true );
      matfn.ri2c = ri2c;
      matfn.fi2c = fi2c;
      matfn.opts = opt1;
      matfn = matfile( strcat(dp,snoLst(ifn).name),'Writable',false );
    else
      fprintf(1,'skip %d\n',ifn)
    end 
  else
    fprintf(1,'too few pairs: %d\n',ifn)
  end
%  fprintf(1,'.');
  clear matfn ri2c fi2c g xri;
end


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


clear x ccx ccs;
unix(['cd ' dp '; find . -noleaf -maxdepth 1 -type f -name ' fpatt ' -printf ''%P\n'' > /tmp/fn.txt;']);
fh = fopen('/tmp/fn.txt');
x  = fgetl(fh);
i  = 1;
while ischar(x)
   ccx{i} = x;
   i = i + 1;
   x = fgetl(fh);
end
fclose(fh);
ccx  = unique(cellstr(ccx));
ccs = sort(ccx);
%fullfile(dp,ccs{i})  Note that i = 1:length(ccs)

%}
