function proc_sno_airs_caf(ans)

%
% Function proc_sno_airs_caf takes AIRS L1b radiance spectrum from AIRS-CrIS or
%   AIRS-IASI SNO MAT files and data from the AIRIBRAD L1b files and produces an 
%   hdf data file suitable for use with the CAF processor.
%
% SYNOPSIS: 
%
% INPUTS:  paired Sensor abbreviation (options are: CRIS or IASI)
%
% OUTPUTS:
%
% DEPENDENCIES:
%
% NOTES:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/chepplew/gitLib/matlib/h4tools         % h4sdwrite() etc

% Check input string
if(~exist('ans','var')) fprintf('Error on input string\n'); return; end
if(~strcmp(ans,'IASI') && ~strcmp(ans,'CRIS')) fprintf('Wrong choice\n'); return; end;
if(strcmp(ans,'IASI')); IASI=1; CRIS=0; fprintf(1,'using AIRS-IASI SNOs\n'); end
if(strcmp(ans,'CRIS')); CRIS=1; IASI=0; fprintf(1,'using AIRS-CRIS SNOs\n'); end

if(IASI) 
  dp   = '/asl/s1/chepplew/projects/sno/airs_iasi/JPL/';
  fstr = 'sno_airs_iasi_*.mat';
end
if(CRIS)
  dp   = '/asl/s1/chepplew/projects/sno/airs_cris/LR/';
  fstr = 'sno_airs_crisLR_clh_2012*_018d600s.mat';
end
flst = dir(strcat(dp, fstr));
fprintf(1,'Found %d sno files to process\n',numel(flst));

for ifn = 1:numel(flst)     % 
  clear g;
  fnam   = flst(ifn).name;
  fprintf(1,'Processing %s, ifn: %d\n',fnam, ifn);
  if(IASI) sdate  = regexp(fnam,'(?<=_)[\d8]+(?=.mat)','match'); end   % for sno_airs_iasi_YYYYMMDD.mat
  if(CRIS) sdate  = regexp(fnam,'(?<=clh_)[\d.]+(?=_)','match'); end
  start_yr  = single(str2num(sdate{1}(1:4)));              % to match type for caf-6.1.0
  start_mn  = single(str2num(sdate{1}(5:6)));              % to match type for caf-6.1.0
  start_dy  = single(str2num(sdate{1}(7:8)));              % to match type for caf-6.1.0
  g      = load(strcat(dp,fnam));
  arad   = single(g.ra);
  xtrack = single(g.axtrack);
  atrack = single(g.aatrack); 
  nSam   = single(size(g.ra,2));

  %%ngran = g.atime(1);
  %%[yr mo dy hr mn sc] = datevec(ngran);
  if(IASI)                      % a month of SNOs
    [yrs mos dys hrs mns scs] = datevec(g.atime);
    [yr0 mo0 dy0 hr0 mn0 sc0] = datevec(g.atime(1)); 
    Dt = g.atime - g.atime(1);
    uDays = unique(dys);
    % find number of SNOs per day
    for i=1:numel(uDays) numDays{i} = find(dys == uDays(i)); end
    % reconstruct date string for uploading AIRIBRAD granules
    for i=1:numel(uDays) rqDates{i} = sprintf('%4d%02d%02d',yr0,mo0,uDays(i)); end
  end
  if(CRIS)                      % a single day of SNOs
    [yr0 mo0 dy0 hr0 mn0 sc0] = datevec(g.atime(1)); 
    t0    = datenum(sprintf('%4d%02d%02d',yr0,mo0,dy0),'yyyymmdd');
    Dt    = g.atime - g.atime(1);
    ngran = single(floor(Dt*239)+1);                        % granule for every SNO pair
    uGran = unique(ngran);   nGrans = single(numel(uGran)); % unique granule nums.
    uDays = unique(dy0);                                    % the day of the SNO
    rqDates{1} = sprintf('%4d%02d%02d',yr0,mo0,uDays(1));
    clear gindx; for i=1:numel(uGran) gindx{i} = find(ngran == uGran(i)); end
  end


  % Loop on each day (for CRIS there is only one day).
  % load the AIRIBRAD granule to extract the calflag, NeN, state 

  s_NeN = [;]; s_calFlag = [;]; s_exChans = [;]; s_nadTAI = [;]; s_state = [;];
  s_ngran = [];
  for jd = 1:numel(uDays)
    fprintf(1,'\t day loop: %d\t granules: ',jd);
    dateP = strcat(rqDates{jd}(1:4),'/',rqDates{jd}(5:6),'/',rqDates{jd}(7:8));
    jday  = datenum(rqDates{jd},'yyyymmdd') - datenum(rqDates{1}(1:4),'yyyy') + 1;

    % Get the granule number for IASI SNOs for this day (jd).
    if(IASI)
      jd_dy = str2num( rqDates{jd}(7:8) );
      t0 = datenum(sprintf('%4d%02d%02d',yr0,mo0,jd_dy),'yyyymmdd');
      Dt = g.atime(numDays{jd}(1):numDays{jd}(end)) - t0; 
      ngran = single(ceil(Dt*240));
      uGran = unique(ngran);   nGrans = single(numel(uGran));
    end
    disp([ngran(1),ngran(end)]);  
    
    % set the sample counters and clear cummulative variables
    ps = 1; pe = 1;
    clear state NeN calFlag exChans nadTAI;

    for nfr = 1:numel(uGran)
      nSamG = numel(find(ngran == uGran(nfr)));     % number of samples from this granule
      pe    = ps + nSamG -1;
      sGran = sprintf('%03d',uGran(nfr));
      fp3   = strcat('/asl/data/airs/AIRIBRAD/',rqDates{jd}(1:4),'/',...
          sprintf('%03d',jday),'/AIRS.',...
          rqDates{jd}(1:4),'.',rqDates{jd}(5:6),'.',rqDates{jd}(7:8),'.',sGran,...
          '.L1B.AIRS_Rad.v5.0.*.hdf');
      % some granules are missing or defective - test first
      hfile = dir(fp3);  
      if(isempty(hfile)) 
        fprintf(1,'Using previous L1b granule\n');
        hfile = dir(fp3_prev); [xpath, xnam, xext] = fileparts(fp3_prev); 
      elseif(~isempty(hfile)) 
        fp3_prev = fp3;
	hfile = dir(fp3); [xpath, xnam, xext] = fileparts(fp3); 
      end  
       
      file_id   = hdfsw('open',strcat(xpath,'/',hfile(1).name),'read');
      swath_id  = hdfsw('attach',file_id,'L1B_AIRS_Science');
      [junk,s]  = hdfsw('readfield',swath_id,'state',[],[],[]);
      if s == -1; disp('Error reading state');end;
      for j = ps:pe  state(j) = single( junk(g.axtrack(j),g.aatrack(j)) );   end  % [90x135]  
      
      [junk,s]  = hdfsw('readfield',swath_id,'NeN',[],[],[]);
      if s == -1; disp('Error reading NeN');end;
      NeN(:,nfr) = single(junk);
    
     [junk,s]   = hdfsw('readfield',swath_id,'CalFlag',[],[],[]);
      if s == -1; disp('Error reading CalFlag');end;
      for j = ps:pe calFlag(:,j) = junk(:,unique(g.aatrack(j))); end

     [junk,s]   = hdfsw('readfield',swath_id,'ExcludedChans',[],[],[]);
      if s == -1; disp('Error reading ExcludedChans');end;
      exChans(:,nfr) = junk; 

      [junk,s]=hdfsw('readfield',swath_id,'Time',[],[],[]);
      if s == -1; disp('Error reading Time');end;
      for j = ps:pe nadTAI(j) = junk(g.axtrack(j),g.aatrack(j) ); end
   
      % Close L1B granule file
      s = hdfsw('detach',swath_id);
      if s == -1; disp('Swatch detach error: L1b');end;   
      s = hdfsw('close',file_id);
      if s == -1; disp('File close error: L1b');end;

      ps = ps + nSamG;
      fprintf('.');
    end
    fprintf('\n');
    %whos exChans calFlag NeN state nadTAI nSam ngran nGrans

    s_NeN     = [s_NeN, NeN];
    s_calFlag = [s_calFlag, calFlag];
    s_exChans = [s_exChans, exChans];
    s_nadTAI  = [s_nadTAI, nadTAI];
    s_state   = [s_state, state];
    s_ngran   = [s_ngran; ngran];
  end
  whos s_NeN s_calFlag s_state s_exChans s_nadTAI s_ngran
    
  % write a temporary hdf-4 file containing the data needed by caf-6.1.0.

  hdfout = char(strcat(dp,'airs_frm_sno_preCAF_',sdate,'.hdf'));

  clear fattr slist;
  fattr = { {'author','Chris Hepplewhite'}, ...
          {'comment','intermediate data between SNO and CAF'},...
	  {'Samples', nSam},...
	  {'Granules',nGrans},... 
	  {'start_year',start_yr},...
	  {'start_month',start_mn},...
	  {'start_day',start_dy} };
  slist = { {'AirsL1b', arad},...
          {'ExcludedChans', s_exChans},...
	  {'calFlag', s_calFlag},...
	  {'NeN', s_NeN},...
	  {'State', s_state},...
	  {'nadirTAI', s_nadTAI},...
	  {'GranuleID', ngran},...
	  {'xtrack', xtrack},...
	  {'atrack', atrack} };
	  
  fprintf(1,'Writing file %s\n',hdfout);
  h4sdwrite(hdfout,slist,fattr);	 

end         % ifn MAT file number
 
%%%%%%%%% STEP 2: Run the proc_frmsno_caf executable and create post CAF files %%%%%%%
%{

%}

%%%%%%%% STEP 3: run the AIRS to CrIS translation tool
%{
%% batch_Airs2Cris.slurm; proc_airs_to_cris_hdf.m; jobs_airs_to_cris_hdf.txt

%}

%%%%%%%% STEP 4: consolidate the AIRS L1C into the original SNO MAT files %%%%%%%%

%{
if(IASI)
  % expected files:
  % /asl/s1/chepplew/projects/sno/airs_iasi/JPL/airs_frm_sno_CAF_20070301.hdf
  % /asl/s1/chepplew/projects/sno/airs_iasi/JPL/sno_airs_iasi_20070301.mat
  % /asl/s1/chepplew/projects/sno/airs_cris/LR/airs2cris_CAF_YYYYMMDD.mat

  dp   = '/asl/s1/chepplew/projects/sno/airs_iasi/JPL/';
  fstr = 'sno_airs_iasi_2011*.mat'; 
  flst = dir(strcat(dp, fstr));
  fprintf(1,'Found %d files to merge\n',numel(flst));

  for ifn=1:numel(flst);
    fnam   = flst(ifn).name; 
    fprintf(1,'Processing: %s\n',fnam); 
    sdate  = regexp(fnam,['\d+'],'match');  
    hdfin  = strcat(dp,fnam);  
    [zpath,znam,zext] = fileparts(hdfin);
    fnsno  = [zpath,'/sno_airs_iasi_',sdate{1},'.mat'];
    fncaf  = [zpath,'/airs_frm_sno_CAF_',sdate{1},'.hdf'];
    fna2c  = [zpath,'/airs2cris_CAF_',sdate{1},'.mat'];
    if(exist(fnsno) ~= 2) fprintf(1,'%s missing\n',fnsno); continue; end
    if(exist(fncaf) ~= 2) fprintf(1,'%s missing\n',fncaf); continue; end
    if(exist(fna2c) ~= 2) fprintf(1,'%s missing\n',fna2c); continue; end

    ss     = h4sdread(fncaf);                      % L1c CAF fields [Nx2645] 3 x cell arrays
    gg     = load(fna2c);                          % a2rc [1185xN], a2cfrq [1185x1]

    matOb1 = matfile(fnsno,'Writable',true);
    junk   = cell2mat(ss{1}(2))';                 % l1cRadiance {single} swap to [2645 x N]
    matOb1.ral1c = junk;      clear junk;
    junk   = cell2mat(ss{2}(2))';                 % l1cProc {uint8}
    matOb1.l1cProc = junk;    clear junk;
    junk   = cell2mat(ss{3}(2))';                 % l1cReason {uint8}
    matOb1.l1cReason = junk;  clear junk;
    junk   = single(gg.a2rc);                     % AIRS to CrIS {double complex} [1185 x N]
    matOb1.ra2c  = junk;
    junk   = single(gg.a2cfrq);                   % must retain grid from decon routines
    matOb1.fa2c  = junk;
    matOb1 = matfile(fnsno,'Writable',false);
    clear junk ss gg matOb1;

  end           % end for ifn


if(CRIS)
  % expected files:
  % /asl/s1/chepplew/projects/sno/airs_cris/LR/airs2cris_CAF_YYYYMMDD.mat
  % /asl/s1/chepplew/projects/sno/airs_cris/LR/sno_airs_crisLR_clh_YYYYMMDD_018d600s.mat
  % /asl/s1/chepplew/projects/sno/airs_cris/LR/airs2cris_CAF_YYYYMMDD.mat

  % if any one of the three are missing then skip to next date.
  
  dp   = '/asl/s1/chepplew/projects/sno/airs_cris/LR/';
  fstr = 'sno_airs_crisLR_clh_2012*_018d600s.mat'; 
  flst = dir(strcat(dp, fstr));
  fprintf(1,'Found %d files to merge\n',numel(flst));

  for ifn = 1:numel(flst)
    fnam   = flst(ifn).name;  
    sdate  = regexp(fnam,'(?<=clh_)[\d.]+(?=_)','match');  
    hdfin  = strcat(dp,fnam);  
    [zpath,znam,zext] = fileparts(hdfin);
    fnsno  = [zpath,'/sno_airs_crisLR_clh_',sdate{1},'_018d600s.mat'];
    fncaf  = [zpath,'/airs_frm_sno_CAF_',sdate{1},'.hdf'];
    fna2c  = [zpath,'/airs2cris_CAF_',sdate{1},'.mat'];
    if(exist(fnsno) ~= 2) fprintf(1,'%s missing\n',fnsno); continue; end
    if(exist(fncaf) ~= 2) fprintf(1,'%s missing\n',fncaf); continue; end
    if(exist(fna2c) ~= 2) fprintf(1,'%s missing\n',fna2c); continue; end

    ss     = h4sdread(fncaf);                      % L1c CAF fields [Nx2645]
    gg     = load(fna2c);                          % drad [1185xN], frq [1185x1]

    matOb1 = matfile(fnsno,'Writable',true);
    junk   = cell2mat(ss{1}(2));                 % l1cRadiance
    matOb1.ral1c = junk;      clear junk;
    junk   = cell2mat(ss{2}(2));                 % l1cProc
    matOb1.l1cProc = junk;    clear junk;
    junk   = cell2mat(ss{3}(2));                 % l1cReason
    matOb1.l1cReason = junk;  clear junk;
    junk   = gg.drad;
    matOb1.ra2c  = junk;
  
    matOb1 = matfile(fnsno,'Writable',false);
    clear ss gg matOb1 fnsno fncaf fna2c;
    fprintf(1,'.');
  end
fprintf(1,'\n');  
%}
