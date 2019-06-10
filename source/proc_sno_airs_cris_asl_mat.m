function proc_sno_airs_cris_asl_mat.m

%
%
%
%
%
%
%

addpath /home/chepplew/gitLib/asl_sno/source

cd /home/chepplew/projects/sno/airs_cris;

spec_res = {'high','medium'};
src      = 1;
vers     = 'v20a';

date_pairs = {'2017/01/01','2017/03/31'; ...
              '2017/04/01','2017/06/30'; ...
	      '2017/07/01','2017/09/30'; ...
	      '2017/10/01','2017/12/31'};

channel_sets = {1:715; 716:1346; 1347:1683};
all_bands = {'LW','MW','SW'};

if(src == 1) ppart1 = 'airs_cris'; end
if(src == 2) ppart1 = 'airs_cris2'; end

if(strcmp(spec_res,{'high','medium'})) cr1 = 'hr'; cr2 = 'mr'; end
if(strcmp(spec_res,{'high','low'}))    cr1 = 'hr'; cr2 = 'lr'; end

k = 1;
  sdate = date_pairs{k,1};
  edate = date_pairs{k,2};;
  % extract information for file name:
  sdnum     = datenum(sdate,'YYYY/mm/dd');
  ednum     = datenum(edate,'YYYY/mm/dd');
  sdtime    = datetime(sdnum,'convertfrom','datenum');
  edtime    = datetime(ednum,'convertfrom','datenum');
  syear     = num2str(year(sdtime));
  smon      = cell2mat(month(sdtime,'shortName'));
  emon      = cell2mat(month(edtime,'shortname'));


  j = 1;

    xchns = channel_sets{j};
    band  = all_bands{j};

    dout  = ['/home/chepplew/data/sno/airs_cris/stats/'] ;
    aout  = ['sno_ac' num2str(src) '_' lower(band) '_' cr2 '_' ...
              syear lower(smon) '_' lower(emon) '.mat'];

    disp(['Doing band: ' band '. Dates: ' sdate ' to ' edate])
    s = load_sno_airs_cris_asl_mat(sdate,edate,xchns,spec_res{2},src,vers);

    r = stats_sno_airs_cris_asl_mat(s,upper(band));

    disp(['Saving data to file: ' strcat(dout,aout)]);
    save([dout aout],'s','r','-v7.3');

  % end   % j band
% end     % dates
disp('Done');


