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


sdate='2018/03/01';
edate='2018/03/31';
xchns=[1:713];
band='lw';
res='high';
src=1;
vers='v20a';

% extract information for file name:
sdnum     = datenum(sdate,'YYYY/mm/dd');
sdtime    = datetime(sdnum,'convertfrom','datenum');
syear     = num2str(year(sdtime));
smon      = cell2mat(month(sdtime,'shortName'));

if(src == 1) ppart1 = 'airs_cris'; end
if(src == 2) ppart1 = 'airs_cris2'; end

if(strcmp(res,'high')) ppart2 = 'HR'; cres = 'hr'; end
if(strcmp(res,'low'))  ppart2 = 'LR'; cres = 'lr'; end

dout  = strcat('/home/chepplew/data/sno/',ppart1,'/ASL/',ppart2,'/stats/');
aout  = ['sno_ac' num2str(src) '_' band '_' cres '_' syear smon];

s = load_sno_airs_cris_asl_mat(sdate,edate,xchns,res,src,vers);

r = stats_sno_airs_cris_asl_mat(s,upper(band));

disp(['Saving data to file: ' strcat(dout,aout)]);
save(strcat(dout,aout),'s','r',-v7.3');

disp('Done');


