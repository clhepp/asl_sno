function run_Airs2Cris_batch()
cd /home/chepplew/projects/airs/AIRStoCRIS/batchJobs

slurmindex  = str2num(getenv('SLURM_ARRAY_TASK_ID'));

%[st, instr] = system(sprintf('sed -n "%dp" ./sno_date_list.txt | tr -d "\n"', slurmindex));
 fh = fopen('./sno_date_list.txt','r');
 fdates = textscan(fh,'%s');                 % cell array
 fclose(fh);
 njob = cell2mat(fdates{1}(slurmindex));
 proc_airs_to_cris_mat(njob);

end

njobs = 2; clear sJobs;
for i=1:njobs sJobs{i} = fLst(i).name; end
%%sJob = '/asl/s1/chepplew/projects/sno/airs_cris/LR//sno_airs_crisLR_clh_20120502_018d600s.mat';

batch = './batch_Airs2Cris.slurm';
FH = fopen(batch,'w+')
fprintf(FH,'#!/bin/bash\n\n')

fprintf(FH,'#SBATCH --job-name=Ar2Cr\n');
fprintf(FH,'#SBATCH --partition=batch\n');
fprintf(FH,'#SBATCH --qos=medium\n');
fprintf(FH,'#SBATCH --account=pi_strow\n');
fprintf(FH,'#SBATCH -N1\n');
fprintf(FH,'#SBATCH --mem-per-cpu=12000\n');
fprintf(FH,'#SBATCH --cpus-per-task 1\n');
fprintf(FH,'#SBATCH --array=1-%d\n\n',njobs);

fprintf(FH,'MATLAB=''/usr/cluster/matlab/2014b/bin/matlab''\n');
fprintf(FH,'MATOPTS='' -nodisplay -nojvm -nosplash''\n\n');

junk = sprintf('srun $MATLAB $MATOPTS -r "proc_airs_to_cris_mat(''%s''); exit"',sJob);
fprintf(FH,'%s\n',junk);

fclose(FH)
