cd /home/chepplew/projects/sno/airs_iasi
run /home/chepplew/myLib/matlib/paths
addpath /asl/packages/iasi_decon

sfile  = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% get IASI user-grid parameters
iasi   = iasi_params;
ifrq1  = iasi.freq;

% get AIRS channel frequencies
load('/asl/data/airs/airs_freq.mat'); fa=freq;  clear freq;

snoDir = '/asl/s1/chepplew/projects/sno/airs_iasi/ASL/';   % '.../{ASL,JPL}/';
snoLst = dir(strcat(snoDir,'sno_airs_iasi_*clh.mat'));     % or 'sno_airs_iasi_*clh.mat'

for fn = 80:numel(snoLst);
  clear i2ra afrq2 vars g;
  vars = whos('-file',strcat(snoDir,snoLst(fn).name));
  if(~ismember('i2ra', {vars.name})) 
    fprintf(1,'Working on %s\t',snoLst(fn).name);
    g = load(strcat(snoDir,snoLst(fn).name));

    % do the translation
    % convert IASI to AIRS
    i2ra   = [;];
    sz     = size(g.ri,2); fprintf(1,'%d samples\n',sz);
    lftovr = mod(sz,1000);
    rnds   = fix( (sz-lftovr)/1000);     % if samples < 5000, rnds = 0.
    tic
    ps = 1; pe = 0;
    if (rnds)
      for jj = 1:rnds
        pe = ps + 1000 - 1;
        fprintf(1,'.');
        [rr5, ff] = iasi2airs(g.ri(:,ps:pe), ifrq1, sfile, fa);
        i2ra = [i2ra, single(rr5)]; 
        ps = ps + 1000;
      end
    end
    if (lftovr >= 1)
      fprintf(1,'.');
      pe = ps + lftovr - 1;
      [rr5, ff] = iasi2airs(g.ri(:,ps:pe), ifrq1, sfile, fa);
      i2ra = [i2ra, single(rr5)];
    end
    fprintf(1,'\n');
    afrq2 = single(ff);

    % timing report
    [m,n] = size(g.ri);
    fprintf(1, 'translated %d obs in %d seconds\n', n, toc)

    matOb1        = matfile(strcat(snoDir,snoLst(fn).name),'Writable',true);
    matOb1.i2ra   = i2ra;    fprintf('.');
    matOb1.fi2a   = afrq2;   fprintf('.\n');

    matOb1        = matfile(strcat(snoDir,snoLst(fn).name),'Writable',false);
    clear matOb1;
  end   % end if(~ismember()...)

end     % for-loop fn 

%{
% some plots
arm = nanmean(gx.ra,2);                    % original AIRS
drm = nanmean(i2ra,2);                    % IASI->AIRS
irm = nanmean(gx.ri,2);                    % original IASI

abt = real(rad2bt(fa,gx.ra));       abm = nanmean(abt,2);
ibt = real(rad2bt(ifrq1,gx.ri));    ibm = nanmean(ibt,2);
dbt = real(rad2bt(afrq2,gx.i2ra));    dbm = nanmean(dbt,2);
bias = dbm - abm;
nig = load('/home/chepplew/projects/airs/master_nig_01_2009.mat');

figure(1);clf;h1=subplot(2,1,1);plot(fa,abm,'b',ifrq1,ibm,'c',afrq2,dbm,'g');
  ylabel('BT K');axis([640 2700 210 270]);grid;
  title('20150101 BT Mean AIRS l1b (b), IASI (c), I2A (g)');
  h2=subplot(2,1,2);plot(fa(nig.nig),bias(nig.nig),'m');axis([640 2700 -1 1]);grid;
  xlabel('wn cm-1');ylabel('Bias K');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');pp=get(h1,'position');
  set(h1,'position',[pp(1) pp(2)-pp(4)*0.1 pp(3) pp(4)*1.1])
  pp=get(h2,'position'); set(h2,'position',[pp(1) pp(2) pp(3) pp(4)*1.1]);
  % aslprint('./figs/201501_airs_iasi_sno_BTmeanSpec.png');
  
figure(1);clf;semilogy(afrq1,(abm-real(ar2m))./ar1m,'m');ylim([2e-5, 0.1]);
  xlabel('frequency (cm-1)');ylabel('Radiance');grid;
  title('2007-Oct Mean orig AIRS l1b - IASI->AIRS');
  % aslprint('200710_MnSpecRadDiff_Airs_Iasi2Airs.png');

figure(1);clf;h1=subplot(2,1,1);
  plot(afrq1,abm,'b',ifrq1,ibm,'g',afrq2,dbm,'c');grid;ylim([210 260]);
  h2=subplot(2,1,2);plot(afrq1,abm-dbm,'m');grid;ylim([-1 1]);

figure(4); plot(afrq1-afrq2,'r'); xlabel('Channel no.'); ylabel('wn (cm-1)');
  title('Difference of frequency scales AIRS l1b and IASI->AIRS grid');grid;
  % aslprint('Airs-Iasi2Airs_freqGrid.png');

%}
