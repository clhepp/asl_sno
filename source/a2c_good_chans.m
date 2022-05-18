function [ig] = a2c_good_chans(fc);

addpath /asl/packages/airs_decon/source          % seq_match
addpath /Users/Hepplewhite/asl.maya/gitLib/airs_decon/source

% ensure supplied vector is a column vector
if(~iscolumn(fc)) fc=fc'; end

% AIRS channel properties file
try
cp = load('/home/strow/Work/Airs/Chan_Prop_Files/chan_prop.2015.03.23.v9.5.3.mat');
catch
  warning('channel property files not here!');
end
try
cp = load('/Users/Hepplewhite/asl.maya/myLib/data/chan_prop.2015.03.23.v9.5.3.mat')
 catch
  warning('channel property file not here!');
end
ig.cp = cp;

mods = unique(cp.cmod);
ig.mods = mods;

% Get indices for each AIRS module
for i=1:length(mods)
   temp = strfind(cp.cmod,mods{i});
   mod(i).mi = find(not(cellfun('isempty',temp)));
end
ig.mod = mod;

% Get channel IDs (AIRS) for each type of channel comment
comments =  unique(cp.comment(:));
comments = comments(2:end);  % First one is always empty comment

clear bad
for i=1:length(comments)
   temp = strfind(cp.comment(:),comments(i));
   bad.(comments{i}) = find(not(cellfun('isempty',temp)));
end
ig.bad = bad;

% Get module starting and ending frequencies
for i=1:length(mods)
   s1(i).k = find(fc >= cp.cfreq(mod(i).mi(1)) & fc <= cp.cfreq(mod(i).mi(end)));
end
% Return array boundaries
ig.s1 = s1;

% Concatenate all CrIS channels that are inside AIRS module boundaries
% Also removes fill channels
gk = [];
for i=1:length(mod)
   gk = [gk; s1(i).k];
end

% Interpolate bad channels to input frequencies (fc)
fci = 1:length(fc);
ig.all = fci;

% Now get the CrIS channels that are outside AIRS detector array boundaries
allk = 1:length(fc);
ig.fill = setdiff(allk,gk);

ftemp = cp.cfreq(bad.Dead);
k = find(ftemp < max(fc) &  ftemp > min(fc));
% Get good CrIS indices
[i,j] = seq_match(sort(ftemp(k)),fc,0.5);
ig.dead = j;
j = setdiff(fci,j);
ig.all = intersect(ig.all,j);

ftemp = cp.cfreq(bad.SRF);
k = find(ftemp < max(fc) &  ftemp > min(fc));
[i,j] = seq_match(sort(ftemp(k)),fc,0.5);
ig.badsrf = j;
j = setdiff(fci,j);
ig.all = intersect(ig.all,j);

ftemp = cp.cfreq(bad.Noise);
k = find(ftemp < max(fc) &  ftemp > min(fc));
[i,j] = seq_match(sort(ftemp(k)),fc,0.5);
ig.noise = j;
j = setdiff(fci,j);
ig.all = intersect(ig.all,j);

ftemp = cp.cfreq(bad.Popping);
k = find(ftemp < max(fc) &  ftemp > min(fc));
[i,j] = seq_match(sort(ftemp(k)),fc,0.5);
ig.pop = j;
j = setdiff(fci,j);
nnnig.all = intersect(ig.all,j);

ftemp = cp.cfreq(bad.Spatial);
k = find(ftemp < max(fc) &  ftemp > min(fc));
[i,j] = seq_match(sort(ftemp(k)),fc,0.5);
ig.spatial = j;
j = setdiff(fci,j);
ig.all = intersect(ig.all,j);

%-------------------------- Remove edge channels -------------------
nig.all = ig.all;
for i=1:17
   if length(s1(i).k) > 0 
   % Only removing 1 channel per edge
   nig.all = setdiff(nig.all,[s1(i).k(1) s1(i).k(end)]);
   end
end
ig.edgechans = setdiff(ig.all,nig.all);
ig.all = nig.all;


