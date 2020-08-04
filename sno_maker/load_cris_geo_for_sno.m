function c = load_cris_geo_for_sno(cFnams,vers);

%
% INPUTS:  cFnams: cell array of file names of CrIS granules w/ full path
%          vers:   string. CrIS data version (used to distinguish mat from nc files).
%
%
% OUTPUTS: c: structure containing CrIS geo data.
%
%

addpath /home/motteler/shome/chirp_test        % read_netcdf_h5
 
c.lat = [];   c.lon = [];  c.tim = [];  c.fov = [];  c.xtrak = []; 
c.atrak = []; c.gran = []; c.sozn = []; c.sazn = [];

switch vers
  case 'nasa'
  % currently no quality flagging
  for fn = 1:numel(cFnams)     % numel(crLst)
    [s, a] = read_netcdf_h5(cFnams{fn});

       junk = s.lat;   %           junk(:,L1a_err) = NaN;
    c.lat   = cat(3,c.lat, junk);
      junk  = tai2dnum(airs2tai(s.obs_time_tai93));
      tmpt  = repmat(junk,[1 1 9]);
      tmpt  = permute(tmpt,[3 1 2]);    clear junk;
    c.tim   = cat(3,c.tim, tmpt);       clear tmpt;
       junk = s.lon;   %           junk(:,L1a_err) = NaN;
    c.lon   = cat(3,c.lon, junk);
      junk  = s.sol_zen;      %junk(:,L1a_err) = NaN;
    c.sozn  = cat(3,c.sozn, junk);
      junk  = s.sat_zen;      %junk(:,L1a_err) = NaN;
    c.sazn  = cat(3,c.sazn, junk);

    % Check granule and generate FOV, a-track, x-track indexes 
    if(ndims(s.lat) ~= 3) fprintf('ERROR: wrong ndims of s.lat\n'); end
    [sz1 sz2 sz3] = size(s.lat);
    disp([num2str(fn) ': ' num2str(sz1) '  ' num2str(sz2) '  ' num2str(sz3)])
    if(sz1 ~= 9 | mod(sz2,30) ~= 0) fprintf('ERROR: granule size is wrong\n'); end
    % granule size is okay:- proceed;
    count = 0; clear tmpv tmpx tmpa;
    for m=1:sz3
      for i=1:sz2 
        for j=1:sz1
          tmpv(j, i, m) = j;          % <- nominal (9 x 30 x 60)
          tmpx(j, i, m) = i;
          tmpa(j, i, m) = m;
        end
      end
    end
    %tmpx(:,L1a_err) = NaN;
    %tmpz(:,L1a_err) = NaN;
    %tmpv(:,L1a_err) = NaN;
    c.xtrak   = cat(3,c.xtrak, tmpx);
    c.atrak   = cat(3,c.atrak, tmpa);
    c.fov     = cat(3,c.fov,   tmpv);

    tmpg = str2num(a.product_name_granule_number(2:4));
    junk = repmat(tmpg,sz1,sz3,sz2); 
    junk = permute(junk,[1 3 2]);
    %junk(:,L1a_err) = NaN; 
    c.gran = cat(3,c.gran, junk);  clear junk;
  end
  
  % ------- all CCAST data versions -----  
  otherwise

  for fn = 1:numel(cFnams)     % numel(crLst)
  load(cFnams{fn});
% Check granule and generate FOV, a-track, x-track indexes - only works with fixed array size
  if(ndims(geo.Latitude) ~= 3) fprintf('ERROR: wrong ndims of geo.Lat\n'); end
  [sz1 sz2 sz3] = size(geo.Latitude);
  disp([num2str(fn) ': ' num2str(sz1) '  ' num2str(sz2) '  ' num2str(sz3)])
  if(sz1 ~= 9 | mod(sz2,30) ~= 0) fprintf('ERROR: granule size is wrong\n'); end
% granule size is okay:- proceed;
  count = 0; clear tmpv tmpx tmpa;
  for m=1:sz3
    for i=1:sz2 
      for j=1:sz1
        tmpv(j, i, m) = j;          % <- nominal (9 x 30 x 60)
        tmpx(j, i, m) = i;
        tmpa(j, i, m) = m;
      end
    end
  end
  tmpx(:,L1a_err) = NaN;
  tmpz(:,L1a_err) = NaN;
  tmpv(:,L1a_err) = NaN;
  c.xtrak   = cat(3,c.xtrak, tmpx);
  c.atrak   = cat(3,c.atrak, tmpa);
  c.fov     = cat(3,c.fov,   tmpv);
  
% Get granule ID, should return one value per aTrack (60), but deal with xxxxxxxxxxxxx
  clear granID junk B C
  k = 1;
  for i=1:length(geo.MidTime)  % length(geo.Granule_ID)
    %if(strcmp(geo.Granule_ID(i,1:3),'NPP'))
      %junk(k) = str2num(geo.Granule_ID(i,4:end));
      %junk(k) = str2num(cFnams{fn}(98:100));
      junk(k) = str2num(cell2mat(regexp(cFnams{fn}, '(?<=_d)\d+(?=_t)','match')));
      %junk(k) = str2num(cell2mat(regexp(cFnams{fn}, '(?<=_g)\d+(?=_v20a)','match')));
      k = k + 1;
    %end
  end
% Expand to full size - take first valid value.
  B = repmat(junk(1),sz1,sz3,sz2); 
  C = permute(B,[1 3 2]);
  C(:,L1a_err) = NaN; 
  c.gran = cat(3,c.gran, C);

  errFlg = reshape(L1a_err,[],1);                     % (30 x nscan) 1: bad, 0: good.
   junk  = geo.FORTime;   junk(L1a_err) = NaN;
  tmpt   = iet2dnum(junk);    clear junk;
  % Need to expand ctim to match every FOV so that sub-setting works
  junk   = repmat(tmpt,[1 1 9]);
  tmpt   = permute(junk,[3 1 2]);    clear junk;
%  c.tim   = [c.tim; reshape(junk,[],1)]; clear junk;           % <- (16200 x 1) all FOVs
  c.tim   = cat(3,c.tim, tmpt);
    junk = geo.Latitude;              junk(:,L1a_err) = NaN;
  c.lat   = cat(3,c.lat, junk);
    junk = geo.Longitude;             junk(:,L1a_err) = NaN;
  c.lon   = cat(3,c.lon, junk);
    junk = geo.SolarZenithAngle;      junk(:,L1a_err) = NaN;
  c.sozn  = cat(3,c.sozn, junk);
    junk = geo.SatelliteZenithAngle;  junk(:,L1a_err) = NaN;
  c.sazn  = cat(3,c.sazn, junk);
  
  %  fprintf('.');
  end
  fprintf('\n');
end

%{
 % need to find bad values and setting those to NaN before can plot
 ibad = find(cCnLat < -90 | cCnLat > 90);
 cCnLat(ibad) = NaN; cCnLon(ibad) = NaN; cCnTim(ibad) = NaN; 
 figure(1);clf(1);simplemap(cCnLat(:),cCnLon(:),cCnTim(:)-cCnTim(1));
%}

