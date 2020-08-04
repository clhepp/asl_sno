function [fx, rx] = load_cris_rad_for_sno(cFnams,upos,npLW,cCn,vers)

%
% INPUTS:   cFnams: 
%           upos: 
%           np:
%           cCn:
%           vers: 
%
% OUTPUTS:  rx:  radiance structure with 4 fields: LW, MW, SW, nbr_rLW
%           fx:  wavenumber structure with 3 fields: LW, MW, SW
%
%

addpath /home/motteler/shome/chirp_test        % read_netcdf_h5
 
fprintf('Loading CIRS SNO robs\n');
ps = 1; pe = 0; cCount = 0;   k = 1;       
tcrLW = [];   tcrMW = [];     tcrSW = [];
tcatrk = [];  tcxtrk = [];    tcFov = [];
tcTim  = [];   tcLat = [];    tcLon = [];  
tcSozn = [];  tcSazn = [];   tcLnfr = [];
snVara = [;];  cnObs = [;]; 
nSNO = size(upos,1); nbr_rLW = zeros(nSNO,4,npLW);

% load CrIS adjacent FOV look-up tables
run cris_neighbor_LUT 

switch vers
  case 'nasa'
  % currently no quality flagging
  for fn = 1:numel(cFnams)-1     % numel(crLst)
    [s, a] = read_netcdf_h5(cFnams{fn});

    % Check granule dimensions and load only nadir values
    if(ndims(s.lat) ~= 3) fprintf('ERROR: wrong ndims of s.lat\n'); end
    [sz1 sz2 sz3] = size(s.lat);

    rLW = s.rad_lw;
    %rLW(:,:,L1a_err) = NaN;  
    rMW = s.rad_mw;
    rSW = s.rad_sw;  
    cradLW = rLW(:,:,15:16,:);          % <- (717 x 9 x 2 x 45) per gran or fewer
    cradMW = rMW(:,:,15:16,:);          % <- (869 x 9 x 2 x 45) per gran or fewer
    cradSW = rSW(:,:,15:16,:);          % <- (637 x 9 x 2 x 45) per gran or fewer

 % grab next granule append the first a-track to current granule create rLWp1 etc. 
    [s2, a2]  = read_netcdf_h5(cFnams{fn+1});
    rLW2      = s2.rad_lw;
    %rLW2(:,:,g2.L1a_err) = NaN;
    rLW2_atr1 = rLW2(:,:,:,1);
    rLWp1     = cat(4,rLW, rLW2_atr1);
%
    ncObs  = sz1 * 2 * sz3;              % 9x2x45 = 810
    pe     = ps + ncObs -1;
    junk   = ismember(upos(:,2),[ps:pe]);
    cSams  = upos(junk,2);                clear junk;
    fprintf('%d: %d %d ,',fn,ps,pe);
    if ( length(cSams) >= 1 )
      fprintf(1,'\t %d ',length(cSams));
      tcrLW   = [tcrLW,  cradLW(:,cSams-ps+1)];
      tcrMW   = [tcrMW,  cradMW(:,cSams-ps+1)];
      tcrSW   = [tcrSW,  cradSW(:,cSams-ps+1)];
%{
      for j = 1:length(cSams)
        clear trLW;
        % fold at granule boundary edges if needed
        map   = cFOR( cCn.Xtr(cSams(j)) ).d{cCn.Fov(cSams(j))};
        trLW(1,:) = rLWp1(:,map(1,3),map(1,2),cCn.Atr(cSams(j))+map(1,1));
        trLW(2,:) = rLWp1(:,map(2,3),map(2,2),cCn.Atr(cSams(j))+map(2,1));
        trLW(3,:) = rLWp1(:,map(3,3),map(3,2),cCn.Atr(cSams(j))+map(3,1));
        if(cCn.Atr(cSams(j)) > 1 | map(4,1) >= 0 )
         trLW(4,:) = rLWp1(:,map(4,3),map(4,2),cCn.Atr(cSams(j))+map(4,1));
        else
         trLW(4,:) = rLWp1(:,map(4,3),map(4,2),cCn.Atr(cSams(j)) );
        end
      % record time stamps for each neighbour for validation
        % TBD
        nbr_rLW(k,:,:) = trLW;
        k = k+1;
      end          % end for cSams
%}
    end            % end: if cSams

    ps = pe + 1;
    cCount = cCount + length(cSams); fprintf('%d \n',cCount);
    
  end             % end: for fn loop

  vLW = s.wnum_lw;
  vMW = s.wnum_mw;
  vSW = s.wnum_sw;
  nbr_rLW = permute(nbr_rLW,[3,2,1]);
  
% ------------------- otherwise CCAST --------------------  
  otherwise
  
  for fn = 1:numel(cFnams)-1
    clear g2 rLW2 rLWp1;
    load(cFnams{fn});
% Check granule dimensions and load only nadir values
    if(ndims(geo.Latitude) ~= 3) fprintf('ERROR: wrong ndims of geo.Lat\n'); end
    [sz1 sz2 sz3] = size(geo.Latitude);

    rLW(:,:,L1a_err) = NaN;  
    cradLW = rLW(:,:,15:16,:);          % <- (717 x 9 x 2 x 60) per gran or fewer
    rMW(:,:,L1a_err) = NaN;
    cradMW = rMW(:,:,15:16,:);          % <- (869 x 9 x 2 x 60) per gran
    rSW(:,:,L1a_err) = NaN;
    cradSW = rSW(:,:,15:16,:);          % <- (637 x 9 x 2 x 60) per gran
 % grab next granule append the first a-track to current granule create rLWp1 etc. 
    g2        = load(cFnams{fn+1});
    rLW2      = g2.rLW;
    rLW2(:,:,g2.L1a_err) = NaN;
    rLW2_atr1 = rLW2(:,:,:,1);
    rLWp1     = cat(4,rLW, rLW2_atr1);
%
    ncObs  = sz1 * 2 * sz3;              % 9x2x60 = 1080
    pe     = ps + ncObs -1;
    junk   = ismember(upos(:,2),[ps:pe]);
    cSams  = upos(junk,2);                clear junk;
    fprintf('%d: %d %d ,',fn,ps,pe);

    if ( length(cSams) >= 1 )
      fprintf(1,'\t %d ',length(cSams));
      tcrLW   = [tcrLW,  cradLW(:,cSams-ps+1)];
      tcrMW   = [tcrMW,  cradMW(:,cSams-ps+1)];
      tcrSW   = [tcrSW,  cradSW(:,cSams-ps+1)];
      for j = 1:length(cSams)
        clear trLW;
        % fold at granule boundary edges if needed
        map   = cFOR( cCn.Xtr(cSams(j)) ).d{cCn.Fov(cSams(j))};
        trLW(1,:) = rLWp1(:,map(1,3),map(1,2),cCn.Atr(cSams(j))+map(1,1));
        trLW(2,:) = rLWp1(:,map(2,3),map(2,2),cCn.Atr(cSams(j))+map(2,1));
        trLW(3,:) = rLWp1(:,map(3,3),map(3,2),cCn.Atr(cSams(j))+map(3,1));
        if(cCn.Atr(cSams(j)) > 1 | map(4,1) >= 0 )
         trLW(4,:) = rLWp1(:,map(4,3),map(4,2),cCn.Atr(cSams(j))+map(4,1));
        else
         trLW(4,:) = rLWp1(:,map(4,3),map(4,2),cCn.Atr(cSams(j)) );
        end
      % record time stamps for each neighbour for validation
        % TBD
        nbr_rLW(k,:,:) = trLW;
        k = k+1;
      end         % end for cSams
    end           % end: if cSams
    ps = pe + 1;
    cCount = cCount + length(cSams); fprintf('%d \n',cCount);
  end             % for fn loop
  nbr_rLW = permute(nbr_rLW,[3,2,1]);
  
end                % end: switch vers

%rx    = [tcrLW; tcrMW; tcrSW];
rx.LW = tcrLW;
rx.MW = tcrMW;
rx.SW = tcrSW;
rx.nbr_rLW = nbr_rLW;

%fx    = [vLW(:); vMW(:); vSW(:)];
fx.LW = vLW;
fx.MW = vMW;
fx.SW = vSW;

% whos tcr* nbr_rLW rc fc

