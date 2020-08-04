function [fx rx] = load_cris_rad_for_ic_sno(crfiles, sno, cCn, vers)

%
% INPUTS:   cFnams:  List of CrIS granule file names with full path.
%           sno:     struct with upos indexes.
%           cCn:     struct of CrIS center track geolocation
%           vers:    string for origin of CrIS SDR data.
%
% OUTPUTS:  fx:      struct of frequencies of the 3 bands
%           rx:      struct of radiances of the 3 bands plus LW neighbours.
%
% Note: This version does not collect neighbour FOVs for LW channels.
%

addpath /home/motteler/shome/chirp_test        % read_netcdf_h5


%  [I, J, K] = ind2sub(size(cCn.Lat), sno.upos(:,1));
ps = 1; pe = 0; cCount = 0;  
tcrLW  = [];    tcrMW = [];   tcrSW = []; 
tcrLat = [];  tcrTim  = [];  tcrLon = []; cnrLat = [];
BB     = []; iCnt=0;
ccntr  = [15,16];

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
    junk   = ismember(sno.upos(:,1),[ps:pe]);
    cSams  = sno.upos(junk,1);                clear junk;
    fprintf('%d: %d %d ,',fn,ps,pe);
    if ( length(cSams) >= 1 )
      fprintf(1,'\t %d ',length(cSams));
      tcrLW   = [tcrLW,  cradLW(:,cSams-ps+1)];
      tcrMW   = [tcrMW,  cradMW(:,cSams-ps+1)];
      tcrSW   = [tcrSW,  cradSW(:,cSams-ps+1)];
      
    end            % end: if cSams

    ps = pe + 1;
    cCount = cCount + length(cSams); fprintf('%d \n',cCount);

  end             % end: for fn loop

  vLW = s.wnum_lw;
  vMW = s.wnum_mw;
  vSW = s.wnum_sw;
  nbr_rLW = [];  % permute(nbr_rLW,[3,2,1]);

% ------------------- otherwise CCAST --------------------
  otherwise

  for fn = 1:numel(crfiles)-1
    clear g2 rLW2 rLWp1;
    load(crfiles{fn});
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
%%    g2        = load(cFnams{fn+1});
%%    rLW2      = g2.rLW;
%%    rLW2(:,:,g2.L1a_err) = NaN;
%%    rLW2_atr1 = rLW2(:,:,:,1);
%%    rLWp1     = cat(4,rLW, rLW2_atr1);
%
    ncObs  = sz1 * 2 * sz3;              % 9x2x60 = 1080
    pe     = ps + ncObs -1;
    junk   = ismember(sno.upos(:,1),[ps:pe]);
    cSams  = sno.upos(junk,1);                clear junk;
    fprintf('%d: %d %d ,',fn,ps,pe);

    if ( length(cSams) >= 1 )
      fprintf(1,'\t %d ',length(cSams));
      tcrLW   = [tcrLW,  cradLW(:,cSams-ps+1)];
      tcrMW   = [tcrMW,  cradMW(:,cSams-ps+1)];
      tcrSW   = [tcrSW,  cradSW(:,cSams-ps+1)];
%%      for j = 1:length(cSams)
%%        clear trLW;
%%        % fold at granule boundary edges if needed
%%        map   = cFOR( cCn.Xtr(cSams(j)) ).d{cCn.Fov(cSams(j))};
%%        trLW(1,:) = rLWp1(:,map(1,3),map(1,2),cCn.Atr(cSams(j))+map(1,1));
%%        trLW(2,:) = rLWp1(:,map(2,3),map(2,2),cCn.Atr(cSams(j))+map(2,1));
%%        trLW(3,:) = rLWp1(:,map(3,3),map(3,2),cCn.Atr(cSams(j))+map(3,1));
%%        if(cCn.Atr(cSams(j)) > 1 | map(4,1) >= 0 )
%%         trLW(4,:) = rLWp1(:,map(4,3),map(4,2),cCn.Atr(cSams(j))+map(4,1));
%%        else
%%         trLW(4,:) = rLWp1(:,map(4,3),map(4,2),cCn.Atr(cSams(j)) );
%%        end
%%      % record time stamps for each neighbour for validation
%%        % TBD
%%        nbr_rLW(k,:,:) = trLW;
%%        k = k+1;
%%      end         % end for cSams

    end           % end: if cSams
    ps = pe + 1;
    cCount = cCount + length(cSams); fprintf('%d \n',cCount);
  end             % for fn loop
%%  nbr_rLW = permute(nbr_rLW,[3,2,1]);
  
end                % end: switch vers



%{
  for fn = 1:numel(crfiles)
    fprintf(1,'%d: ', fn);
    load(crfiles{fn});

    if(ndims(geo.Latitude) ~= 3) fprintf('ERROR: wrong ndims of geo.Lat\n'); end
    [nFOV, nFOR, nscan] = size(geo.Latitude);
    if(nFOV ~= 9 | mod(nFOR,30) ~= 0) fprintf('ERROR: granule size is wrong\n'); end

      junk = geo.FORTime;   junk(L1a_err) = NaN;
    junk   = junk(ccntr,:);
    tmpt   = iet2dnum(junk);
    tmpt2  = reshape(ones(9,1)*tmpt(:)', nFOV, 2, nscan);
    cnrTim = permute(tmpt2,[3 2 1]);
      junk = geo.Latitude;         junk(:,L1a_err) = NaN;
    cnrLat = permute(junk(:,ccntr,:),[3 2 1]);
      junk = geo.Longitude;             junk(:,L1a_err) = NaN;
    cnrLon = permute(junk(:,ccntr,:),[3 2 1]);
    %iic = find(cnrLat(:)' == sno.clat & cnrLon(:)' == sno.clon & cnrTim(:)' == sno.ctim);
    %     iic = intersect(cnrLat(:), sno.clat); 
    %     iid = intersect(cnrLon(:), sno.clon); 
    %     iie = intersect(cnrTim(:), sno.ctim);
    %     iCnt = iCnt+numel(iic);
    % disp([num2str(numel(iic)) ' ' num2str(numel(iid)) ' ' num2str(numel(iie)) ...
    %    ' ' num2str(iCnt)])
      junk = rLW;           junk(:,:,L1a_err) = NaN;
    cnrLW = permute(junk(:,:,ccntr,:),[1 4 3 2]);
%    tcrLW = cat(2, tcrLW, permute(junk(:,:,ccntr,:),[1 4 3 2]));
      junk = rMW;           junk(:,:,L1a_err) = NaN;
    cnrMW = permute(junk(:,:,ccntr,:),[1 4 3 2]);
%    tcrMW = cat(2, tcrMW, permute(junk(:,:,ccntr,:),[1 4 3 2]));
      junk = rSW;           junk(:,:,L1a_err) = NaN;
    cnrSW = permute(junk(:,:,ccntr,:),[1 4 3 2]);
%    tcrSW = cat(2, tcrSW, permute(junk(:,:,ccntr,:),[1 4 3 2]));
    
%%    cradLW = reshape(rLW(:,:,~L1a_err),int16(npLW),[],1);             % <- (717 x 16200) per gran
%%    cradMW = reshape(rMW(:,:,~L1a_err),int16(npMW),[],1);             % <- (869 x 16200) per gran
%%    cradSW = reshape(rSW(:,:,~L1a_err),int16(npSW),[],1);             % <- (637 x 16200) per gran

% Record sample size, accumulate index counter, find if SNO sample present.
    sz   = size(cnrLat);                         % sz [{45 or 60} x 2 x 9]
    ncsz = sz(1)*sz(2)*sz(3);
    pe   = ps + ncsz-1;

    junk   = ismember(sno.upos(:,1),[ps:pe]);
    cSams  = sno.upos(junk,1);                clear junk;
    fprintf('\t %d %d ,',ps,pe);

    if ( length(cSams) >= 1 )
      fprintf(1,'\t %d ',length(cSams));
      tcrLW   = [tcrLW,  cnrLW(:,cSams-ps+1)];
      tcrMW   = [tcrMW,  cnrMW(:,cSams-ps+1)];
      tcrSW   = [tcrSW,  cnrSW(:,cSams-ps+1)];

    end            % end: if cSams
    ps = pe + 1;
    cCount = cCount + length(cSams); fprintf('%d \n',cCount);

  end             % end: for fn loop

end                % end: switch vers
%}

%{
    fprintf(1,'%d %d ',ps,pe);
    iix = ismember(I, [ps:pe]);
    iiy = find(iix);
    if( length(iiy) >= 1)
      fprintf('%d ',length(iiy));
      IS  = I(iiy); JS = J(iiy); KS = K(iiy);
      for ii = 1:numel(iiy)
        inds = num2cell([IS(ii)-ps+1,JS(ii),KS(ii)]);
        BB   = [BB cnrLat( sub2ind(size(cnrLat), inds{:}) )];
    tcrLW = [tcrLW cnrLW(:, sub2ind(size(cnrLat), inds{:}) )];
    tcrMW = [tcrMW cnrMW(:, sub2ind(size(cnrLat), inds{:}) )];
    tcrSW = [tcrSW cnrSW(:, sub2ind(size(cnrLat), inds{:}) )];
      end
    end
    %[iix, iiy] = intersect(upos(:,1),[ps:pe]);  cSams = upos(iiy,1);  clear junk;
    %if ( length(cSams) >= 1 )
    %  tcrLat = [tcrLat; cnrLat(cSams-ps+1)];
    %  tcrLW  = [tcrLW, cnrLW(:,cSams-ps+1)];
    %  tcrMW  = [tcrMW, cnrMW(:,cSams-ps+1)];
    %  tcrSW  = [tcrSW, cnrSW(:,cSams-ps+1)];
    %end
    ps = pe + 1;
    cCount = cCount + length(iiy); fprintf('%d \n',cCount);  
  end      % <- end: for fn=1:numel(CrLst) 
%}

%%  fx = [vLW(:); vMW(:); vSW(:)];
%%  rx = [tcrLW; tcrMW; tcrSW];

rx.LW = tcrLW;
rx.MW = tcrMW;
rx.SW = tcrSW;
rx.nbr_rLW = []; % nbr_rLW;

%fx    = [vLW(:); vMW(:); vSW(:)];
fx.LW = vLW;
fx.MW = vMW;
fx.SW = vSW;
