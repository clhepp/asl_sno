function [fax ta] = load_airs_rad_for_sno(arLst, upos, aCn);

%
% INPUTS:   
%
% OUTPUTS: ta: struct. Variables: rad, Prc, Reas 
%
%
%

  fprintf('loading AIRS robs\n');
  
  ps = 1; pe = 0; aCount = 0;  k = 1;
  ta.rad  = [;]; ta.Atrk = [];  ta.Xtrk = [];  ta.sozn = []; ta.sazn  = []; 
  ta.lnfr = [];   ta.tim = [];  ta.lat  = [];  ta.lon  = []; ta.Prc   = [];
  ta.Reas = [];
  aCntr  = [43:48];
  nSNO   = size(upos,1);
  nbr_ra = zeros(nSNO,2,2645);

  for fn = 1:numel(arLst); % fn = 1;
    l1c      = read_airs_l1c(strcat(arLst(fn).folder,'/',arLst(fn).name));
    atim     = airs2dnum(l1c.Time);  % conv to matlab time.
    aatrak   = l1c.atrack;
    axtrak   = l1c.xtrack;
    tCnInd   = find(axtrak(:,aCntr)); 
    pe       = ps + numel(tCnInd) - 1;
    %  tests which data points we want to keep match to appropriate FOR
    fprintf('%d: %d %d ,',fn,ps,pe);
    junk   = ismember(upos(:,1),[ps:pe]);       aSams = upos(junk,1);    clear junk;
    if ( length(aSams >= 1)  )
      fprintf('%d \t',length(aSams));
      [nx ny nz] = size(l1c.radiances);
      % subset onto center track b4 using upos index.
      junk    = l1c.radiances(:,[aCntr],:);
      aCnRad  = permute(junk,[3,2,1]);            clear junk;
      ta.rad   = [ta.rad,  aCnRad(:,aSams-ps+1)];
      junk    = l1c.L1cProc(:,[aCntr],:);
      aCnPrc  = permute(junk,[3,2,1]);
      ta.Prc   = [ta.Prc, aCnPrc(:,aSams-ps+1)];    clear junk;
      junk    = l1c.L1cSynthReason(:,[aCntr],:);
      aCnReas = permute(junk,[3,2,1]);
      ta.Reas  = [ta.Reas, aCnReas(:,aSams-ps+1)];  clear junk;
      
      for j = 1:length(aSams)
        xxt = aSams(j)-ps+1;
        nbr_ra(k,1,:) = l1c.radiances(aCn.atr(xxt),aCn.xtr(xxt)-1,:);
        nbr_ra(k,2,:) = l1c.radiances(aCn.atr(xxt),aCn.xtr(xxt)+1,:);
        k = k+1;
      end
    end
    ps = pe+1;
    aCount = aCount + length(aSams);  fprintf('%d \n',aCount);
%   fprintf('.');
  end
  fprintf('\n');
  nbr_ra = permute(nbr_ra, [3,2,1]);              % match tarad.
  fax = [];
