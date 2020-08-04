function [iObs] = load_iasi_rad_for_ic_sno(iaLst, sno, iCn)

%
%
%
%
%

addpath /home/chepplew/gitLib/asl_sno/utilities     % readl1c_epsflip_all, iasi2mattime

fprintf('\n loading IASI obs\n');

%%[I, J] = ind2sub(size(iCn.Lat), sno.upos(:,2));
%%  whos I J 

%{
ps = 1; pe = 0; iSams = []; BD = []; tirad = []; 
iObs   = zeros(size(sno.upos,1),8461); 
iCount = 0;

for fn = 1:numel(iaLst);
    s = readl1c_epsflip_all([iaLst(fn).folder '/' ,iaLst(fn).name]);

    junk     = s.IASI_Radiances;
    iarad    = permute(junk,[3 1 2]);
    tilat    = s.Latitude;
    tiscln   = s.Scan_Line;

    tixtrk = []; tiatrk = [];
    nscans = size(tiscln,1)/30;
    if( abs(nscans - fix(size(tiscln,1)/30)) > 0.01)
      error('Problem with size of Scan_Line'); return;
    end
    for j = 1:nscans              % retain row x col structure
      is = (j-1)*30 +1;  ie = j*30;
      tiatrk(is:ie,[1 2 3 4]) = j;
      tixtrk(is:ie,[1 2 3 4]) = [1:30; 1:30; 1:30; 1:30]';    
    end

% subset on center track (FORs 15 and 16) to update counters
    ticntr   = find(tixtrk(:,1) == 15 | tixtrk(:,1) == 16);
    icnLat   = tilat(ticntr,:);
    icnrad   = iarad(:,ticntr,:);

    sz       = size(icnLat);             % [ 44 x 4]
    nisz     = sz(1)*sz(2);
    pe       = ps + nisz - 1;
    fprintf('%d  %d %d ',fn,ps,pe);

    junk   = ismember(sno.upos(:,2),[ps:pe]);
    iSams  = sno.upos(junk,2);                clear junk;

   if ( length(iSams) >= 1 )
      fprintf(1,'\t %d ',length(iSams));
      tirad   = [tirad,  icnrad(:,iSams-ps+1)];

    end            % end: if iSams
    ps = pe + 1;
    iCount = iCount + length(iSams); fprintf('%d \n',iCount);

  end             % end: for fn loop
%}

  ps = 1; pe = 0; iSams = []; iObs = zeros(size(sno.upos,1),8461); iCount = 0;

  for fn = 1:numel(iaLst);
    s = readl1c_epsflip_all([iaLst(fn).folder '/' ,iaLst(fn).name]);
    for j = 1:numel(sno.ilat)
     xv = ismember(s.Latitude,  sno.ilat(j));
     yv = ismember(s.Longitude, sno.ilon(j));
     zv = ismember(iasi2mattime(s.Time2000),  sno.itim(j));
     iwant = intersect(find(xv),intersect(find(yv),find(zv)));
       if(iwant)
         iCount = iCount + 1;
	 disp(['fn ' num2str(fn) '   iCount ' num2str(iCount) ' j ' num2str(j)]) 
         %disp(['fn  numel(iwant): ' num2str(fn) ': ' num2str(numel(iwant))]);
         [I,J] = ind2sub(size(s.Latitude),iwant(1));
         iObs(iCount,:)  = squeeze(s.IASI_Radiances(I,J,:))';
       end
     end
     fprintf('.')
   end


%{
    iix = ismember(I, [ps:pe]);
    iiy = find(iix);
    if( ~isempty(iiy))
      fprintf('%d ',length(iiy));
      IS  = I(iiy); JS = J(iiy);
      for ii = 1:numel(iiy)
        inds = num2cell([IS(ii)-ps+1,JS(ii)])
        BD   = [BD icnLat( sub2ind(size(icnLat), inds{:}) )];
        tirad = [tirad icnrad(:, sub2ind(size(icnLat), inds{:}) )];
      end
    end
    ps = pe + 1;
    iCount = iCount + length(iiy); fprintf('%d \n',iCount);  

  end                             % <- end: for fn=1:numel(iaLst) 
%}







%{
    [nx ny ~]  = size(s.IASI_Radiances);
    [mx my]    = size(s.Latitude);
    if(nx ~= mx | ny ~= my) error('IASI radiance and latitude array mis-match'); end
    t_ilat     = reshape(s.Latitude,[],1);
    t_ilon     = reshape(s.Longitude,[],1);
    t_itim     = iasi2mattime(reshape(s.Time2000,[],1));
    t_iatrak   = reshape(s.Scan_Line,[],1);
    t_ixtrak   = zeros(size(t_iatrak,1)/4,4);
    for j = 1:fix(size(t_iatrak,1)/120)              % retain row x col structure
     is = (j-1)*30 +1;  ie = j*30;
     t_ixtrak(is:ie,[1 2 3 4]) = [1:30; 1:30; 1:30; 1:30]';
    end
    t_ixtrak  = reshape(t_ixtrak,[],1);
% subset on center track (FORs 15 and 16) to update counters
    tCnInd   = find(t_ixtrak == 15 | t_ixtrak == 16);
    mx = size(s.Scan_Line,1)/30;
    xCnInd = []; for i=1:mx xCnInd(i,:)   = [15*i,16*i]; end; xCnInd = sort(xCnInd(:));
    tCnLat   = t_ialat(tCnInd);
    num_iCn  = size(s.Latitude,1)/15 * 4;
    pe = ps + length(tCnInd) - 1;
    junk = ismember(upos(:,2),[ps:pe]);   iSams = upos(junk,2);  clear junk;
    fprintf('%d\t %d %d %d\t',fn,ps,pe,length(iSams));
    if (length(iSams) >= 1)
      t_ilat   = t_ilat(tCnInd);
      t_ilon   = t_ilon(tCnInd);
      t_itim   = t_itim(tCnInd);
      x_ilat   = [x_ilat; t_ilat(iSams-ps+1)];
      x_ilon   = [x_ilon; t_ilon(iSams-ps+1)];
      x_itim   = [x_itim; t_itim(iSams-ps+1)];
      [nx ny ns] = size(s.IASI_Radiances);
      %iCnObs = s.IASI_Radiances(tCnInd,:,:);
      junk   = permute(s.IASI_Radiances(xCnInd,:,:),[1 2 3]);
      junk2  = reshape(junk(:,:,ich),[],1);
      %t_irad = junk2(tCnInd,:);
      iObs   = [iObs; junk2(iSams-ps+1,:)]; clear junk junk2 t_irad;
        end
  % Update counters:
    ps = pe+1;
    iCount = iCount + length(iSams);
    fprintf(1,'%d\n', iCount);
  end
%}

%{
 figure(1);clf;plot(sno.alat, sno.alon,'+', sno.ilat,sno.ilon,'o')
 figure(2);clf;scatter(sno.dist,sno.tdiff*24*60,[])
 %[nx nz] = size(aCnObs)
 figure(3);clf;plot([1:77],aCnObs(pos(:,1)),'.-', [1:77],iCnObs(pos(:,2)),'o-')

figure(4);clf;plot([1:17], aObs(:,ach),'.-', [1:17], iObs(:,ich),'o-')
figure(3);clf;plot([1:77], x_aCnObs(pos(:,1)),'.-', [1:77],x_iCnObs(pos(:,2)),'o-')
figure(3);clf;plot([1:17], x_aCnObs(upos(:,1)),'.-', [1:17],x_iCnObs(upos(:,2)),'o-')

%}
      
