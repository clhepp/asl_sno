function [frq tcr] = translate_cris(freq, rad, res, hapod)

% INPUTS:
%         freq.  structure. 3 Bands of frequency: (cm-1). LW, MW, SW.
%         rad.   structure. 3 bands of radiance:          LW; MW; SW.
%         res.   translation resolution. one of {'mid-res','low_res'}
%         hapod: hamming apodization. one of [0, 1];
%
% OUTPUTS:
%         frq:  Concatenated frequency after translation
%         tcr:  Concatenated radiance after translation.
%
% Based on: HM chirp_test (available on github)
%  C Hepplewhite
%  Version: draft: jun 2020

addpath /asl/packages/ccast/source           % inst_params
addpath /asl/packages/airs_decon/source      % hamm_app       

% Check inputs: freq
if(~isstruct(freq)) error('supply frequency struct with 3 members'); return; end
%if(~all(isfield(freq,{'LW','MW','SW'}))) 
if( numel(fieldnames(freq)) ~= 3 ) 
  error('Enter 3 valid members'); 
  return;
end

% Check inputs: rad
if(~isstruct(rad)) error('supply radiance struct with 3 members'); return; end
%if(~all(isfield(rad,{'LW','MW','SW'}))) 
if( numel(fieldnames(rad)) ~= 3 ) 
  error('Enter 3 valid members'); 
  return;
end
%%%strcmp(cell2mat(fieldnames(rad)),{'LW','MW','SW'})

% Check inputs: res
if(~all(ismember(res,{'midres','lowres'}))); 
   error('Supply valid user resolution');
   return;
end

% Check inputs: hapod
if(~ismember(hapod,[0 1])); error('SUpply valid hamming logic'); return; end

disp(['hapod = ', num2str(hapod)]);

% arguments for inst_params
opt2 = struct;             % inst_params opts
opt2.user_res = res;       % pass along user_res
wlaser = 773.1301;         % nominal wlaser value

%-----------------------------
% CrIS to CHIRP interpolation
%-----------------------------

% trim the LW user grid
[~, user_lw] = inst_params('LW', wlaser, opt2);
rad.LW = double(rad.LW);
if hapod, rad.LW = hamm_app(rad.LW); end
ix_lw = find(user_lw.v1 <= freq.LW & freq.LW <= user_lw.v2);
vtmp_lw = freq.LW(ix_lw);
rtmp_lw = rad.LW(ix_lw, :);
rad.LW = [];

% interpolate and trim the MW user grid
[~, user_mw] = inst_params('MW', wlaser, opt2);
[rtmp_mw, vtmp_mw] = finterp(rad.MW, freq.MW, user_mw.dv);
rtmp_mw = double(real(rtmp_mw));
if hapod, rtmp_mw = hamm_app(rtmp_mw); end
ix_mw = find(user_mw.v1 <= vtmp_mw & vtmp_mw <= user_mw.v2);
vtmp_mw = vtmp_mw(ix_mw);
rtmp_mw = rtmp_mw(ix_mw, :);
rad.MW = [];

% interpolate and trim the SW user grid
[~, user_sw] = inst_params('SW', wlaser, opt2);
[rtmp_sw, vtmp_sw] = finterp(rad.SW, freq.SW, user_sw.dv);
rtmp_sw = double(real(rtmp_sw));
if hapod, rtmp_sw = hamm_app(rtmp_sw); end
ix_sw = find(user_sw.v1 <= vtmp_sw & vtmp_sw <= user_sw.v2);
vtmp_sw = vtmp_sw(ix_sw);
rtmp_sw = rtmp_sw(ix_sw, :);
rad.SW = [];

% concatenate the bands
tcr  = [rtmp_lw; rtmp_mw; rtmp_sw];
frq  = [vtmp_lw; vtmp_mw; vtmp_sw];
clear rtmp_lw rtmp_mw rtmp_sw;


