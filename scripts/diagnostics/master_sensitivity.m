%MASTER_SENSITIVITY Does the choice of main (master) pass change the answer?
%
%   THE CLEANEST REPRODUCIBILITY TEST THIS DATASET ADMITS. Two builds of
%   GL3 from the SAME combine_passes input, the SAME 13 passes, the SAME
%   frozen calibration vectors - differing in exactly one processing
%   choice, baseline_master_idx (11, near the end in time, vs 6, near the
%   middle). Everything the master touches is exercised: the coregistration
%   target, the resampling grid, the surface reference, the ref_z origin.
%   Any disagreement IS master-induced systematic; nothing else differs.
%
%   Both builds are network-inverted over all 78 pairs first, so the
%   comparison is of reference-free quantities and the master enters only
%   through the PROCESSING, not through the pairing or the fit.
%
%   Blocks are matched by LATITUDE/LONGITUDE, not along-track index: the
%   two builds are resampled onto different masters' axes, which need not
%   share an origin or even a direction (the legs are walked out-and-back).
%
%   Run on the server, after both builds exist:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../master_sensitivity.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef

root   = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
mp_arch = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
mp_scr  = fullfile(root,'CSARP_multipass');
REF_DEPTH = 100; MAX_BASELINE = 10;

BUILDS = { ...
  'EAGER_2022_GL3',     fullfile(root,'CSARP_vvel_net'),    mp_arch; ...
  'EAGER_2022_GL3_m06', fullfile(root,'CSARP_vvel_netm06'), mp_scr};

S = [];
for q = 1:size(BUILDS,1)
  [pn, vdir, mdir] = BUILDS{q,:};
  B = one_build(pn, vdir, mdir, REF_DEPTH, MAX_BASELINE);
  assert(~isempty(B), 'build %s failed to load', pn);
  if isempty(S), S = B; else, S(end+1) = B; end %#ok<AGROW>
  fprintf('%-22s master %2d: %2d pairs kept, closure %.2f mm, %d blocks\n', ...
    pn, B.master, B.n_used, B.closure_mm, numel(B.lat));
end

%% Match blocks by position and compare
A = S(1); B = S(2);
da = []; ea = []; ds = []; es = []; d_m = [];
for k = 1:numel(A.lat)
  dist = hypot((B.lat - A.lat(k))*110540, ...
               (B.lon - A.lon(k))*111320*cosd(A.lat(k)));
  [dmin, m] = min(dist);
  if dmin > 100, continue; end
  d_m(end+1) = dmin; %#ok<AGROW>
  if isfinite(A.adm(k)) && isfinite(B.adm(m))
    da(end+1) = A.adm(k) - B.adm(m); %#ok<AGROW>
    ea(end+1) = A.adm_std(k)^2 + B.adm_std(m)^2; %#ok<AGROW>
  end
  if isfinite(A.sec(k)) && isfinite(B.sec(m))
    ds(end+1) = A.sec(k) - B.sec(m); %#ok<AGROW>
    es(end+1) = A.sec_std(k)^2 + B.sec_std(m)^2; %#ok<AGROW>
  end
end
assert(numel(da) >= 4, 'only %d matched blocks', numel(da));

fprintf('\n===== master 11 vs master 6, same input, same passes =====\n');
fprintf('matched blocks: %d (median position offset %.0f m)\n', numel(da), median(d_m));
fprintf('tidal admittance:  rms diff %.2f mm, expected from sigmas %.2f (%.1fx), mean diff %+.2f\n', ...
  sqrt(mean(da.^2)), sqrt(mean(ea)), sqrt(mean(da.^2))/sqrt(mean(ea)), mean(da));
fprintf('secular over win:  rms diff %.2f mm, expected from sigmas %.2f (%.1fx), mean diff %+.2f\n', ...
  sqrt(mean(ds.^2)), sqrt(mean(es)), sqrt(mean(ds.^2))/sqrt(mean(es)), mean(ds));
fprintf(['\nReading: at ~1x expected, the master choice is benign and the\n' ...
  'between-build floor must come from something else. Well above 1x, the\n' ...
  'master IS the undiagnosed systematic, and a middle master (balanced\n' ...
  'intervals, shortest mean baseline) should become the standard.\n']);

%% ========================================================================
function B = one_build(pn, vdir, mdir, REF_DEPTH, MAX_BASELINE)
B = [];
L = load(fullfile(mdir, sprintf('%s_multipass03.mat', pn)), 'pass','param_multipass');
Np = numel(L.pass); elev = nan(1,Np); ptime = nan(1,Np);
for k = 1:Np
  elev(k)  = mean(L.pass(k).elev,'omitnan');
  ptime(k) = mean(L.pass(k).gps_time,'omitnan');
end
master = L.param_multipass.multipass.baseline_master_idx;
clear L;

f = dir(fullfile(vdir, [pn '_vvel_*.mat']));
P = []; D = []; W = []; lat = []; lon = [];
for q = 1:numel(f)
  tok = regexp(f(q).name, ['^' regexptranslate('escape',pn) '_vvel_(\d+)_(\d+)\.mat$'], ...
    'tokens','once');
  if isempty(tok), continue; end
  o = load(fullfile(vdir, f(q).name));
  if isfield(o,'coalign_applied') && ~o.coalign_applied, continue; end
  if max(abs(o.baseline_y)) > MAX_BASELINE, continue; end
  Nblk = numel(o.S1); sv = nan(Nblk,1);
  for b = 1:Nblk
    d = o.depth_blk(:,b); ok = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(ok) || max(d(ok)) < REF_DEPTH, continue; end
    sv(b) = interp1(d(ok), o.dh_blk(ok,b), REF_DEPTH,'linear',NaN)/REF_DEPTH;
  end
  if all(~isfinite(sv)), continue; end
  if isempty(D)
    D = sv; lat = o.Latitude(:); lon = o.Longitude(:);
  else
    D(:,end+1) = sv; %#ok<AGROW>
  end
  P(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
  W(end+1) = max(mean(o.coh_blk(:),'omitnan'),1e-3); %#ok<AGROW>
end
if isempty(P), return; end

N = vdef.invertNetwork(P, D, struct('n_sigma',3,'weights',W));
tday = (ptime - min(ptime))/86400;
tide = elev - mean(elev);
A = vdef.fitTideAdmittance(N.x, tday, tide);
span = max(tday) - min(tday);

B = struct('name',pn,'master',master,'lat',lat,'lon',lon, ...
  'adm', 1e3*REF_DEPTH*A.admittance(:), ...
  'adm_std', 1e3*REF_DEPTH*A.admittance_std(:), ...
  'sec', 1e3*REF_DEPTH*A.trend(:)*span, ...
  'sec_std', 1e3*REF_DEPTH*A.trend_std(:)*span, ...
  'n_used', round(median(N.n_used)), ...
  'closure_mm', 1e3*REF_DEPTH*median(N.rms,'omitnan'));
end
