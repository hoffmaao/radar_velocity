%MASTER_SENSITIVITY Does the choice of main (master) pass change the answer?
%
%   THE CLEANEST REPRODUCIBILITY TEST THIS DATASET ADMITS, run for EVERY
%   line: two builds from the SAME combine_passes input, the SAME passes,
%   the SAME frozen calibration vectors - differing in exactly one
%   processing choice, baseline_master_idx (the production choice vs a
%   middle-of-the-record pass). Everything the master touches is
%   exercised: the coregistration target, the resampling grid, the surface
%   reference, the ref_z origin. Any disagreement IS master-induced
%   systematic; nothing else differs.
%
%   First verified on GL3 alone (masters 11 vs 6: admittance 1.84 mm rms =
%   0.7x expected, secular 0.5x - benign); this generalisation asks the
%   same question of all five lines, since a master pathology could easily
%   afflict one product and not another.
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

% product, production master (for the report), override master. Override
% masters are the middle of each record in time; pass index tracks time
% order in these products. The m-build vvel outputs all live in one dir,
% CSARP_vvel_netm (names are unique: <product>_mNN_vvel_i_j.mat).
% EAGER_2022 is ABSENT deliberately: its archived combine_passes input
% now holds 10 passes against the 13 in the archived multipass03 product,
% so the product cannot be rebuilt from the archive at all - with any
% master. Discovered 2026-08-17; see opr_vvel/server/README.md.
LINES = { ...
  'EAGER_2022_GL1', 7; ...
  'EAGER_2022_GL2', 8; ...
  'EAGER_2022_GL3', 6; ...
  'EAGER_2022_GL4', 7};

fprintf('\n%-18s %9s %9s | %9s %9s %7s | %9s %7s\n', 'line', ...
  'masters','blocks','adm rms','adm expct','excess','sec rms','excess');
for li = 1:size(LINES,1)
  [pn, mo] = LINES{li,:};
  pn_m = sprintf('%s_m%02d', pn, mo);
  vdir_m = fullfile(root,'CSARP_vvel_netm');
  A = one_build(pn,   fullfile(root,'CSARP_vvel_net'), mp_arch, REF_DEPTH, MAX_BASELINE);
  B = one_build(pn_m, vdir_m,                          mp_scr,  REF_DEPTH, MAX_BASELINE);
  if isempty(A) || isempty(B)
    fprintf('%-18s SKIP (missing %s)\n', pn, ternary(isempty(A),'standard','override'));
    continue;
  end
  [da, ea, ds, es] = match_blocks(A, B);
  if numel(da) < 4
    fprintf('%-18s SKIP (only %d matched blocks)\n', pn, numel(da));
    continue;
  end
  fprintf('%-18s %4d/%-4d %9d | %9.2f %9.2f %6.1fx | %9.2f %6.1fx\n', ...
    pn, A.master, B.master, numel(da), ...
    sqrt(mean(da.^2)), sqrt(mean(ea)), sqrt(mean(da.^2))/sqrt(mean(ea)), ...
    sqrt(mean(ds.^2)), sqrt(mean(ds.^2))/sqrt(mean(es)));
end
fprintf(['\nReading: at ~1x expected the master choice is benign for that\n' ...
  'line. Well above 1x on any line, that line has a master-keyed\n' ...
  'systematic the others do not.\n']);

%% ========================================================================
function [da, ea, ds, es] = match_blocks(A, B)
da = []; ea = []; ds = []; es = [];
for k = 1:numel(A.lat)
  dist = hypot((B.lat - A.lat(k))*110540, ...
               (B.lon - A.lon(k))*111320*cosd(A.lat(k)));
  [dmin, m] = min(dist);
  if dmin > 100, continue; end
  if isfinite(A.adm(k)) && isfinite(B.adm(m))
    da(end+1) = A.adm(k) - B.adm(m); %#ok<AGROW>
    ea(end+1) = A.adm_std(k)^2 + B.adm_std(m)^2; %#ok<AGROW>
  end
  if isfinite(A.sec(k)) && isfinite(B.sec(m))
    ds(end+1) = A.sec(k) - B.sec(m); %#ok<AGROW>
    es(end+1) = A.sec_std(k)^2 + B.sec_std(m)^2; %#ok<AGROW>
  end
end
end

%% ========================================================================
function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end

%% ========================================================================
function B = one_build(pn, vdir, mdir, REF_DEPTH, MAX_BASELINE)
B = [];
mp_fn = fullfile(mdir, sprintf('%s_multipass03.mat', pn));
if ~exist(mp_fn,'file'), return; end     % missing build -> caller SKIPs
L = load(mp_fn, 'pass','param_multipass');
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
