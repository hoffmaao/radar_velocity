%LEG1_MERGE_CHECK Acceptance test for the coalignment fix: same leg, two builds.
%
%   HISTORY AND CURRENT MEANING. This comparison is what exposed the
%   tide-proportional artefact (the two builds' strain series were
%   anti-correlated at r = -0.944), and it is the regression test for
%   vdef.coalignPair: after the fix the two r(along-track) profiles
%   correlate at 0.807 (from 0.125). The residual offset in absolute r is
%   a limitation of the plain correlation itself - the builds reference
%   different epochs, and the secular strain trend aliases in with a
%   reference-dependent sign (demonstrated deterministically in
%   scripts/test_tide_admittance.m) - so the acceptance metric is now the
%   TIDE ADMITTANCE of the joint strain = a + b*t + c*tide fit
%   (vdef.fitTideAdmittance), which that aliasing cannot touch. The plain
%   r comparison is retained below for continuity with the history.
%
%   THE PASS CRITERION. EAGER_2022 and GL1-minus-20221209_01 use the SAME
%   13 passes, so their admittance profiles measure the same ice with the
%   same sampling and should agree within the fit uncertainties: profile
%   correlation high AND rms admittance difference comparable to the
%   quoted 1-sigma. GL1 with all 14 passes adds one pair, so it may
%   differ by that pair's leverage but not more.
%
%   EAGER_2022 and EAGER_2022_GL1 are the SAME out-and-back leg of the same
%   line: identical pass mid-times for every shared segment, cross-track
%   medians -213 m both. GL1's pass list is a strict superset - EAGER_2022's
%   13 segments plus 20221209_01. So the merge of the two IS GL1, and
%   EAGER_2022 is a duplicate that should be dropped rather than counted as
%   an independent line.
%
%   That also makes the fragility test cheap. The two products gave hinge
%   positions 1.8 km apart, and they differ by exactly one pass, so
%   recomputing GL1 without 20221209_01 isolates the cause: if it then
%   matches EAGER_2022, one pass was driving the difference and the
%   correlations are pair-count-limited, not measuring different ice.
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../leg1_merge_check.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef

if ~exist('VVEL_SUFFIX','var'), VVEL_SUFFIX = ''; end
vvel_dir = ['/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel' VVEL_SUFFIX];
mp_dir   = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
REF_DEPTH    = 100;
MAX_BASELINE = 10;

cases = { ...
  'EAGER_2022',     {},              'leg 1, 13 passes'; ...
  'EAGER_2022_GL1', {'20221209_01'}, 'leg 1, GL1 minus the extra pass'; ...
  'EAGER_2022_GL1', {},              'leg 1, GL1 all 14 passes (the merge)'};

R = cell(1,size(cases,1));
for c = 1:size(cases,1)
  R{c} = leg_response(cases{c,1}, vvel_dir, mp_dir, ...
    REF_DEPTH, MAX_BASELINE, cases{c,2});
  [xc, latc] = hinge(R{c}.rp, R{c}.along, R{c}.lat);
  fprintf('%-16s %-38s %2d pairs  hinge (r_partial) %s km (lat %s)\n', ...
    cases{c,1}, cases{c,3}, R{c}.npair, numstr(xc/1e3,'%.2f'), numstr(latc,'%.4f'));
end

%% The acceptance metric: do the two 13-pass builds give the same admittance?
% Blocks are matched by along-track position, not by array index: the two
% builds need not start in the same place or hold the same block count
% (same convention as strain_rates.m's same_leg_floor).
adm_compare('EAGER_2022 vs GL1-minus-20221209_01 (same 13 passes)', R{1}, R{2});
adm_compare('EAGER_2022 vs GL1 with all 14 passes', R{1}, R{3});

%% The historical plain-r comparison, kept for continuity
[i1, i2] = match_blocks(R{1}, R{2});
a = R{1}.r(i1); b = R{2}.r(i2);
ok = isfinite(a) & isfinite(b);
fprintf('\nPlain r (NOT reference-invariant, see header) over %d common blocks:\n', nnz(ok));
fprintf('  EAGER_2022 vs GL1 minus one: max |dr| = %.3f, rms = %.3f, profile corr = %.3f\n', ...
  max(abs(a(ok)-b(ok))), sqrt(mean((a(ok)-b(ok)).^2)), corr2(a(ok), b(ok)));
[j1, j3] = match_blocks(R{1}, R{3});
a3 = R{1}.r(j1); c14 = R{3}.r(j3);
ok2 = isfinite(a3) & isfinite(c14);
fprintf('  EAGER_2022 vs GL1 all 14:    max |dr| = %.3f, rms = %.3f, profile corr = %.3f\n', ...
  max(abs(a3(ok2)-c14(ok2))), sqrt(mean((a3(ok2)-c14(ok2)).^2)), corr2(a3(ok2), c14(ok2)));

% per-block table, R{1}'s grid; the other builds looked up by position
m2 = nan(numel(R{1}.adm),1); m2(i1) = i2;
m3 = nan(numel(R{1}.adm),1); m3(j1) = j3;
fprintf('\nblock  along[km] %11s %11s %11s %8s %8s %8s\n', ...
  'adm[ue/m]', 'adm[ue/m]', 'adm[ue/m]', 'rp', 'rp', 'rp');
fprintf('%17s %11s %11s %11s %8s %8s %8s\n', '', 'EAGER', 'GL1-1', 'GL1-14', ...
  'EAGER', 'GL1-1', 'GL1-14');
for k = 1:numel(R{1}.adm)
  a2v = NaN; r2v = NaN; a3v = NaN; r3v = NaN;
  if isfinite(m2(k)), a2v = R{2}.adm(m2(k)); r2v = R{2}.rp(m2(k)); end
  if isfinite(m3(k)), a3v = R{3}.adm(m3(k)); r3v = R{3}.rp(m3(k)); end
  fprintf('%5d %10.2f %11s %11s %11s %8s %8s %8s\n', k, R{1}.along(k)/1e3, ...
    numstr(1e6*R{1}.adm(k),'%.1f'), numstr(1e6*a2v,'%.1f'), ...
    numstr(1e6*a3v,'%.1f'), numstr(R{1}.rp(k),'%.2f'), ...
    numstr(r2v,'%.2f'), numstr(r3v,'%.2f'));
end

%% ========================================================================
function [ia, ib] = match_blocks(Ra, Rb)
% Match blocks between two builds by along-track position: blocks closer
% than TOL are the same piece of ice (strain_rates.m uses the same rule).
TOL = 100;   % [m]
ia = []; ib = [];
for k = 1:numel(Ra.along)
  [gap, m] = min(abs(Rb.along - Ra.along(k)));
  if gap <= TOL
    ia(end+1) = k; ib(end+1) = m; %#ok<AGROW>
  end
end
end

%% ========================================================================
function adm_compare(label, Ra, Rb)
[ia, ib] = match_blocks(Ra, Rb);
va = Ra.adm(ia); vb = Rb.adm(ib);
sa = Ra.adm_std(ia); sb = Rb.adm_std(ib);
ok = isfinite(va) & isfinite(vb);
d  = va(ok) - vb(ok);
% The quadrature sum of the two builds' quoted sigmas, averaged over the
% common blocks: the scale the rms difference should sit at if the two
% builds differ only by their noise
sig = sqrt(mean(sa(ok).^2 + sb(ok).^2));
fprintf('\n%s, %d common blocks:\n', label, nnz(ok));
fprintf('  admittance: max |diff| = %.1f ue/m, rms = %.1f ue/m (expected from fit sigma: %.1f), profile corr = %.3f\n', ...
  1e6*max(abs(d)), 1e6*sqrt(mean(d.^2)), 1e6*sig, corr2(va(ok), vb(ok)));
end

%% ========================================================================
function R = leg_response(pass_name, vvel_dir, mp_dir, ...
    REF_DEPTH, MAX_BASELINE, drop_segs)
f = dir(fullfile(vvel_dir, [pass_name '_vvel_*.mat']));
keep = false(1,numel(f));
for i = 1:numel(f)
  keep(i) = ~isempty(regexp(f(i).name, ...
    ['^' regexptranslate('escape',pass_name) '_vvel_\d+_\d+\.mat$'], 'once'));
end
f = f(keep);

L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pass_name)), 'pass');
Np = numel(L.pass);
pass_elev = nan(1,Np); pass_seg = cell(1,Np);
for k = 1:Np
  pass_elev(k) = mean(L.pass(k).elev,'omitnan');
  pass_seg{k}  = L.pass(k).param_pass.day_seg;
end
clear L;

sec_per_year = 365.25*86400;
np = numel(f);
sec_idx = nan(1,np); ref_idx = nan(1,np); maxbl = nan(1,np); t_sec = nan(1,np);
strain = []; along = []; lat = []; lon = [];
for i = 1:np
  o = load(fullfile(vvel_dir, f(i).name));
  sec_idx(i) = o.pass_idx_sec; ref_idx(i) = o.pass_idx_ref;
  if ref_idx(i) ~= ref_idx(1)
    error('%s: %s is referenced to pass %d but %s is referenced to pass %d. This script assumes the "main" pairing, where every product shares one reference pass; a sequential-pairing product must not be co-located in %s.', ...
      pass_name, f(i).name, ref_idx(i), f(1).name, ref_idx(1), vvel_dir);
  end
  maxbl(i) = max(abs(o.baseline_y));
  t_sec(i) = mean(o.GPS_time + o.delta_t_blk*sec_per_year,'omitnan');
  if isempty(strain)
    Nblk = numel(o.S1); strain = nan(Nblk, np);
    along = o.Along_track(:); lat = o.Latitude(:); lon = o.Longitude(:);
  end
  for b = 1:Nblk
    d = o.depth_blk(:,b);
    okd = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(okd) || max(d(okd)) < REF_DEPTH, continue; end
    strain(b,i) = interp1(d(okd), o.dh_blk(okd,b), REF_DEPTH, 'linear', NaN)/REF_DEPTH;
  end
end

use = ~(isfinite(maxbl) & maxbl > MAX_BASELINE);
for k = 1:numel(drop_segs)
  use = use & ~strcmp(pass_seg(sec_idx), drop_segs{k});
end
sec_idx = sec_idx(use); strain = strain(:,use); t_sec = t_sec(use);
tide   = pass_elev(sec_idx) - pass_elev(ref_idx(1));
t_days = (t_sec - min(t_sec))/86400;
npair  = numel(sec_idx);

Nblk = size(strain,1);
r = nan(Nblk,1);
for b = 1:Nblk
  s = strain(b,:);
  okb = isfinite(s) & isfinite(tide);
  if nnz(okb) < 5, continue; end
  rm = corrcoef(tide(okb), s(okb));
  r(b) = rm(1,2);
end

A = vdef.fitTideAdmittance(strain, t_days, tide);

R = struct('r',r,'adm',A.admittance(:),'adm_std',A.admittance_std(:), ...
  'rp',A.r_partial(:),'along',along,'lat',lat,'lon',lon,'npair',npair);
end

%% ========================================================================
function [xc, latc] = hinge(r, along, lat)
MIN_CONTRAST = 0.5; MIN_SIDE = 3;
xc = NaN; latc = NaN;
ok = find(isfinite(r)); n = numel(ok);
if n < 2*MIN_SIDE, return; end
rv = r(ok); xx = along(ok);
best = -inf; kb = [];
for k = MIN_SIDE:(n-MIN_SIDE)
  d = mean(rv(1:k)) - mean(rv(k+1:end));
  if d > best, best = d; kb = k; end
end
if isempty(kb) || best < MIN_CONTRAST, return; end
xc = 0.5*(xx(kb)+xx(kb+1));
latc = 0.5*(lat(ok(kb))+lat(ok(kb+1)));
end

%% ========================================================================
function s = numstr(v, fmt)
if isfinite(v), s = sprintf(fmt, v); else, s = 'none'; end
end
