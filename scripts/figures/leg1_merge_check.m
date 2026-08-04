%LEG1_MERGE_CHECK Merge the two leg-1 products, and test what separated them.
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
  [r, along, lat, lon, nseg] = leg_response(cases{c,1}, vvel_dir, mp_dir, ...
    REF_DEPTH, MAX_BASELINE, cases{c,2});
  R{c} = struct('r',r,'along',along,'lat',lat,'lon',lon);
  [xc, latc] = hinge(r, along, lat);
  fprintf('%-16s %-38s %2d pairs  hinge %s km (lat %s)\n', cases{c,1}, cases{c,3}, ...
    nseg, numstr(xc/1e3,'%.2f'), numstr(latc,'%.4f'));
end

%% Do the two leg-1 products agree once the extra pass is removed?
a = R{1}.r; b = R{2}.r;
ok = isfinite(a) & isfinite(b);
fprintf('\nEAGER_2022 vs GL1-minus-20221209_01 over %d common blocks:\n', nnz(ok));
fprintf('  max |dr| = %.3f, rms |dr| = %.3f, correlation of the two r profiles = %.3f\n', ...
  max(abs(a(ok)-b(ok))), sqrt(mean((a(ok)-b(ok)).^2)), corr2(a(ok), b(ok)));

c14 = R{3}.r;
ok2 = isfinite(a) & isfinite(c14);
fprintf('EAGER_2022 vs GL1 with all 14 passes over %d common blocks:\n', nnz(ok2));
fprintf('  max |dr| = %.3f, rms |dr| = %.3f, correlation of the two r profiles = %.3f\n', ...
  max(abs(a(ok2)-c14(ok2))), sqrt(mean((a(ok2)-c14(ok2)).^2)), corr2(a(ok2), c14(ok2)));

fprintf('\nblock  along[km]   EAGER_2022   GL1 minus one   GL1 all 14\n');
for k = 1:numel(a)
  fprintf('%5d %10.2f %12s %15s %12s\n', k, R{1}.along(k)/1e3, ...
    numstr(a(k),'%.2f'), numstr(b(k),'%.2f'), numstr(c14(k),'%.2f'));
end

%% ========================================================================
function [r, along, lat, lon, npair] = leg_response(pass_name, vvel_dir, mp_dir, ...
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
sec_idx = nan(1,np); ref_idx = nan(1,np); maxbl = nan(1,np);
strain = []; along = []; lat = []; lon = [];
for i = 1:np
  o = load(fullfile(vvel_dir, f(i).name));
  sec_idx(i) = o.pass_idx_sec; ref_idx(i) = o.pass_idx_ref;
  maxbl(i) = max(abs(o.baseline_y));
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
sec_idx = sec_idx(use); strain = strain(:,use);
tide = pass_elev(sec_idx) - pass_elev(ref_idx(1));
npair = numel(sec_idx);

Nblk = size(strain,1);
r = nan(Nblk,1);
for b = 1:Nblk
  s = strain(b,:);
  okb = isfinite(s) & isfinite(tide);
  if nnz(okb) < 5, continue; end
  rm = corrcoef(tide(okb), s(okb));
  r(b) = rm(1,2);
end
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
