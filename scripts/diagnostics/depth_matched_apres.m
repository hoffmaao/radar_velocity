%DEPTH_MATCHED_APRES Compare radar and ApRES over the SAME depth interval.
%
%   WHY THIS EXISTS. Locating the ApRES sites put GA04 at 4.70 km along
%   track, where the radar admittance over its standard 0-100 m column is
%   +2.22 +/- 0.86 mm per metre of tide while ApRES reports -1.24 +/-
%   0.04: opposite signs, ~4 sigma. That looked like a hard disagreement
%   between two instruments over the same ice.
%
%   IT MAY NOT BE A DISAGREEMENT AT ALL. The two quantities are not the
%   same measurement. The radar number is the mean strain over 0-100 m.
%   The ApRES vsr is the slope of displacement against depth over the
%   interval its fit actually uses, which is 101-224 m (read from
%   fit_used_mask in the pair npz files). GA04's bed sits at 251 m, so a
%   floating plate's neutral plane is near 126 m - which puts the radar
%   interval ENTIRELY ABOVE it and about 80% of the ApRES interval BELOW
%   it. Bending strain reverses sign across the neutral plane, so two
%   instruments straddling it are EXPECTED to report opposite signs.
%
%   The lever-arm arithmetic already supports that reading: mean bending
%   strain over [a,b] scales as ((a+b)/2 - z_n), giving -75.5 m for the
%   radar and +37.0 m for ApRES, so bending predicts a ratio of -0.49
%   against the -0.56 observed - 14% apart, and the observed ratio implies
%   a neutral plane at 122 m against 126 m expected. But that is an
%   argument from a model, not a measurement.
%
%   WHAT THIS SCRIPT DOES - the decisive test. The radar measures dh(z) at
%   every depth, so it can be asked for the SAME quantity ApRES reports:
%   the mean strain over 101-224 m, instead of over 0-100 m, at GA04's own
%   position. If that comes out near -1.24, the tension dissolves and the
%   two instruments corroborate each other rather than contradicting.
%
%   It also computes the admittance in depth SLICES so the neutral plane
%   can be located from the radar alone, with no plate model assumed. A
%   sign change in that profile near 120-130 m would be the direct
%   observation the lever-arm argument only infers.
%
%   CAVEATS kept in view. The radar far-end magnitude is not resolved
%   (+2.22 at 500 m bins against +0.30 in the 1 km bin containing the
%   site, see far_end_check.m), so the ratio test inherits that
%   uncertainty. Pure bending is assumed - a depth-uniform membrane strain
%   would offset both intervals. And 224 m is inside the trusted column
%   (the coherent ice ends near 295 m and fits past ~250 m are known to
%   produce a spurious sign reversal), so the interval is usable but its
%   deep end is the weakest part of it.
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../depth_matched_apres.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef

root    = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
mp_dir  = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
net_dir = fullfile(root,'CSARP_vvel_net');
gis_dir = '/kucresis/scratch/hoffmana_sta/vvel/gis';
PASS_NAMES = {'EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};
MAX_BASELINE = 10;

% ApRES facts for GA04, from the pair npz files (fit_used_mask, bed_depth_m).
% CONVENTION: APRES_MM is the vsr slope over 101-224 m quoted in this
% project's 100 m convention - mean strain x 100 m, mm per metre of tide -
% the same convention the radar's standard 0-100 m number uses.
AP_TOP = 101; AP_BOT = 224; AP_BED = 251;
APRES_MM = -1.24; APRES_SE = 0.04;
GA04_XY = [276157.9055130380, -1318017.617914339];   % EPSG:3031 m

% Depth slices for the profile. Kept to the trusted column: the coherent
% ice ends near 295 m and fits past ~250 m are known to reverse sign
% spuriously, so nothing below 250 m is asked for.
SLICES = [ 20 60; 60 100; 100 140; 140 180; 180 224 ];

% Surface-referenced ladder. adm over [0, z] IS dh(z) in mm per metre of
% tide, so this traces the displacement profile directly instead of
% reconstructing it by differencing slices. It also separates two things
% the slices cannot: a CONSTANT dh offset (a surface-reference or residual
% coalignment error) shows up in every [0,z] row but cancels out of every
% slice, so a ladder that is large and flat while the slices are small is
% the signature of an offset rather than of strain.
LADDER = (40:40:200).';

%% Load every line at both intervals plus the slices
INTERVALS = [0 100; AP_TOP AP_BOT; SLICES; [zeros(size(LADDER)) LADDER]];
LBL = {'radar standard 0-100 m', 'ApRES-matched 101-224 m'};
for k = 1:size(SLICES,1)
  LBL{end+1} = sprintf('slice %.0f-%.0f m', SLICES(k,1), SLICES(k,2));
end
for k = 1:numel(LADDER)
  LBL{end+1} = sprintf('dh at %.0f m', LADDER(k));
end

R = [];
for n = 1:numel(PASS_NAMES)
  S = one_line(PASS_NAMES{n}, net_dir, mp_dir, MAX_BASELINE, INTERVALS);
  if isempty(S), continue; end
  if isempty(R), R = S; else, R(end+1) = S; end %#ok<AGROW>
end
assert(~isempty(R), 'no lines loaded');

%% How many blocks each interval could actually be measured over
fprintf('\n=== per-interval block coverage ===\n');
for j = 1:numel(LBL)
  nb = 0;
  for i = 1:numel(R), nb = nb + nnz(isfinite(R(i).adm(:,j))); end
  fprintf('%-26s %4d blocks across %d lines\n', LBL{j}, nb, numel(R));
end

%% Where is GA04 on the track
ps = projcrs(3031);
[along_m, off_m] = site_along(R, ps, GA04_XY(1), GA04_XY(2));
fprintf('\nGA04 at %.2f km along, %.2f km off the line\n', along_m/1e3, off_m/1e3);

%% THE TEST: the same place, two depth intervals
% ONE quantity throughout this comparison: mean strain over the interval,
% quoted in the 100 m convention (strain x 100 m, mm per metre of tide).
% one_line returns thickness change over the interval - strain times the
% interval's own thickness - so each radar row is rescaled by
% 100/(interval thickness) before anything is compared; APRES_MM is
% already quoted in this convention (see its definition above).
fprintf('\n=== the decisive comparison, at GA04''s position ===\n');
fprintf('(every row: mean strain x 100 m, in mm per metre of tide)\n');
fprintf('%-26s %10s %9s\n','quantity','mm/m tide','sigma');
TH = [100, AP_BOT - AP_TOP];
v = nan(1,2); e = nan(1,2);
for j = 1:2
  [v(j), e(j)] = stack_at(R, j, along_m/1e3);
  v(j) = v(j) * 100/TH(j); e(j) = e(j) * 100/TH(j);
  fprintf('%-26s %+10.2f %9.2f\n', LBL{j}, v(j), e(j));
end
fprintf('%-26s %+10.2f %9.2f\n','ApRES GA04 (101-224 m)', APRES_MM, APRES_SE);

fprintf('\nradar vs ApRES over the SAME interval (101-224 m): ');
if isfinite(v(2))
  nsig = abs(v(2)-APRES_MM)/hypot(e(2), APRES_SE);
  fprintf('%+.2f vs %+.2f -> %.1f sigma\n', v(2), APRES_MM, nsig);
  if nsig < 2
    fprintf('  => CONSISTENT. The 0-100 m disagreement was a depth-interval\n');
    fprintf('     mismatch, not an instrument conflict.\n');
  else
    fprintf('  => STILL DISCREPANT over matched depths; the depth-interval\n');
    fprintf('     explanation does NOT account for it on its own.\n');
  end
else
  fprintf('no coverage over 101-224 m at that position\n');
end
if all(isfinite(v))
  fprintf('mean-strain ratio (ApRES-matched / standard) = %+.2f; pure bending\n', v(2)/v(1));
  fprintf('  about a neutral plane at %.0f m predicts %+.2f (both sides are\n', ...
    AP_BED/2, ((AP_TOP+AP_BOT)/2 - AP_BED/2) / ((0+100)/2 - AP_BED/2));
  fprintf('  mean-strain ratios, so the lever arms compare like with like)\n');
end

%% The profile: locate the neutral plane from the radar alone
fprintf('\n=== admittance by depth slice at GA04''s position ===\n');
fprintf('(a sign change locates the neutral plane with NO plate model)\n');
fprintf('%-16s %10s %9s\n','slice','mm/m','sigma');
prev = NaN; zc = NaN;
for k = 1:size(SLICES,1)
  [vv, ee] = stack_at(R, 2+k, along_m/1e3);
  fprintf('%-16s %+10.2f %9.2f\n', LBL{2+k}, vv, ee);
  if isfinite(prev) && isfinite(vv) && sign(prev) ~= sign(vv)
    zc = mean([mean(SLICES(k-1,:)) mean(SLICES(k,:))]);
  end
  prev = vv;
end
if isfinite(zc)
  fprintf('sign change between slice midpoints, near %.0f m depth\n', zc);
  fprintf('  (mid-depth of a %.0f m column is %.0f m)\n', AP_BED, AP_BED/2);
else
  fprintf('no sign change across the slices covered\n');
end

%% Surface-referenced displacement profile dh(z)
fprintf('\n=== dh(z), surface-referenced, at GA04''s position ===\n');
fprintf('(flat and large here while the slices are small => a constant\n');
fprintf(' offset near the surface reference, not strain)\n');
fprintf('%-16s %10s %9s\n','depth','dh mm/m','sigma');
j0 = 2 + size(SLICES,1);
for k = 1:numel(LADDER)
  [vv, ee] = stack_at(R, j0+k, along_m/1e3);
  fprintf('%-16s %+10.2f %9.2f\n', LBL{j0+k}, vv, ee);
end

%% ========================================================================
function [v, e] = stack_at(R, j, along_km)
% Inverse-variance stack of interval j in 500 m bins, read at along_km.
xmax = 0;
for i = 1:numel(R), xmax = max(xmax, max(R(i).along)/1e3); end
edges = 0:0.5:ceil(xmax*2)/2; xc = edges(1:end-1) + 0.25;
sm = nan(size(xc)); ss = nan(size(xc));
for k = 1:numel(xc)
  a = []; q = [];
  for i = 1:numel(R)
    al = R(i).along(:)/1e3; av = R(i).adm(:,j); as_ = R(i).adm_std(:,j);
    inb = al >= edges(k) & al < edges(k+1) & isfinite(av) & isfinite(as_) & as_ > 0;
    a = [a av(inb).']; q = [q 1./as_(inb).'.^2]; %#ok<AGROW>
  end
  if numel(a) >= 2, sm(k) = sum(q.*a)/sum(q); ss(k) = sqrt(1/sum(q)); end
end
v = NaN; e = NaN;
ok = isfinite(sm);
if ~any(ok), return; end                 % this interval has no coverage
xo = xc(ok); so = sm(ok); eo = ss(ok);
% the containing bin if the point is past the outermost centre, else interpolate
if along_km > max(xo) || along_km < min(xo)
  [~, kk] = min(abs(xo - along_km));
  v = so(kk); e = eo(kk);
else
  v = interp1(xo, so, along_km, 'linear', NaN);
  e = interp1(xo, eo, along_km, 'linear', NaN);
end
end

%% ========================================================================
function [along_m, off_m] = site_along(R, ps, X, Y)
% Along-track position of an EPSG:3031 point, projected onto the block track.
along_m = NaN; off_m = NaN; best = inf;
for i = 1:numel(R)
  [bx, by] = projfwd(ps, R(i).lat, R(i).lon);
  [~, k] = min(hypot(bx - X, by - Y));
  for kk = [k-1, k+1]
    if kk < 1 || kk > numel(bx), continue; end
    vx = bx(kk)-bx(k); vy = by(kk)-by(k); L2 = vx^2 + vy^2;
    if L2 == 0, continue; end
    t = max(0, min(1, ((X-bx(k))*vx + (Y-by(k))*vy)/L2));
    d = hypot(bx(k) + t*vx - X, by(k) + t*vy - Y);
    if d < best
      best = d; off_m = d;
      along_m = R(i).along(k) + t*(R(i).along(kk) - R(i).along(k));
    end
  end
end
end

%% ========================================================================
function S = one_line(pn, net_dir, mp_dir, MAX_BASELINE, INTERVALS)
% Per-block tidal admittance for EVERY depth interval in INTERVALS.
%
% The strain over [a,b] is (dh(b) - dh(a))/(b-a) - the slope of
% displacement against depth, which is the quantity ApRES's vsr reports.
% For the 0-100 m row a = 0, so it reduces to dh(100)/100, the standard
% radar column strain, and the two are directly comparable.
S = [];
f = dir(fullfile(net_dir, [pn '_vvel_*.mat']));
keep = ~cellfun('isempty', regexp({f.name}, ...
  ['^' regexptranslate('escape',pn) '_vvel_\d+_\d+\.mat$'], 'once'));
f = f(keep);
if isempty(f), return; end
L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass');
Np = numel(L.pass); elev = nan(1,Np); ptime = nan(1,Np);
for k = 1:Np
  elev(k)  = mean(L.pass(k).elev,'omitnan');
  ptime(k) = mean(L.pass(k).gps_time,'omitnan');
end
NI = size(INTERVALS,1);

P = []; D = []; W = []; along = []; lat = []; lon = []; Nblk = 0;
for q = 1:numel(f)
  tok = regexp(f(q).name, '_vvel_(\d+)_(\d+)\.mat$', 'tokens','once');
  o = load(fullfile(net_dir, f(q).name));
  if isfield(o,'coalign_applied') && ~o.coalign_applied, continue; end
  if max(abs(o.baseline_y)) > MAX_BASELINE, continue; end
  if Nblk == 0
    Nblk = numel(o.S1); along = o.Along_track(:);
    lat = o.Latitude(:); lon = o.Longitude(:);
  end
  sv = nan(Nblk, NI);
  for b = 1:Nblk
    d = o.depth_blk(:,b); okd = isfinite(d) & isfinite(o.dh_blk(:,b));
    if nnz(okd) < 4, continue; end
    dd = d(okd); hh = o.dh_blk(okd,b);
    for j = 1:NI
      a = INTERVALS(j,1); bb = INTERVALS(j,2);
      if max(dd) < bb, continue; end
      if a == 0
        % dh is referenced to the SURFACE, so dh(0) = 0 by construction and
        % there is never a sample at exactly 0 (fits start near 20 m). This
        % is the standard radar column strain, dh(100)/100.
        ha = 0;
      else
        if min(dd) > a, continue; end
        ha = interp1(dd, hh, a, 'linear', NaN);
      end
      hb = interp1(dd, hh, bb, 'linear', NaN);
      if isfinite(ha) && isfinite(hb), sv(b,j) = (hb - ha)/(bb - a); end
    end
  end
  if all(~isfinite(sv(:))), continue; end
  if isempty(D), D = sv; else, D(:,:,end+1) = sv; end %#ok<AGROW>
  P(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
  W(end+1) = max(mean(o.coh_blk(:),'omitnan'),1e-3); %#ok<AGROW>
end
if size(P,1) < 10, return; end

adm = nan(Nblk, NI); adm_std = nan(Nblk, NI);
for j = 1:NI
  Dj = reshape(D(:,j,:), size(D,1), []);
  if all(~isfinite(Dj(:))), continue; end
  N = vdef.invertNetwork(P, Dj, struct('n_sigma',3,'weights',W,'n_epoch',Np));
  tday = (ptime - min(ptime))/86400; tide = elev - mean(elev);
  A = vdef.fitTideAdmittance(N.x, tday, tide);
  % mm of thickness change over the interval, per metre of tide: the
  % strain times the interval thickness, matching how the 0-100 m number
  % is quoted (strain x 100 m, in mm)
  th = INTERVALS(j,2) - INTERVALS(j,1);
  adm(:,j)     = 1e3*th*A.admittance(:);
  adm_std(:,j) = 1e3*th*A.admittance_std(:);
end
S = struct('name',pn,'along',along,'lat',lat,'lon',lon, ...
  'adm',adm,'adm_std',adm_std);
end
