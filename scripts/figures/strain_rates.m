%STRAIN_RATES Vertical strain rates from the repeat-pass products.
%
%   What this shows, and why it is built this way. Each pair measures the
%   vertical displacement profile dh(z) between two epochs, so
%   strain(z) = dh(z)/z is the mean vertical strain of the column between
%   the surface and depth z. Fitting all the pairs of a line jointly,
%
%       strain(z,t) = a + b*t + c*tide
%
%   (vdef.fitTideAdmittance) separates the two things worth reporting:
%     b  the SECULAR vertical strain rate, b*365.25 in 1/yr - firn
%        compaction plus any dynamic thinning
%     c  the TIDE ADMITTANCE, strain per metre of tidal heave
%   Both are invariant to which pass is the reference; a plain correlation
%   of strain against tide is not.
%
%   THE FLOOR IS DRAWN BECAUSE IT IS THE POINT. EAGER_2022 and
%   EAGER_2022_GL1 are the SAME leg over the SAME ice, so their
%   disagreement measures systematic error directly, with no appeal to a
%   noise model: 83 microstrain per metre of tide rms per block
%   (scripts/figures/leg1_merge_check.m). Any admittance inside that band
%   is not a measurement. It is drawn as a shaded band rather than left to
%   the reader, because the previous version of this analysis reported a
%   flexure hinge that sat entirely inside it.
%
%   Panels:
%     (a) secular strain rate against DEPTH - the strain-rate profile
%     (b) secular strain rate of the 0-100 m column ALONG TRACK
%     (c) tide admittance of the 0-100 m column along track, over the floor
%
%   Run on the server, where the products live:
%     /opt/sw/matlab/2024b/bin/matlab -batch "VVEL_SUFFIX='_v3'; run('.../strain_rates.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef

if ~exist('VVEL_SUFFIX','var'), VVEL_SUFFIX = '_v3'; end
vvel_dir = ['/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel' VVEL_SUFFIX];
mp_dir   = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
out_dir  = '/kucresis/scratch/hoffmana_sta/vvel/figures';

PASS_NAMES = {'EAGER_2022','EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};
% Depths start at 50 m, not at the surface. strain(z) = dh(z)/z divides by
% the depth, so the shallowest bins amplify their own noise without adding
% information - at 25 m the block-to-block scatter is several times the
% whole signal and it dominates any axis it shares.
DEPTHS       = 50:25:250;   % [m] depths at which the joint fit is run
ALONG_DEPTH  = 100;         % [m] column depth driving panels (b) and (c)
MAX_BASELINE = 10;          % [m] cross-track baseline cut
MIN_PAIRS    = 6;

% The systematic floor for BOTH quantities is measured, not assumed, and is
% computed below from the two builds of leg 1: EAGER_2022 and
% EAGER_2022_GL1 are the same traverse over the same ice, so their
% per-block disagreement is systematic error with no noise model involved.
FLOOR_A = 'EAGER_2022';
FLOOR_B = 'EAGER_2022_GL1';

if ~exist(out_dir,'dir'), mkdir(out_dir); end

% Categorical palette, fixed order, never cycled. Validated: worst adjacent
% CVD dE 9.1 (protan), normal-vision floor 19.6. Three slots sit under 3:1
% against white, so identity is carried by a legend AND a distinct marker
% per line, never by colour alone, and the numbers are printed as a table.
PAL.cat = [0.165 0.471 0.839;    % #2a78d6 blue
           0.922 0.408 0.204;    % #eb6834 orange
           0.106 0.686 0.478;    % #1baf7a aqua
           0.929 0.631 0.000;    % #eda100 yellow
           0.910 0.482 0.643];   % #e87ba4 magenta
PAL.cat_mk  = {'o','s','^','d','v'};
PAL.ink      = [0.20 0.20 0.20];
PAL.ink_soft = [0.45 0.45 0.45];
PAL.band     = [0.90 0.90 0.88];

%% Measure
res = [];
for ip = 1:numel(PASS_NAMES)
  R = analyse(PASS_NAMES{ip}, vvel_dir, mp_dir, DEPTHS, MAX_BASELINE, MIN_PAIRS);
  if isempty(R), continue; end
  if isempty(res), res = R; else, res(end+1) = R; end %#ok<SAGROW>
end
assert(~isempty(res), 'no products could be analysed');

jz = find(DEPTHS == ALONG_DEPTH, 1);
assert(~isempty(jz), 'ALONG_DEPTH %g is not one of DEPTHS', ALONG_DEPTH);

%% The systematic floor, from two builds of the same leg
[RATE_FLOOR_Z, ADM_FLOOR_Z, nfl] = same_leg_floor(res, FLOOR_A, FLOOR_B, numel(DEPTHS));
RATE_FLOOR = RATE_FLOOR_Z(jz);
ADM_FLOOR  = ADM_FLOOR_Z(jz);
if isfinite(ADM_FLOOR)
  fprintf('\nSystematic floor from %s vs %s over %d common blocks:\n', ...
    FLOOR_A, FLOOR_B, nfl);
  fprintf('%8s %16s %16s\n','depth[m]','rate rms[1/yr]','adm rms[ue/m]');
  for z = 1:numel(DEPTHS)
    fprintf('%8.0f %16.2e %16.0f\n', DEPTHS(z), RATE_FLOOR_Z(z), 1e6*ADM_FLOOR_Z(z));
  end
else
  warning('could not measure the same-leg floor; the bands will not be drawn');
end

%% Report - the table that relieves the sub-3:1 contrast on three hues
fprintf('\n===== secular vertical strain rate, 0-%.0f m column =====\n', ALONG_DEPTH);
fprintf('%-16s %8s %14s %12s %10s\n','product','blocks','mean[1/yr]','se[1/yr]','|t|');
for i = 1:numel(res)
  v = res(i).rate(:,jz); s = res(i).rate_std(:,jz);
  ok = isfinite(v) & isfinite(s) & s > 0;
  if ~any(ok), continue; end
  mu = mean(v(ok)); se = std(v(ok))/sqrt(nnz(ok));
  fprintf('%-16s %8d %14.3e %12.3e %10.2f\n', res(i).name, nnz(ok), mu, se, abs(mu/se));
end

fprintf('\n===== tide admittance, 0-%.0f m column (floor %.0f ue/m) =====\n', ...
  ALONG_DEPTH, 1e6*ADM_FLOOR);
fprintf('%-16s %8s %14s %12s %10s\n','product','blocks','mean[ue/m]','se[ue/m]','vs floor');
for i = 1:numel(res)
  v = res(i).adm(:,jz);
  ok = isfinite(v);
  if ~any(ok), continue; end
  mu = mean(v(ok)); se = std(v(ok))/sqrt(nnz(ok));
  fprintf('%-16s %8d %14.1f %12.1f %10s\n', res(i).name, nnz(ok), 1e6*mu, 1e6*se, ...
    ternary(abs(mu) > ADM_FLOOR, 'ABOVE', 'inside'));
end

%% Figure
h = figure('Visible','off','Position',[100 100 1000 1160],'Color','w');
set(0,'CurrentFigure',h);   % Octave: imagesc/plot ignore 'parent' otherwise
axstyle = {'GridAlpha',0.15,'XColor',PAL.ink_soft,'YColor',PAL.ink_soft,'Box','off'};

% (a) strain-rate profile against depth
ax1 = axes('parent',h,'Position',[0.10 0.715 0.62 0.235]);
set(0,'CurrentFigure',h); hold(ax1,'on');
% The same-leg floor, depth by depth. Without it five disagreeing profiles
% read as structure; with it they read as what they are.
if any(isfinite(RATE_FLOOR_Z))
  zf = DEPTHS(isfinite(RATE_FLOOR_Z)); ff = RATE_FLOOR_Z(isfinite(RATE_FLOOR_Z));
  fill(ax1, [ff fliplr(-ff)], [zf fliplr(zf)], PAL.band, 'EdgeColor','none');
end
plot(ax1, [0 0], [0 max(DEPTHS)], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
hleg = []; lbl = {}; prof_all = [];
for i = 1:numel(res)
  col = PAL.cat(i,:); mk = PAL.cat_mk{i};
  prof = nan(1,numel(DEPTHS));
  for z = 1:numel(DEPTHS)
    v = res(i).rate(:,z); prof(z) = mean(v(isfinite(v)));
  end
  prof_all = [prof_all prof]; %#ok<AGROW>
  hleg(end+1) = plot(ax1, prof, DEPTHS, '-', 'Color', col, 'LineWidth', 2); %#ok<AGROW>
  plot(ax1, prof, DEPTHS, mk, 'MarkerSize', 7, 'MarkerFaceColor', col, ...
    'MarkerEdgeColor','w','LineWidth',1);
  lbl{end+1} = res(i).name; %#ok<AGROW>
end
grid(ax1,'on'); set(ax1, axstyle{:}); set(ax1,'YDir','reverse');
ylim(ax1,[0 max(DEPTHS)]);
% Explicit limits with padding: autoscale was clipping the shallowest
% points against the axis edge, which hides how far they actually swing
pa = prof_all(isfinite(prof_all));
if ~isempty(pa)
  lim = max(abs([pa(:); RATE_FLOOR_Z(isfinite(RATE_FLOOR_Z))'])) * 1.12;
  xlim(ax1, [-lim lim]);
end
xlabel(ax1,'Secular vertical strain rate (1/yr)','Color',PAL.ink);
ylabel(ax1,'Depth below surface (m)','Color',PAL.ink);
title(ax1,'Strain-rate profile: mean over blocks, over the same-leg floor (shaded)', ...
  'Color',PAL.ink);
lg = legend(ax1, hleg, lbl, 'Location','eastoutside','Interpreter','none');
set(lg,'TextColor',PAL.ink,'Box','off');

% (b) secular rate along track
ax2 = axes('parent',h,'Position',[0.10 0.395 0.80 0.235]);
set(0,'CurrentFigure',h); hold(ax2,'on');
xmax = 0;
for i = 1:numel(res), xmax = max(xmax, max(res(i).along)/1e3); end
if isfinite(RATE_FLOOR)
  fill(ax2, [0 xmax xmax 0], [-RATE_FLOOR -RATE_FLOOR RATE_FLOOR RATE_FLOOR], ...
    PAL.band, 'EdgeColor','none');
end
plot(ax2, [0 xmax], [0 0], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
for i = 1:numel(res)
  col = PAL.cat(i,:); mk = PAL.cat_mk{i};
  v = res(i).rate(:,jz); s = res(i).rate_std(:,jz); x = res(i).along/1e3;
  ok = isfinite(v);
  draw_errbars(ax2, x(ok), v(ok), s(ok), col);
  plot(ax2, x(ok), v(ok), mk, 'MarkerSize', 7, 'MarkerFaceColor', col, ...
    'MarkerEdgeColor','w','LineWidth',1);
end
grid(ax2,'on'); set(ax2, axstyle{:}); xlim(ax2,[0 xmax]);
xlabel(ax2,'Along track (km)','Color',PAL.ink);
ylabel(ax2, sprintf('Secular strain rate 0-%.0f m (1/yr)', ALONG_DEPTH),'Color',PAL.ink);
title(ax2, sprintf('Secular vertical strain rate, over the %.1e /yr same-leg floor (shaded); bars are 1 sigma', ...
  RATE_FLOOR), 'Color', PAL.ink);

% (c) tide admittance along track, against the measured systematic floor
ax3 = axes('parent',h,'Position',[0.10 0.075 0.80 0.235]);
set(0,'CurrentFigure',h); hold(ax3,'on');
if isfinite(ADM_FLOOR)
  fill(ax3, [0 xmax xmax 0], 1e6*[-ADM_FLOOR -ADM_FLOOR ADM_FLOOR ADM_FLOOR], ...
    PAL.band, 'EdgeColor','none');
end
plot(ax3, [0 xmax], [0 0], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
for i = 1:numel(res)
  col = PAL.cat(i,:); mk = PAL.cat_mk{i};
  v = 1e6*res(i).adm(:,jz); s = 1e6*res(i).adm_std(:,jz); x = res(i).along/1e3;
  ok = isfinite(v);
  draw_errbars(ax3, x(ok), v(ok), s(ok), col);
  plot(ax3, x(ok), v(ok), mk, 'MarkerSize', 7, 'MarkerFaceColor', col, ...
    'MarkerEdgeColor','w','LineWidth',1);
end
grid(ax3,'on'); set(ax3, axstyle{:}); xlim(ax3,[0 xmax]);
xlabel(ax3,'Along track (km)','Color',PAL.ink);
ylabel(ax3, sprintf('Tide admittance 0-%.0f m (\\mu\\epsilon/m)', ALONG_DEPTH), ...
  'Color',PAL.ink);
title(ax3, sprintf('Tidal response, over the %.0f \\mu\\epsilon/m same-leg floor (shaded)', ...
  1e6*ADM_FLOOR), 'Color', PAL.ink);
% Label the band at its own edge, clear of the markers, rather than on the
% zero line where it reads as a label for zero and collides with the data
text(ax3, 0.985*xmax, 1e6*ADM_FLOOR, 'same-leg disagreement ', ...
  'Color', PAL.ink_soft, 'FontSize', 9, ...
  'HorizontalAlignment','right', 'VerticalAlignment','bottom');

out_fn = fullfile(out_dir, sprintf('EAGER_2022_strain_rates%s.png', VVEL_SUFFIX));
print(h, out_fn, '-dpng', '-r120');
close(h);
fprintf('\nWrote %s\n', out_fn);

%% ========================================================================
function [rate_floor, adm_floor, n] = same_leg_floor(res, nameA, nameB, nz)
% Per-block rms disagreement between two builds of the SAME leg, at EVERY
% depth. Blocks are matched by along-track position rather than by index,
% because the two products need not start at the same place or hold the
% same block count.
rate_floor = nan(1,nz); adm_floor = nan(1,nz); n = 0;
ia = find(strcmp({res.name}, nameA), 1);
ib = find(strcmp({res.name}, nameB), 1);
if isempty(ia) || isempty(ib), return; end
A = res(ia); B = res(ib);
TOL = 100;   % [m] two blocks closer than this are the same piece of ice
pairs = [];
for k = 1:numel(A.along)
  [gap, m] = min(abs(B.along - A.along(k)));
  if gap > TOL, continue; end
  pairs(end+1,:) = [k m]; %#ok<AGROW>
end
if isempty(pairs), return; end
for z = 1:nz
  da = A.adm(pairs(:,1),z)  - B.adm(pairs(:,2),z);
  dr = A.rate(pairs(:,1),z) - B.rate(pairs(:,2),z);
  da = da(isfinite(da)); dr = dr(isfinite(dr));
  if ~isempty(da), adm_floor(z)  = sqrt(mean(da.^2)); end
  if ~isempty(dr), rate_floor(z) = sqrt(mean(dr.^2)); end
end
d = A.adm(pairs(:,1),1) - B.adm(pairs(:,2),1);
n = nnz(isfinite(d));
end

%% ========================================================================
function R = analyse(pass_name, vvel_dir, mp_dir, DEPTHS, MAX_BASELINE, MIN_PAIRS)
R = [];
f = dir(fullfile(vvel_dir, [pass_name '_vvel_*.mat']));
keep = ~cellfun('isempty', regexp({f.name}, ...
  ['^' regexptranslate('escape',pass_name) '_vvel_\d+_\d+\.mat$'], 'once'));
f = f(keep);
if isempty(f), warning('no products for %s', pass_name); return; end

L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pass_name)), 'pass');
Np = numel(L.pass); elev = nan(1,Np);
for k = 1:Np, elev(k) = mean(L.pass(k).elev,'omitnan'); end
clear L;

spy = 365.25*86400;
np = numel(f); nz = numel(DEPTHS);
strain = []; tide = nan(1,np); tday = nan(1,np); mb = nan(1,np);
along = []; lat = []; lon = []; ref0 = NaN;
for i = 1:np
  o = load(fullfile(vvel_dir, f(i).name));
  if isnan(ref0), ref0 = o.pass_idx_ref; end
  if o.pass_idx_ref ~= ref0
    error('%s: %s is referenced to pass %d but the first product uses %d; this script assumes the "main" pairing.', ...
      pass_name, f(i).name, o.pass_idx_ref, ref0);
  end
  if isempty(strain)
    Nblk = numel(o.S1);
    strain = nan(Nblk, np, nz);
    along = o.Along_track(:); lat = o.Latitude(:); lon = o.Longitude(:);
  end
  for b = 1:size(strain,1)
    d = o.depth_blk(:,b); ok = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(ok), continue; end
    for z = 1:nz
      if max(d(ok)) < DEPTHS(z), continue; end
      strain(b,i,z) = interp1(d(ok), o.dh_blk(ok,b), DEPTHS(z), 'linear', NaN) ...
        / DEPTHS(z);
    end
  end
  tide(i) = elev(o.pass_idx_sec) - elev(ref0);
  tday(i) = mean(o.GPS_time + o.delta_t_blk*spy,'omitnan')/86400;
  mb(i)   = max(abs(o.baseline_y));
end

use = ~(isfinite(mb) & mb > MAX_BASELINE);
if nnz(use) < MIN_PAIRS
  warning('%s: only %d pairs survive the baseline cut', pass_name, nnz(use));
  return;
end
strain = strain(:,use,:); tide = tide(use); tday = tday(use) - min(tday(use));

Nblk = size(strain,1);
R = struct('name', pass_name, 'along', along, 'lat', lat, 'lon', lon, ...
  'npair', nnz(use), ...
  'rate', nan(Nblk,nz), 'rate_std', nan(Nblk,nz), ...
  'adm',  nan(Nblk,nz), 'adm_std',  nan(Nblk,nz));
for z = 1:nz
  A = vdef.fitTideAdmittance(strain(:,:,z), tday, tide);
  R.rate(:,z)     = A.trend(:)     * 365.25;    % strain/day -> 1/yr
  R.rate_std(:,z) = A.trend_std(:) * 365.25;
  R.adm(:,z)      = A.admittance(:);
  R.adm_std(:,z)  = A.admittance_std(:);
end
end

%% ========================================================================
function draw_errbars(ax, x, y, s, col)
% Drawn as explicit segments rather than errorbar(): the gnuplot backend
% Octave falls back to on a headless box renders errorbar caps
% inconsistently, and this keeps the styling identical to the markers.
ok = isfinite(s);
for k = 1:numel(x)
  if ~ok(k), continue; end
  plot(ax, [x(k) x(k)], [y(k)-s(k) y(k)+s(k)], '-', 'Color', col, 'LineWidth', 1);
end
end

%% ========================================================================
function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
