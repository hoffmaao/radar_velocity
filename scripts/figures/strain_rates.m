%STRAIN_RATES How much does the ice column actually change thickness?
%
%   METHOD (2026-08-17): everything here comes from the NETWORK inversion
%   (vdef.invertNetwork) over ALL pairs of each product - 78-105 pairs
%   against the 12-14 a single-reference pairing used - with robust
%   rejection, then the joint x = a + b*t + c*tide fit on the per-pass
%   displacements. No reference pass exists anywhere in the chain, and
%   the same estimator produced every number and every reference line on
%   the figure. Requires the CSARP_vvel_net products (pairs='all').
%
%   EVERYTHING HERE IS IN MILLIMETRES OF COLUMN THICKNESS CHANGE, because
%   that is the thing being measured and the thing a reader can picture.
%   Earlier versions reported "tide admittance" in microstrain per metre
%   and "secular strain rate" in 1/yr, which are correct, jargon, and
%   impossible to sanity-check by eye.
%
%   The conversion is just multiplication by the column depth:
%     dh(z) = strain(z) * z
%   so over the top 100 m, 1 mm of thickness change IS 10 microstrain.
%   Negative means the column SHORTENED - the ice squashed.
%
%   WHAT IS PLOTTED
%     (a) tidal response against depth: mm the column changes thickness for
%         each metre the shelf rises on the tide, against the flexure model
%     (b) the same at one depth, along the line, against the ApRES
%         rate-method value (-1.24 +/- 0.04 mm over the top 100 m,
%         scripts/diagnostics/apres_rate_check.py)
%     (c) the secular change: mm the column changed over the observation
%         window itself, with no extrapolation to a year - a 3-day baseline
%         reported as a yearly rate multiplies both signal and error by
%         ~110 and stops being interpretable
%
%   THREE THINGS ARE DRAWN ON EVERY PANEL
%     the measurements, one series per product;
%     the FLOOR (shaded) - what the method's own systematic error is,
%       measured from two builds of the SAME leg over the SAME ice, so it
%       needs no noise model to justify;
%     the EXPECTED signal (dashed) - what the physics predicts, so the
%       reader can see the gap without doing arithmetic.
%
%   Run on the server, where the products live:
%     /opt/sw/matlab/2024b/bin/matlab -batch "VVEL_SUFFIX='_v3'; run('.../strain_rates.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef
addpath(fileparts(mfilename('fullpath')));            % grl_figure
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'diagnostics'));  % pass_tide

if ~exist('VVEL_SUFFIX','var'), VVEL_SUFFIX = '_net'; end
vvel_dir = ['/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel' VVEL_SUFFIX];
mp_dir   = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
out_dir  = vdef.figureDir();

% The four lines, plus the duplicate build of GL1. EAGER_2022 is NOT a
% fifth profile - it is a second build of GL1 sharing thirteen of its
% fourteen passes (vdef.surveyLines) - but the floor computed below needs
% two builds of the SAME ice, so it is loaded here deliberately and must
% not be read as another line.
[LINES, DUP] = vdef.surveyLines();
PASS_NAMES = [{DUP.name}, LINES];
% strain(z) = dh(z)/z divides by depth, so the shallowest bins amplify
% their own noise; start at 50 m.
DEPTHS       = 50:25:250;
ALONG_DEPTH  = 100;
MAX_BASELINE = 10;
MIN_PAIRS    = 6;
FLOOR_A = 'EAGER_2022';      % two builds of the SAME leg: their
FLOOR_B = 'EAGER_2022_GL1';  % disagreement IS the systematic error
APRES_MM = -1.24;            % mm per m of tide over the top 100 m, ApRES
APRES_SE = 0.04;             % GA04 rate method, apres_rate_check.py

% Expected TIDAL signal, from thin-plate flexure. A plate of thickness H
% bent to curvature kappa = 1/L^2 per metre of tide has surface bending
% strain eps_xx = (H/2)*kappa, falling linearly to zero at the neutral
% plane H/2 and reversing below it. The vertical strain follows by
% Poisson, eps_zz = -nu/(1-nu)*eps_xx. Averaged over the top z:
%   eps_zz_avg(z) = -nu/(1-nu) * (H/2)*kappa * (1 - z/H)
%   dh(z)         = eps_zz_avg(z) * z
% SENSITIVE TO L: the curvature goes as 1/L^2, so halving the flexure
% wavelength quadruples the expected signal. 2 km is a mid-range guess for
% this setting, not a fitted value - read the dashed curve as an order of
% magnitude, not a prediction.
H_ICE  = 300;    % [m] ice thickness
L_FLEX = 2000;   % [m] flexure length scale
NU     = 0.33;   % Poisson ratio

% Expected SECULAR signal: steady-state firn compaction. An accumulation
% rate b_dot of ice-equivalent thickness requires the firn column to
% compact at b_dot*(rho_ice/rho_surf - 1) to stay in steady state.
%
% That compaction happens in the FIRN, between the surface and bubble
% close-off - NOT spread through the full ice thickness. So a column
% reaching past close-off already contains essentially all of it, and
% deepening the column further adds no more signal. vdef.defaultParams
% puts close-off at 60 m here, so every depth plotted from 75 m down
% captures the whole of it.
ACC_ICE_EQ = 0.20;   % [m ice/yr] Windless Bight, order of magnitude
RHO_SURF   = 350; RHO_ICE = 917;
BCO_DEPTH  = vdef.defaultParams().bco_depth;   % [m] bubble close-off

if ~exist(out_dir,'dir'), mkdir(out_dir); end

% Categorical palette, fixed order, never cycled. Validated: worst adjacent
% CVD dE 9.1 (protan), normal-vision floor 19.6. Three slots sit under 3:1
% against white, so identity is carried by a legend AND a distinct marker
% per series, never by colour alone, and the numbers are printed as a table.
PAL.cat = [0.165 0.471 0.839;    % #2a78d6 blue
           0.922 0.408 0.204;    % #eb6834 orange
           0.106 0.686 0.478;    % #1baf7a aqua
           0.929 0.631 0.000;    % #eda100 yellow
           0.910 0.482 0.643];   % #e87ba4 magenta
PAL.cat_mk  = {'o','s','^','d','v'};
PAL.ink      = [0 0 0];
PAL.ink_soft = [0.45 0.45 0.45];
PAL.band     = [0.90 0.90 0.88];
PAL.expect   = [0.35 0.35 0.35];

%% Measure
res = [];
for ip = 1:numel(PASS_NAMES)
  R = analyse(PASS_NAMES{ip}, vvel_dir, mp_dir, DEPTHS, MAX_BASELINE, MIN_PAIRS);
  if isempty(R), continue; end
  % pi is the PASS_NAMES index, so palette colour/marker stay tied to the
  % product identity even when an earlier product was skipped
  R.pi = ip;
  if isempty(res), res = R; else, res(end+1) = R; end %#ok<SAGROW>
end
assert(~isempty(res), 'no products could be analysed');

jz = find(DEPTHS == ALONG_DEPTH, 1);
assert(~isempty(jz), 'ALONG_DEPTH %g is not one of DEPTHS', ALONG_DEPTH);
span_days = median([res.span]);

%% Expected curves, in mm of column thickness change
kappa   = 1/L_FLEX^2;
epsx_s  = (H_ICE/2)*kappa;
exp_tide_mm = 1e3 * (-NU/(1-NU)) * epsx_s * (1 - DEPTHS/H_ICE) .* DEPTHS;
% Firn compaction is distributed through the firn column; over the top z it
% scales with how much of that column is included, taken as linear in z.
compact_rate = ACC_ICE_EQ * (RHO_ICE/RHO_SURF - 1);            % [m/yr] firn column
exp_sec_mm   = -1e3 * compact_rate * (span_days/365.25) ...
  * min(DEPTHS/BCO_DEPTH, 1);   % saturates once the column passes close-off

%% Floor, from two builds of the same leg
[FL_SEC_MM, FL_TIDE_MM, nfl] = same_leg_floor(res, FLOOR_A, FLOOR_B, numel(DEPTHS));

%% Report
fprintf('\nObservation window: %.1f days\n', span_days);
fprintf('\n===== TIDAL: column thickness change per metre of tide [mm] =====\n');
fprintf('%-16s %10s | expected %.2f mm, floor %.2f mm (at %.0f m)\n', ...
  'product','mean', exp_tide_mm(jz), FL_TIDE_MM(jz), ALONG_DEPTH);
for i = 1:numel(res)
  v = res(i).tide_mm(:,jz); v = v(isfinite(v));
  if isempty(v), continue; end
  fprintf('%-16s %10.2f\n', res(i).name, mean(v));
end
fprintf('\n===== SECULAR: column thickness change over %.1f days [mm] =====\n', span_days);
fprintf('%-16s %10s | expected %.2f mm, floor %.2f mm (at %.0f m)\n', ...
  'product','mean', exp_sec_mm(jz), FL_SEC_MM(jz), ALONG_DEPTH);
for i = 1:numel(res)
  v = res(i).sec_mm(:,jz); v = v(isfinite(v));
  if isempty(v), continue; end
  fprintf('%-16s %10.2f\n', res(i).name, mean(v));
end
fprintf('\nFloor / expected ratio at %.0f m: tidal %.1fx, secular %.1fx\n', ...
  ALONG_DEPTH, FL_TIDE_MM(jz)/abs(exp_tide_mm(jz)), ...
  FL_SEC_MM(jz)/abs(exp_sec_mm(jz)));

%% Figure
[h, GRL] = grl_figure(140, 162.0);
set(0,'CurrentFigure',h);
axstyle = {'GridAlpha',0.15,'XColor',PAL.ink,'YColor',PAL.ink,'Box','off'};

% (a) tidal squeeze against depth
ax1 = axes('parent',h,'Position',[0.16 0.730 0.54 0.200]);
set(0,'CurrentFigure',h); hold(ax1,'on');
band_x(ax1, FL_TIDE_MM, DEPTHS, PAL.band);
plot(ax1, [0 0], [0 max(DEPTHS)], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
hexp = plot(ax1, exp_tide_mm, DEPTHS, '--', 'Color', PAL.expect, 'LineWidth', 2);
hleg = []; lbl = {}; allv = [];
for i = 1:numel(res)
  prof = nanmean_cols(res(i).tide_mm);
  allv = [allv prof]; %#ok<AGROW>
  hleg(end+1) = plot(ax1, prof, DEPTHS, '-', 'Color', PAL.cat(res(i).pi,:), 'LineWidth', 2); %#ok<AGROW>
  plot(ax1, prof, DEPTHS, PAL.cat_mk{res(i).pi}, 'MarkerSize', 7, ...
    'MarkerFaceColor', PAL.cat(res(i).pi,:), 'MarkerEdgeColor','w','LineWidth',1);
  lbl{end+1} = res(i).name; %#ok<AGROW>
end
grid(ax1,'on'); set(ax1, axstyle{:}); set(ax1,'YDir','reverse');
ylim(ax1,[0 max(DEPTHS)]); set_sym_xlim(ax1, [allv FL_TIDE_MM]);
xlabel(ax1,'Change per m of tide (mm)','Color',PAL.ink);
ylabel(ax1,'Depth below surface (m)','Color',PAL.ink);
% Titles are two short lines: the panel is 78 mm wide at 8 pt, and the
% one-line versions ran off both edges of the page.
title(ax1, {'Tidal response', 'column change per metre of tide'}, ...
  'Color',PAL.ink, 'FontWeight','normal');
lg = legend(ax1, [hleg hexp], [lbl {'flexure model'}], ...
  'Location','eastoutside','Interpreter','none');
set(lg,'TextColor',PAL.ink,'Box','off');

% (b) tidal squeeze along track
ax2 = axes('parent',h,'Position',[0.16 0.405 0.54 0.200]);
set(0,'CurrentFigure',h); hold(ax2,'on');
xmax = 0;
for i = 1:numel(res), xmax = max(xmax, max(res(i).along)/1e3); end
fill(ax2, [0 xmax xmax 0], [-1 -1 1 1]*FL_TIDE_MM(jz), PAL.band, 'EdgeColor','none');
plot(ax2, [0 xmax], [0 0], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
plot(ax2, [0 xmax], [1 1]*exp_tide_mm(jz), '--', 'Color', PAL.expect, 'LineWidth', 2);
plot(ax2, [0 xmax], [1 1]*APRES_MM, '-', 'Color', PAL.expect, 'LineWidth', 3);
for i = 1:numel(res)
  v = res(i).tide_mm(:,jz); s = res(i).tide_mm_std(:,jz); x = res(i).along/1e3;
  ok = isfinite(v);
  draw_errbars(ax2, x(ok), v(ok), s(ok), PAL.cat(res(i).pi,:));
  plot(ax2, x(ok), v(ok), PAL.cat_mk{res(i).pi}, 'MarkerSize', 7, ...
    'MarkerFaceColor', PAL.cat(res(i).pi,:), 'MarkerEdgeColor','w','LineWidth',1);
end
grid(ax2,'on'); set(ax2, axstyle{:}); xlim(ax2,[0 xmax]);
xlabel(ax2,'Along track (km)','Color',PAL.ink);
ylabel(ax2, 'Change per m of tide (mm)','Color',PAL.ink);
title(ax2, {sprintf('Along track, top %.0f m: floor %.1f mm (shaded)', ALONG_DEPTH, FL_TIDE_MM(jz)), ...
  sprintf('ApRES %+.2f mm (solid), flexure %+.1f mm (dashed)', APRES_MM, exp_tide_mm(jz))}, ...
  'Color', PAL.ink, 'FontWeight','normal');
label_band(ax2, xmax, FL_TIDE_MM(jz), PAL.ink_soft);

% (c) secular change over the observation window
ax3 = axes('parent',h,'Position',[0.16 0.080 0.54 0.200]);
set(0,'CurrentFigure',h); hold(ax3,'on');
fill(ax3, [0 xmax xmax 0], [-1 -1 1 1]*FL_SEC_MM(jz), PAL.band, 'EdgeColor','none');
plot(ax3, [0 xmax], [0 0], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
plot(ax3, [0 xmax], [1 1]*exp_sec_mm(jz), '--', 'Color', PAL.expect, 'LineWidth', 2);
for i = 1:numel(res)
  v = res(i).sec_mm(:,jz); s = res(i).sec_mm_std(:,jz); x = res(i).along/1e3;
  ok = isfinite(v);
  draw_errbars(ax3, x(ok), v(ok), s(ok), PAL.cat(res(i).pi,:));
  plot(ax3, x(ok), v(ok), PAL.cat_mk{res(i).pi}, 'MarkerSize', 7, ...
    'MarkerFaceColor', PAL.cat(res(i).pi,:), 'MarkerEdgeColor','w','LineWidth',1);
end
grid(ax3,'on'); set(ax3, axstyle{:}); xlim(ax3,[0 xmax]);
xlabel(ax3,'Along track (km)','Color',PAL.ink);
ylabel(ax3, 'Change over the window (mm)', 'Color',PAL.ink);
title(ax3, {sprintf('Steady thinning over %.1f d: floor %.0f mm (shaded)', span_days, FL_SEC_MM(jz)), ...
  sprintf('firn compaction %.1f mm (dashed)', abs(exp_sec_mm(jz)))}, ...
  'Color', PAL.ink, 'FontWeight','normal');
label_band(ax3, xmax, FL_SEC_MM(jz), PAL.ink_soft);

out_fn = fullfile(out_dir, sprintf('EAGER_2022_strain_rates%s.png', VVEL_SUFFIX));
print(h, out_fn, '-dpng', sprintf('-r%d', GRL.dpi));
close(h);
fprintf('\nWrote %s\n', out_fn);

%% ========================================================================
function band_x(ax, fl, z, col)
ok = isfinite(fl);
if ~any(ok), return; end
fill(ax, [fl(ok) fliplr(-fl(ok))], [z(ok) fliplr(z(ok))], col, 'EdgeColor','none');
end

function label_band(ax, xmax, fl, col)
if ~isfinite(fl), return; end
text(ax, 0.99*xmax, fl, 'method floor ', 'Color', col, 'FontSize', 8, ...
  'HorizontalAlignment','right','VerticalAlignment','bottom');
end

% NOTE on labelling the expected line: an inline annotation was tried and
% removed. At this scale the expected signal sits almost on the zero line -
% which is the finding - so any label long enough to be useful runs across
% several kilometres of the axis and collides with the markers. The panel
% titles carry the number instead, and the legend carries the dash style.

function set_sym_xlim(ax, v)
v = v(isfinite(v));
if isempty(v), return; end
lim = max(abs(v))*1.12;
xlim(ax, [-lim lim]);
end

function m = nanmean_cols(A)
m = nan(1, size(A,2));
for k = 1:size(A,2)
  v = A(:,k); v = v(isfinite(v));
  if ~isempty(v), m(k) = mean(v); end
end
end

%% ========================================================================
function [fl_sec, fl_tide, n] = same_leg_floor(res, nameA, nameB, nz)
% Per-block rms disagreement between two builds of the SAME leg, at every
% depth, in mm. Blocks are matched by along-track position, not by index:
% the two products need not start in the same place or hold the same count.
fl_sec = nan(1,nz); fl_tide = nan(1,nz); n = 0;
ia = find(strcmp({res.name}, nameA), 1);
ib = find(strcmp({res.name}, nameB), 1);
if isempty(ia) || isempty(ib), return; end
A = res(ia); B = res(ib);
TOL = 100;   % [m] blocks closer than this are the same piece of ice
pairs = [];
for k = 1:numel(A.along)
  [gap, m] = min(abs(B.along - A.along(k)));
  if gap <= TOL, pairs(end+1,:) = [k m]; end %#ok<AGROW>
end
if isempty(pairs), return; end
for z = 1:nz
  dt = A.tide_mm(pairs(:,1),z) - B.tide_mm(pairs(:,2),z);
  ds = A.sec_mm(pairs(:,1),z)  - B.sec_mm(pairs(:,2),z);
  dt = dt(isfinite(dt)); ds = ds(isfinite(ds));
  if ~isempty(dt), fl_tide(z) = sqrt(mean(dt.^2)); end
  if ~isempty(ds), fl_sec(z)  = sqrt(mean(ds.^2)); end
end
d = A.tide_mm(pairs(:,1),1) - B.tide_mm(pairs(:,2),1);
n = nnz(isfinite(d));
end

%% ========================================================================
function R = analyse(pass_name, vvel_dir, mp_dir, DEPTHS, MAX_BASELINE, MIN_PAIRS)
% Network version: load EVERY pair, build (block x depth) strain rows,
% invert the pass network once, then fit time+tide per depth. MIN_PAIRS
% here gates the PAIR COUNT of the network, not a per-block sample count.
R = [];
f = dir(fullfile(vvel_dir, [pass_name '_vvel_*.mat']));
keep = ~cellfun('isempty', regexp({f.name}, ...
  ['^' regexptranslate('escape',pass_name) '_vvel_\d+_\d+\.mat$'], 'once'));
f = f(keep);
if isempty(f), warning('no products for %s', pass_name); return; end

L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pass_name)), 'pass');
Np = numel(L.pass); elev = nan(1,Np); ptime = nan(1,Np);
for k = 1:Np
  elev(k)  = mean(L.pass(k).elev,'omitnan');
  ptime(k) = mean(L.pass(k).gps_time,'omitnan');
end
clear L;

nz = numel(DEPTHS);
P = []; D = []; W = []; along = []; lat = []; lon = []; Nblk = 0;
for q = 1:numel(f)
  tok = regexp(f(q).name, ['^' regexptranslate('escape',pass_name) '_vvel_(\d+)_(\d+)\.mat$'], ...
    'tokens','once');
  if isempty(tok), continue; end
  o = load(fullfile(vvel_dir, f(q).name));
  if ~vdef.pairAligned(o), continue; end
  if max(abs(o.baseline_y)) > MAX_BASELINE, continue; end
  if Nblk == 0
    Nblk = numel(o.S1);
    along = o.Along_track(:); lat = o.Latitude(:); lon = o.Longitude(:);
  end
  sv = nan(Nblk*nz, 1);          % rows ordered (b,z) = (b-1)*nz + z
  for b = 1:Nblk
    d = o.depth_blk(:,b); ok = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(ok), continue; end
    for z = 1:nz
      if max(d(ok)) < DEPTHS(z), continue; end
      sv((b-1)*nz + z) = interp1(d(ok), o.dh_blk(ok,b), DEPTHS(z),'linear',NaN) ...
        / DEPTHS(z);
    end
  end
  if all(~isfinite(sv)), continue; end
  if isempty(D), D = sv; else, D(:,end+1) = sv; end %#ok<AGROW>
  P(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
  W(end+1) = max(mean(o.coh_blk(:),'omitnan'),1e-3); %#ok<AGROW>
end
if size(P,1) < MIN_PAIRS
  warning('%s: only %d usable pairs', pass_name, size(P,1));
  return;
end

N = vdef.invertNetwork(P, D, struct('n_sigma',3,'weights',W,'n_epoch',Np));
tday = (ptime - min(ptime))/86400;
% The tide is CATS2008 at each pass mid-time, with the pass gate, from the
% same helper the flexure inversion uses (scripts/diagnostics/pass_tide.m)
% - NOT the line-mean GPS height, which carries a per-pass height error
% that attenuated this admittance by 15-25% and let a 1.1 m bad pass
% through (see vdef.surfaceAdmittance).
tide = pass_tide(pass_name, mp_dir);
span = max(tday) - min(tday);

R = struct('name', pass_name, 'along', along, 'lat', lat, 'lon', lon, ...
  'npair', size(P,1), 'span', span, ...
  'tide_mm', nan(Nblk,nz), 'tide_mm_std', nan(Nblk,nz), ...
  'sec_mm',  nan(Nblk,nz), 'sec_mm_std',  nan(Nblk,nz));
for z = 1:nz
  X = N.x((0:Nblk-1)*nz + z, :);          % Nblk x Nep at this depth
  A = vdef.fitTideAdmittance(X, tday, tide);
  R.tide_mm(:,z)     = 1e3 * A.admittance(:)     * DEPTHS(z);
  R.tide_mm_std(:,z) = 1e3 * A.admittance_std(:) * DEPTHS(z);
  R.sec_mm(:,z)      = 1e3 * A.trend(:)     * span * DEPTHS(z);
  R.sec_mm_std(:,z)  = 1e3 * A.trend_std(:) * span * DEPTHS(z);
end
end

%% ========================================================================
function draw_errbars(ax, x, y, s, col)
% Explicit segments rather than errorbar(): the gnuplot backend Octave
% falls back to on a headless box renders caps inconsistently, and this
% keeps the styling identical to the markers.
ok = isfinite(s);
for k = 1:numel(x)
  if ~ok(k), continue; end
  plot(ax, [x(k) x(k)], [y(k)-s(k) y(k)+s(k)], '-', 'Color', col, 'LineWidth', 1);
end
end
