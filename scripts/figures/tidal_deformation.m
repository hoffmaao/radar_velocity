%TIDAL_DEFORMATION Vertical column strain against the tide at Windless Bight.
%   Builds the analysis that motivates the whole project, for every EAGER
%   2022 repeat-pass line: Windless Bight is floating, so if tidal flexure
%   strains the ice column then the vertical strain measured between repeat
%   passes should track the tide.
%
%   Every product in the 'main' pairing measures the SAME reference pass
%   against a different pass, so together they are a time series of the
%   column's state relative to that one epoch - which is exactly what is
%   wanted here, and is why this reads the main pairing rather than the
%   sequential one.
%
%   Two observables, both relative to the main pass:
%     tide   mean platform elevation of each pass (the passes sit on a
%            floating shelf, so the antenna rides the tide directly)
%     strain dh at REF_DEPTHS divided by that depth - the mean vertical
%            strain of the column between the surface and that depth
%
%   They are NOT plotted on a shared pair of axes: metres of heave and
%   dimensionless strain are different scales, and a twin y-axis invites
%   reading a correlation off the plot that the axis scaling created. They
%   get stacked panels over a common time axis, and the correlation gets
%   its own panels where it can be judged honestly.
%
%   THE ROBUSTNESS TESTS.
%   1. The tide series is IDENTICAL for every block of a line, so a
%      common-mode artefact would correlate with the same sign everywhere.
%      A sign REVERSAL along track has to come from the strain data.
%   2. The far-end blocks are also the ones whose coherence thins with
%      depth, so the correlation is computed over two column depths: 100 m,
%      where every block has near-full coverage, and 200 m, where the far
%      end is patchy. A reversal present at both is a property of the ice.
%   3. Pairs whose cross-track baseline exceeds MAX_BASELINE are dropped:
%      surface referencing does not remove topographic phase, and
%      EAGER_2022_GL2 in particular has passes tens of metres off track.
%   4. Five independent lines (EAGER_2022 and GL1-GL4). A hinge that is real
%      should appear in all of them, at consistent GEOGRAPHIC positions rather
%      than at the same distance along each line.
%
%   Outputs one 4-panel figure per line plus a cross-line summary figure.
%
%   Run on the server, where the products live:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../tidal_deformation.m')"

% VVEL_SUFFIX selects which processing run to read: '' for the standing
% 500 m blocks, '_fine' for the shorter along-track window. Set it before
% running to switch.
if ~exist('VVEL_SUFFIX','var'), VVEL_SUFFIX = ''; end
vvel_dir = ['/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel' VVEL_SUFFIX];
mp_dir   = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
out_dir  = '/kucresis/scratch/hoffmana_sta/vvel/figures';
gis_dir  = '/kucresis/scratch/hoffmana_sta/vvel/gis';

PASS_NAMES = {'EAGER_2022','EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};
REF_DEPTHS    = [100 200];  % [m] column depths over which strain is measured
MAIN_DEPTH    = 2;          % which of those drives the per-line panels (b), (c)
SUMMARY_DEPTH = 1;          % which drives the cross-line summary (the robust one)
MIN_OBS       = 10;         % below this a block's correlation is drawn hollow
MAX_BASELINE  = 10;         % [m] cross-track baseline above which a pair is dropped

if ~exist(out_dir,'dir'), mkdir(out_dir); end

% Palette. Blocks along a line are ORDERED -> single-hue sequential ramp.
% Lines are CATEGORIES -> validated categorical hues (5 slots, worst
% adjacent CVD dE 9.1; the sub-3:1 contrast warning on three of them is
% relieved by the legend and the printed table). The correlation map is
% POLARITY -> diverging, two poles with a neutral gray midpoint and no hue
% at the middle.
PAL = struct();
PAL.seq_anchors = [0.776 0.859 0.937; 0.419 0.682 0.839; 0.129 0.443 0.710; ...
                   0.032 0.271 0.580; 0.031 0.188 0.420];
PAL.cat = [0.165 0.471 0.839;    % #2a78d6 blue
           0.922 0.408 0.204;    % #eb6834 orange
           0.106 0.686 0.478;    % #1baf7a aqua
           0.929 0.631 0.000;    % #eda100 yellow
           0.910 0.482 0.643];   % #e87ba4 magenta
PAL.cat_mk  = {'o','s','^','d','v'};
PAL.div_neg = [0.698 0.094 0.169];   % red pole:  strain falls as the tide rises
PAL.div_mid = [0.941 0.937 0.925];   % neutral gray midpoint
PAL.div_pos = [0.165 0.471 0.839];   % blue pole: strain rises with the tide
PAL.ink      = [0.20 0.20 0.20];
PAL.ink_soft = [0.45 0.45 0.45];

%% Analyse every line
res = [];
for ip = 1:numel(PASS_NAMES)
  r1 = analyse_line(PASS_NAMES{ip}, vvel_dir, mp_dir, REF_DEPTHS, MAX_BASELINE);
  if isempty(r1), continue; end
  if isempty(res), res = r1; else, res(end+1) = r1; end %#ok<SAGROW>
end
assert(~isempty(res),'no lines could be analysed');

%% Per-line report and figure
for i = 1:numel(res)
  R = res(i);
  fprintf('\n===== %s =====\n', R.pass_name);
  fprintf('main pass %d (%s), %d of %d pairs used (%d dropped on cross-track baseline)\n', ...
    R.main_pass, R.main_seg, R.np_used, R.np_all, R.np_all - R.np_used);
  fprintf('tide range %.3f m\n', R.tide_range);
  fprintf('%-6s %10s %9s %9s', 'block','along[km]','lat','lon');
  for j = 1:numel(REF_DEPTHS)
    fprintf(' %7s%3.0f %7s%3.0f', 'n@', REF_DEPTHS(j), 'r@', REF_DEPTHS(j));
  end
  fprintf('\n');
  for b = 1:R.Nblk
    fprintf('%-6d %10.2f %9.4f %9.4f', b, R.along(b)/1e3, R.lat(b), R.lon(b));
    for j = 1:numel(REF_DEPTHS)
      fprintf(' %10d %10.2f', R.nobs(b,j), R.r(b,j));
    end
    fprintf('\n');
  end
  res(i).suffix = VVEL_SUFFIX; R.suffix = VVEL_SUFFIX;
  plot_line_figure(R, REF_DEPTHS, MAIN_DEPTH, MIN_OBS, PAL, out_dir);
end

%% Cross-line summary
fprintf('\n===== hinge position per line (change point in r at %.0f m) =====\n', ...
  REF_DEPTHS(SUMMARY_DEPTH));
fprintf('%-18s %10s %9s %9s %9s %8s %8s\n', ...
  'line','along[km]','lat','lon','contrast','r_near','r_far');
for i = 1:numel(res)
  R = res(i);
  [xc, latc, lonc, contrast] = zero_crossing(R, SUMMARY_DEPTH);
  near = R.r(1:min(3,R.Nblk), SUMMARY_DEPTH);
  far  = R.r(max(1,R.Nblk-3):R.Nblk, SUMMARY_DEPTH);
  fprintf('%-18s %10s %9s %9s %9s %8.2f %8.2f\n', R.pass_name, ...
    numstr(xc/1e3,'%.2f'), numstr(latc,'%.4f'), numstr(lonc,'%.4f'), ...
    numstr(contrast,'%.2f'), mean(near,'omitnan'), mean(far,'omitnan'));
  res(i).xc = xc; res(i).latc = latc; res(i).lonc = lonc;
end

plot_summary_figure(res, REF_DEPTHS, SUMMARY_DEPTH, MIN_OBS, PAL, out_dir, gis_dir, VVEL_SUFFIX);

%% ========================================================================
function R = analyse_line(pass_name, vvel_dir, mp_dir, REF_DEPTHS, MAX_BASELINE)
R = [];
f = dir(fullfile(vvel_dir, [pass_name '_vvel_*.mat']));
% dir() prefix-matches, so the EAGER_2022_GL* files would also answer to
% EAGER_2022; keep only the exact line
keep = false(1,numel(f));
for i = 1:numel(f)
  keep(i) = ~isempty(regexp(f(i).name, ...
    ['^' regexptranslate('escape',pass_name) '_vvel_\d+_\d+\.mat$'], 'once'));
end
f = f(keep);
if isempty(f)
  warning('no products for %s', pass_name);
  return;
end

mp_fn = fullfile(mp_dir, sprintf('%s_multipass03.mat', pass_name));
L = load(mp_fn,'pass');
Npass = numel(L.pass);
pass_elev = nan(1,Npass); pass_seg = cell(1,Npass);
for k = 1:Npass
  pass_elev(k) = mean(L.pass(k).elev,'omitnan');
  pass_seg{k}  = L.pass(k).param_pass.day_seg;
end
clear L;

sec_per_year = 365.25*86400;
np = numel(f); nd = numel(REF_DEPTHS);
sec_idx = nan(1,np); ref_idx = nan(1,np); t_sec = nan(1,np); maxbl = nan(1,np);
strain = []; along = []; lat = []; lon = [];
for i = 1:np
  o = load(fullfile(vvel_dir, f(i).name));
  sec_idx(i) = o.pass_idx_sec;
  ref_idx(i) = o.pass_idx_ref;
  maxbl(i)   = max(abs(o.baseline_y));
  if isempty(strain)
    Nblk = numel(o.S1);
    strain = nan(Nblk, np, nd);
    along = o.Along_track(:); lat = o.Latitude(:); lon = o.Longitude(:);
  end
  t_sec(i) = mean(o.GPS_time + o.delta_t_blk*sec_per_year,'omitnan');
  for b = 1:Nblk
    d  = o.depth_blk(:,b);
    ok = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(ok), continue; end
    for j = 1:nd
      if max(d(ok)) < REF_DEPTHS(j), continue; end
      strain(b,i,j) = interp1(d(ok), o.dh_blk(ok,b), REF_DEPTHS(j), 'linear', NaN) ...
        / REF_DEPTHS(j);
    end
  end
end

% Drop pairs that carry topographic phase from a long cross-track baseline
use = ~(isfinite(maxbl) & maxbl > MAX_BASELINE);
np_all = np;
sec_idx = sec_idx(use); ref_idx = ref_idx(use); t_sec = t_sec(use);
strain = strain(:,use,:);
if numel(t_sec) < 5
  warning('%s: only %d usable pairs after the baseline cut', pass_name, numel(t_sec));
  R = []; return;
end

[t_sec, ord] = sort(t_sec);
sec_idx = sec_idx(ord);
strain  = strain(:,ord,:);
main_pass = ref_idx(1);
tide   = pass_elev(sec_idx) - pass_elev(main_pass);
t_days = (t_sec - min(t_sec))/86400;

Nblk = size(strain,1);
r = nan(Nblk,nd); slope = nan(Nblk,nd); nobs = zeros(Nblk,nd);
for j = 1:nd
  for b = 1:Nblk
    s  = strain(b,:,j);
    ok = isfinite(s) & isfinite(tide);
    nobs(b,j) = nnz(ok);
    if nobs(b,j) < 5, continue; end
    rm = corrcoef(tide(ok), s(ok));
    r(b,j) = rm(1,2);
    p = polyfit(tide(ok), s(ok), 1);
    slope(b,j) = p(1);
  end
end

R = struct('pass_name',pass_name,'main_pass',main_pass,'main_seg',{pass_seg{main_pass}}, ...
  'np_all',np_all,'np_used',numel(t_sec),'Nblk',Nblk,'along',along,'lat',lat,'lon',lon, ...
  'tide',tide,'t_days',t_days,'t_sec',t_sec,'strain',strain,'r',r,'slope',slope, ...
  'nobs',nobs,'tide_range',max(tide)-min(tide),'xc',NaN,'latc',NaN,'lonc',NaN, ...
  'suffix','');
end

%% ========================================================================
function [xc, latc, lonc, contrast] = zero_crossing(R, j)
% Along-track position of the strongest positive-to-negative transition in
% the tidal correlation - the hinge.
%
% NOT the first zero crossing. With the short along-track window there are
% 39 blocks per line and the per-block correlation is noisy, so the first
% crossing latches onto whichever wiggle happens to dip below zero near the
% start. This is a change-point instead: the split that maximises
% mean(r) before minus mean(r) after, which is driven by the whole profile
% rather than by one block, and which reports nothing when there is no real
% step (MIN_CONTRAST) instead of inventing a hinge.
MIN_CONTRAST = 0.5;   % the drop in mean r either side must be at least this
MIN_SIDE     = 3;     % blocks required on each side

xc = NaN; latc = NaN; lonc = NaN; contrast = NaN;
rr = R.r(:,j);
ok = find(isfinite(rr));
n  = numel(ok);
if n < 2*MIN_SIDE, return; end
rv = rr(ok); xx = R.along(ok);

best = -inf; kbest = [];
for k = MIN_SIDE:(n-MIN_SIDE)
  d = mean(rv(1:k)) - mean(rv(k+1:end));
  if d > best, best = d; kbest = k; end
end
if isempty(kbest) || best < MIN_CONTRAST, return; end

contrast = best;
xc   = 0.5*(xx(kbest)        + xx(kbest+1));
latc = 0.5*(R.lat(ok(kbest)) + R.lat(ok(kbest+1)));
lonc = 0.5*(R.lon(ok(kbest)) + R.lon(ok(kbest+1)));
end

%% ========================================================================
function s = numstr(v, fmt)
if isfinite(v), s = sprintf(fmt, v); else, s = 'none'; end
end

%% ========================================================================
function plot_line_figure(R, REF_DEPTHS, MAIN_DEPTH, MIN_OBS, PAL, out_dir)
cmap = interp1(linspace(0,1,size(PAL.seq_anchors,1)), PAL.seq_anchors, ...
  linspace(0,1,max(R.Nblk,2)));
ink = PAL.ink; ink_soft = PAL.ink_soft;
depth_col = {PAL.cat(1,:), PAL.cat(2,:)};
depth_mk  = {'o','s'};

h = figure('Visible','off','Position',[100 100 950 1180],'Color','w');
axstyle = {'GridAlpha',0.15,'XColor',ink_soft,'YColor',ink_soft,'Box','off'};

ax1 = axes('parent',h,'Position',[0.10 0.775 0.73 0.175]);
plot(ax1, R.t_days, R.tide, '-', 'Color', ink_soft, 'LineWidth', 1.5); hold(ax1,'on');
plot(ax1, R.t_days, R.tide, 'o', 'MarkerSize', 7, 'MarkerFaceColor', ink, ...
  'MarkerEdgeColor','w','LineWidth',1);
grid(ax1,'on'); set(ax1, axstyle{:});
ylabel(ax1,'Platform elevation (m)','Color',ink);
title(ax1, sprintf('%s  -  tide from pass GPS elevation, relative to the main pass', ...
  R.pass_name), 'Interpreter','none','Color',ink);

ax2 = axes('parent',h,'Position',[0.10 0.535 0.73 0.175]);
hold(ax2,'on');
for b = 1:R.Nblk
  if R.nobs(b,MAIN_DEPTH) < 5, continue; end
  plot(ax2, R.t_days, 1e6*R.strain(b,:,MAIN_DEPTH), '-o', 'Color', cmap(b,:), ...
    'MarkerSize', 5, 'MarkerFaceColor', cmap(b,:), 'MarkerEdgeColor','w','LineWidth',1.5);
end
grid(ax2,'on'); set(ax2, axstyle{:});
ylabel(ax2, sprintf('Strain 0-%.0f m (\\mu\\epsilon)', REF_DEPTHS(MAIN_DEPTH)),'Color',ink);
xlabel(ax2, sprintf('Days from %s UTC', ...
  datestr(epoch_to_datenum(min(R.t_sec)),'yyyy-mm-dd HH:MM')), 'Color', ink);
title(ax2,'Vertical strain of the column, relative to the main pass','Color',ink);
cb = colorbar(ax2,'Position',[0.855 0.535 0.020 0.175]);
colormap(ax2, cmap); caxis(ax2, [R.along(1) R.along(end)]/1e3);
set(get(cb,'ylabel'),'string','Along track (km)','Color',ink);
set(cb,'XColor',ink_soft,'YColor',ink_soft);

ax3 = axes('parent',h,'Position',[0.10 0.295 0.73 0.175]);
hold(ax3,'on');
for b = 1:R.Nblk
  if R.nobs(b,MAIN_DEPTH) < 5, continue; end
  s  = R.strain(b,:,MAIN_DEPTH);
  ok = isfinite(s) & isfinite(R.tide);
  plot(ax3, R.tide, 1e6*s, 'o', 'MarkerSize', 6, 'MarkerFaceColor', cmap(b,:), ...
    'MarkerEdgeColor','w','LineWidth',0.75);
  xf = linspace(min(R.tide(ok)), max(R.tide(ok)), 2);
  b0 = mean(s(ok)) - R.slope(b,MAIN_DEPTH)*mean(R.tide(ok));
  plot(ax3, xf, 1e6*(R.slope(b,MAIN_DEPTH)*xf + b0), '-', 'Color', cmap(b,:), 'LineWidth',1.5);
end
grid(ax3,'on'); set(ax3, axstyle{:});
xlabel(ax3,'Platform elevation relative to main pass (m)','Color',ink);
ylabel(ax3, sprintf('Strain 0-%.0f m (\\mu\\epsilon)', REF_DEPTHS(MAIN_DEPTH)),'Color',ink);
title(ax3,'Strain against tide, coloured by position along the line','Color',ink);

ax4 = axes('parent',h,'Position',[0.10 0.055 0.73 0.175]);
hold(ax4,'on');
plot(ax4, [R.along(1) R.along(end)]/1e3, [0 0], '-', 'Color', [0.75 0.75 0.75], 'LineWidth',1);
hleg = []; lbl = {};
for j = 1:numel(REF_DEPTHS)
  solid = R.nobs(:,j) >= MIN_OBS;
  hleg(end+1) = plot(ax4, R.along/1e3, R.r(:,j), '-', 'Color', depth_col{j}, 'LineWidth',2); %#ok<AGROW>
  plot(ax4, R.along(solid)/1e3, R.r(solid,j), depth_mk{j}, 'MarkerSize',8, ...
    'MarkerFaceColor', depth_col{j}, 'MarkerEdgeColor','w','LineWidth',1);
  plot(ax4, R.along(~solid)/1e3, R.r(~solid,j), depth_mk{j}, 'MarkerSize',8, ...
    'MarkerFaceColor','w','MarkerEdgeColor', depth_col{j}, 'LineWidth',1.5);
  lbl{end+1} = sprintf('column 0-%.0f m', REF_DEPTHS(j)); %#ok<AGROW>
end
grid(ax4,'on'); set(ax4, axstyle{:}); ylim(ax4,[-1 1]);
xlabel(ax4,'Along track (km)','Color',ink);
ylabel(ax4,'Correlation of strain with tide','Color',ink);
title(ax4, sprintf('Tidal response along the line (hollow: fewer than %d pairs)', MIN_OBS), ...
  'Color', ink);
lg = legend(ax4, hleg, lbl, 'Location','SouthWest');
set(lg,'TextColor',ink,'Box','off');

out_fn = fullfile(out_dir, sprintf('%s_tidal_deformation%s.png', R.pass_name, R.suffix));
print(h, out_fn, '-dpng', '-r120');
close(h);
fprintf('Wrote %s\n', out_fn);
end

%% ========================================================================
function plot_summary_figure(res, REF_DEPTHS, j, MIN_OBS, PAL, out_dir, gis_dir, suffix)
ink = PAL.ink; ink_soft = PAL.ink_soft;
axstyle = {'GridAlpha',0.15,'XColor',ink_soft,'YColor',ink_soft,'Box','off'};

h = figure('Visible','off','Position',[100 100 980 900],'Color','w');

% (a) Correlation against along-track distance, one series per line
ax1 = axes('parent',h,'Position',[0.09 0.60 0.60 0.33]);
hold(ax1,'on');
xmax = 0;
for i = 1:numel(res), xmax = max(xmax, res(i).along(end)/1e3); end
plot(ax1, [0 xmax], [0 0], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
hleg = []; lbl = {};
for i = 1:numel(res)
  R = res(i);
  col = PAL.cat(mod(i-1,size(PAL.cat,1))+1,:);
  mk  = PAL.cat_mk{mod(i-1,numel(PAL.cat_mk))+1};
  hleg(end+1) = plot(ax1, R.along/1e3, R.r(:,j), '-', 'Color', col, 'LineWidth', 2); %#ok<AGROW>
  solid = R.nobs(:,j) >= MIN_OBS;
  plot(ax1, R.along(solid)/1e3, R.r(solid,j), mk, 'MarkerSize',7, ...
    'MarkerFaceColor',col,'MarkerEdgeColor','w','LineWidth',1);
  plot(ax1, R.along(~solid)/1e3, R.r(~solid,j), mk, 'MarkerSize',7, ...
    'MarkerFaceColor','w','MarkerEdgeColor',col,'LineWidth',1.5);
  lbl{end+1} = R.pass_name; %#ok<AGROW>
end
grid(ax1,'on'); set(ax1, axstyle{:}); ylim(ax1,[-1 1]); xlim(ax1,[0 xmax]);
xlabel(ax1,'Along track (km)','Color',ink);
ylabel(ax1,'Correlation of strain with tide','Color',ink);
title(ax1, sprintf('Tidal response of the 0-%.0f m column, all five lines', REF_DEPTHS(j)), ...
  'Color', ink);
lg = legend(ax1, hleg, lbl, 'Location','eastoutside','Interpreter','none');
set(lg,'TextColor',ink,'Box','off');

% (b) The same correlation in map view - diverging, gray where no response
ax2 = axes('parent',h,'Position',[0.09 0.07 0.60 0.42]);
hold(ax2,'on');
% Local tangent-plane km rather than degrees: at 77.7 S a degree of
% longitude is ~4.7x shorter than a degree of latitude, so a lon/lat axis
% would misrepresent both the line spacing and the hinge geometry.
lat0 = mean(cellfun(@(v) mean(v,'omitnan'), {res.lat}));
lon0 = mean(cellfun(@(v) mean(v,'omitnan'), {res.lon}));
xy = @(lon,lat) deal((lon-lon0)*111320*cosd(lat0)/1e3, (lat-lat0)*110540/1e3);

ndiv = 256; half = round(ndiv/2);
dmap = [interp1([0 1],[PAL.div_neg; PAL.div_mid], linspace(0,1,half)); ...
        interp1([0 1],[PAL.div_mid; PAL.div_pos], linspace(0,1,ndiv-half))];
for i = 1:numel(res)
  R = res(i);
  mk = PAL.cat_mk{mod(i-1,numel(PAL.cat_mk))+1};
  ok = isfinite(R.r(:,j));
  [xk, yk] = xy(R.lon, R.lat);
  plot(ax2, xk(ok), yk(ok), '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
  for b = 1:R.Nblk
    if ~isfinite(R.r(b,j)), continue; end
    ci = max(1, min(ndiv, round((R.r(b,j)+1)/2*(ndiv-1))+1));
    plot(ax2, xk(b), yk(b), mk, 'MarkerSize', 11, ...
      'MarkerFaceColor', dmap(ci,:), 'MarkerEdgeColor', ink_soft, 'LineWidth', 0.5);
  end
  if isfinite(R.lonc)
    [xc2, yc2] = xy(R.lonc, R.latc);
    plot(ax2, xc2, yc2, 'p', 'MarkerSize', 17, 'MarkerFaceColor', 'w', ...
      'MarkerEdgeColor', ink, 'LineWidth', 1.5);
  end
end

% Grounding zone last: it needs the axis limits the survey markers set, so
% that it can be clipped to this window instead of dragging the limits out
% to the whole continent.
gz_drawn = overlay_grounding(ax2, gis_dir, xy);
grid(ax2,'on'); set(ax2, axstyle{:}); axis(ax2,'equal');
xlabel(ax2, sprintf('East of %.4f deg (km)', lon0),'Color',ink);
ylabel(ax2, sprintf('North of %.4f deg (km)', lat0),'Color',ink);
if gz_drawn
  title(ax2,'Map view, equal scale - stars mark the sign change, black line the grounding zone','Color',ink);
else
  title(ax2,'Map view, equal scale - stars mark where the response changes sign','Color',ink);
end
colormap(ax2, dmap); caxis(ax2,[-1 1]);
cb = colorbar(ax2,'Position',[0.715 0.07 0.020 0.42]);
set(get(cb,'ylabel'),'string','Correlation of strain with tide','Color',ink);
set(cb,'XColor',ink_soft,'YColor',ink_soft);

out_fn = fullfile(out_dir, sprintf('EAGER_2022_tidal_summary%s.png', suffix));
print(h, out_fn, '-dpng', '-r120');
close(h);
fprintf('\nWrote %s\n', out_fn);
end

%% ========================================================================
function drawn = overlay_grounding(ax, gis_dir, xy)
% Draw whatever grounding-line geometry is in gis_dir, converted into the
% map's local km frame. Returns false (with a note) when nothing is there,
% so the figure still builds before the data has been fetched.
%
% Accepts either a .mat holding lat/lon vectors (NaN-separated parts) or a
% shapefile. MEaSUREs Antarctic Boundaries ships EPSG:3031 metres, so
% anything that does not look like degrees is inverse-projected from the
% Antarctic Polar Stereographic grid.
drawn = false;
if ~exist(gis_dir,'dir')
  fprintf('No grounding-zone data in %s - skipping that overlay.\n', gis_dir);
  return;
end
f = [dir(fullfile(gis_dir,'*.shp')); dir(fullfile(gis_dir,'*.mat'))];
if isempty(f)
  fprintf('No .shp or .mat in %s - skipping the grounding-zone overlay.\n', gis_dir);
  return;
end

fn = fullfile(gis_dir, f(1).name);
[~,~,ext] = fileparts(fn);
lat = []; lon = [];
try
  if strcmpi(ext,'.mat')
    G = load(fn);
    flds = fieldnames(G);
    latf = flds(~cellfun('isempty', regexpi(flds,'^lat')));
    lonf = flds(~cellfun('isempty', regexpi(flds,'^lon')));
    if isempty(latf) || isempty(lonf)
      fprintf('%s has no lat/lon fields - skipping the grounding-zone overlay.\n', fn);
      return;
    end
    lat = G.(latf{1})(:); lon = G.(lonf{1})(:);
  else
    S = shaperead(fn);
    X = []; Y = [];
    for k = 1:numel(S)
      X = [X; S(k).X(:); NaN]; %#ok<AGROW>
      Y = [Y; S(k).Y(:); NaN]; %#ok<AGROW>
    end
    if max(abs(X(isfinite(X)))) <= 180 && max(abs(Y(isfinite(Y)))) <= 90
      lon = X; lat = Y;                       % already degrees
    else
      [lat, lon] = projinv(projcrs(3031), X, Y);   % EPSG:3031 metres
    end
  end
catch ME
  fprintf('Could not read %s (%s) - skipping the grounding-zone overlay.\n', fn, ME.message);
  return;
end

[gx, gy] = xy(lon, lat);
% Keep only what falls near the survey, else the whole continent sets the
% axis limits
xl = xlim(ax); yl = ylim(ax);
pad = 5;
keep = ~(gx < xl(1)-pad | gx > xl(2)+pad | gy < yl(1)-pad | gy > yl(2)+pad);
gx(~keep) = NaN; gy(~keep) = NaN;
if all(isnan(gx))
  fprintf('The grounding-zone geometry in %s does not reach this survey.\n', fn);
  return;
end
plot(ax, gx, gy, '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 2);
xlim(ax, xl); ylim(ax, yl);
fprintf('Overlaid grounding zone from %s\n', fn);
drawn = true;
end
