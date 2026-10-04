function tidal_stack_compare_map(opts)
%TIDAL_STACK_COMPARE_MAP Two coherent-stack builds compared in map view.
%
%   tidal_stack_compare_map() draws the dropped-acquisition test (all
%   passes the gates allow against the flagged acquisitions dropped).
%   tidal_stack_compare_map(opts) compares any two tidal_stack outputs on
%   the same block grid:
%     opts.files     {first, second} .mat names in out_dir
%     opts.fields    {first, second} response fields to draw (default
%                    {'a','a'}; 'a_t' is the trend-controlled estimate),
%                    with opts.sd_fields their errors ({'a_sd','a_sd'})
%     opts.labels    {first, second} panel titles
%     opts.out_name  output .png name
%     opts.dclim     colour range of the change panel in sigma (default 1)
%     opts.change    false draws (a) and (b) only, for two fields that are
%                    not two estimates of one quantity (default true)
%   e.g. the standard build against the surface-coupled rebuild:
%     tidal_stack_compare_map(struct('files', {{'tidal_stack.mat', 'tidal_stack_nozc.mat'}}, ...
%       'labels', {{'Standard build (z-motion compensated)', 'Surface-coupled rebuild'}}, ...
%       'out_name', 'EAGER_2022_tidal_stack_nozc_map.png'))
%   or, within one build, the tide-only scan against the trend-controlled one:
%     tidal_stack_compare_map(struct('files', {{'tidal_stack_nozc.mat', 'tidal_stack_nozc.mat'}}, ...
%       'fields', {{'a', 'a_t'}}, 'sd_fields', {{'a_sd', 'a_t_sd'}}, ...
%       'labels', {{'Tide alone', 'Tide with the secular trend co-estimated'}}, ...
%       'out_name', 'EAGER_2022_tidal_stack_trend_map.png'))
%   or the in-phase and quadrature responses of the phase-lag fit:
%     tidal_stack_compare_map(struct('files', {{'tidal_stack_nozc.mat', 'tidal_stack_nozc.mat'}}, ...
%       'fields', {{'a_q', 'b_q'}}, 'sd_fields', {{'a_q_sd', 'b_q_sd'}}, 'change', false, ...
%       'labels', {{'In phase with the tide', 'In phase with the tide rate (quadrature)'}}, ...
%       'out_name', 'EAGER_2022_tidal_phase_map.png'))
%
%   WHY A MAP. tidal_stack_dropped_figure.m draws each line against its own
%   seaward distance, and that axis is not shared: each line's zero is its
%   own flexure-fit hinge, and the lines are 120-213 m apart and offset
%   along track, so the same x on two panels is not the same ice. Here every
%   block sits at its own position on one common frame, so lines can be
%   compared where they actually run side by side.
%
%   WHAT IS DRAWN. Three panels on identical axes:
%     (a) the first build's coherent stack at 100 m (by default, every pass
%         the gates allow),
%     (b) the second build (by default, the acquisitions pass_quality.m
%         flagged removed), with each line's pair count first -> second,
%     (c) the per-block change (b - a) in units of the combined bootstrap
%         sigma.
%   (a) and (b) use the project's diverging tidal-response scale (+/-CLIM)
%   with saturation carrying significance exactly as in
%   tidal_response_map.m: full colour only at NSIG sigma, washed toward the
%   neutral midpoint below it (ribbon_line). Each block is one segment
%   spanning its own along-track extent; nothing is smoothed or
%   interpolated between blocks.
%
%   THE FRAME. EPSG:3031 km, ROTATED so the survey runs left to right with
%   the grounding end on the left. True scale in both directions (axis
%   equal); a north arrow gives the orientation. The survey is a ~5 km by
%   ~0.9 km sliver, so north-up panels side by side would each need to be
%   over 200 mm tall. The x origin is the mean position of the four lines'
%   seaward-frame zeros projected onto the survey axis, so x reads as
%   distance along the survey on ONE axis shared by every line. There is no
%   hillshade: REMA is a north-up raster and an image cannot be rotated in
%   place without resampling it, which a comparison figure does not need.
%
%   BLOCK POSITIONS. tidal_stack.m now stores each block's lat/lon. Products
%   saved before that carry only x_sea, and for those the position is
%   recovered exactly from flexure_fit_cats.mat, whose x_flip_sign and
%   x_flip_ref are what tidal_stack.m used to make x_sea from the main pass
%   along-track coordinate, and whose full-resolution track maps that
%   coordinate to lat/lon.
%
%   Reads the two builds from vdef.figureDir (tidal_stack.mat and
%   tidal_stack_dropped.mat by default) and writes opts.out_name there
%   (EAGER_2022_tidal_stack_dropped_map.png by default).

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));             % +vdef
addpath(here);                                   % grl_figure, ribbon_line, map helpers
if ~isfield(opts,'out_dir') || isempty(opts.out_dir), opts.out_dir = vdef.figureDir(); end
if ~isfield(opts,'zref') || isempty(opts.zref), opts.zref = 100; end
if ~isfield(opts,'files') || isempty(opts.files), opts.files = {'tidal_stack.mat', 'tidal_stack_dropped.mat'}; end
if ~isfield(opts,'labels') || isempty(opts.labels)
  opts.labels = {'All passes the gates allow', 'Flagged acquisitions dropped'};
end
if ~isfield(opts,'out_name') || isempty(opts.out_name), opts.out_name = 'EAGER_2022_tidal_stack_dropped_map.png'; end
if ~isfield(opts,'fields') || isempty(opts.fields), opts.fields = {'a', 'a'}; end
if ~isfield(opts,'sd_fields') || isempty(opts.sd_fields), opts.sd_fields = {'a_sd', 'a_sd'}; end
if ~isfield(opts,'change') || isempty(opts.change), opts.change = true; end
NP = 2 + logical(opts.change);   % map panels
A = load(fullfile(opts.out_dir, opts.files{1}));  BASE = A.OUT;
B = load(fullfile(opts.out_dir, opts.files{2}));  DROP = B.OUT;
F = [];
ff = fullfile(opts.out_dir, 'flexure_fit_cats.mat');
if exist(ff, 'file'), F = load(ff); F = F.(char(fieldnames(F))); end

order = {'EAGER_2022_GL1','EAGER_2022_GL2','EAGER_2022_GL3','EAGER_2022_GL4'};
lab   = {'GL1','GL2','GL3','GL4'};
CLIM = 4;          % mm/m, same scale as tidal_response_map
NSIG = 2;          % sigma for full saturation, same rule as the map
DCLIM = 1;         % sigma, colour range of the change panel (opts.dclim)
if isfield(opts,'dclim') && ~isempty(opts.dclim), DCLIM = opts.dclim; end
LW = 4.5;          % block segment width, pt
PAL.ink = [0 0 0]; PAL.soft = [0.45 0.45 0.45]; PAL.case = [0.30 0.30 0.30];
PAL.div_neg = [0.698 0.094 0.169]; PAL.div_mid = [0.941 0.937 0.925];
PAL.div_pos = [0.165 0.471 0.839];
ndiv = 256; half = ndiv/2;
dmap = [interp1([0 1],[PAL.div_neg; PAL.div_mid], linspace(0,1,half)); ...
        interp1([0 1],[PAL.div_mid; PAL.div_pos], linspace(0,1,ndiv-half))];
ps = projcrs(3031);

%% Per line: block positions, both estimates, and the change
L = struct('lab',{},'X',{},'Y',{},'x0',{},'y0',{},'ab',{},'sb',{},'ad',{},'sd',{},'z',{},'npa',{},'npb',{});
for k = 1:numel(order)
  O = BASE(strcmp({BASE.name}, order{k})); P = DROP(strcmp({DROP.name}, order{k}));
  if isempty(O) || isempty(P), continue; end
  assert(isequal(size(O.x_sea), size(P.x_sea)) && max(abs(O.x_sea - P.x_sea)) < 1e-6, ...
    '%s: baseline and dropped runs are on different block grids', lab{k});
  [lat, lon, lat0, lon0] = block_latlon(O, F);
  [X, Y] = ps_km(ps, lon, lat); [x0, y0] = ps_km(ps, lon0, lat0);
  [~, ib] = min(abs(O.zsel - opts.zref)); [~, id] = min(abs(P.zsel - opts.zref));
  ab = 1e3*O.(opts.fields{1})(ib,:).'; sb = 1e3*O.(opts.sd_fields{1})(ib,:).';
  ad = 1e3*P.(opts.fields{2})(id,:).'; sd = 1e3*P.(opts.sd_fields{2})(id,:).';
  z = (ad - ab) ./ max(sqrt(sb.^2 + sd.^2), eps);
  [~, s] = sort(O.x_sea);           % blocks in seaward order
  L(end+1) = struct('lab', lab{k}, 'X', X(s), 'Y', Y(s), 'x0', x0, 'y0', y0, ...
    'ab', ab(s), 'sb', sb(s), 'ad', ad(s), 'sd', sd(s), 'z', z(s), ...
    'npa', O.npair, 'npb', P.npair); %#ok<AGROW>
end
assert(numel(L) >= 2, 'need at least two lines in both runs');

%% Rotated frame: survey axis left to right, grounding end on the left
allX = vertcat(L.X); allY = vertcat(L.Y);
C = cov([allX allY]); [V, D] = eig(C); [~, imax] = max(diag(D)); u = V(:, imax);
sea = [mean(arrayfun(@(l) l.X(end) - l.X(1), L)); mean(arrayfun(@(l) l.Y(end) - l.Y(1), L))];
if u.'*sea < 0, u = -u; end               % +x' points seaward
nrm = [-u(2); u(1)];                       % +y' is +90 deg from the survey axis
org = [mean([L.x0]); mean([L.y0])];
rot = @(x, y) deal((x - org(1))*u(1) + (y - org(2))*u(2), (x - org(1))*nrm(1) + (y - org(2))*nrm(2));
for i = 1:numel(L), [L(i).xr, L(i).yr] = rot(L(i).X, L(i).Y); end %#ok<AGROW>
fprintf('survey axis %.1f deg from grid east; origin at the mean seaward-frame zero\n', ...
  atan2d(u(2), u(1)));

%% ApRES sites near the survey, as orientation marks only
gis_dir = '/kucresis/scratch/hoffmana_sta/vvel/gis';
if ~exist(gis_dir, 'dir'), gis_dir = fullfile(fileparts(fileparts(here)), 'data', 'gis'); end
AP = read_apres_xy(fullfile(gis_dir, 'eastwind_2022_2023_apres_xy.txt'));
for q = 1:numel(AP), [AP(q).xr, AP(q).yr] = rot(AP(q).x/1e3, AP(q).y/1e3); end

%% Figure
% Laid out in millimetres: each map panel's height follows from its width
% and the survey's true aspect, so the canvas is exactly as tall as the
% three panels need. The x range reserves room for the line names at the
% landward end and the north arrow beyond the seaward end, so neither sits
% on a ribbon or runs into a colour bar.
pad = 0.10; lab_room = 0.30; arrow_room = 0.42;
xr_all = vertcat(L.xr); yr_all = vertcat(L.yr);
half_blk = 0.5*median(diff(L(1).xr));
xl = [min(xr_all) - half_blk - lab_room, max(xr_all) + half_blk + arrow_room];
yl = [min(yr_all) - pad, max(yr_all) + pad];
if ~isempty(AP)
  near = [AP.xr] >= xl(1) + lab_room & [AP.xr] <= xl(2) & [AP.yr] > yl(1) - 0.4 & [AP.yr] < yl(2) + 0.4;
  AP = AP(near);
  if ~isempty(AP), yl = [min([yl(1), [AP.yr] - 0.10]), max([yl(2), [AP.yr] + 0.06])]; end
end
fprintf('%d ApRES site(s) inside the survey frame\n', numel(AP));

W = 170; ML = 13; MR = 19; TT = 4.6; GAP = 2.2; MB = 10; CAP = 8.5; MT = 1.0;
aw = W - ML - MR; ah = aw*diff(yl)/diff(xl);
H = MT + NP*(TT + ah) + (NP-1)*GAP + MB + CAP;
[h, GRL] = grl_figure(W, H); set(0,'CurrentFigure',h);
nx = @(mm) mm/W; ny = @(mm) mm/H;
ybot = @(p) H - MT - p*TT - (p-1)*GAP - p*ah;          % panel bottom, mm
med = @(fn) strjoin(arrayfun(@(l) sprintf('%s %s', l.lab, fn(l)), L, 'uni', 0), ',  ');
if strcmp(opts.files{1}, opts.files{2})
  tb = sprintf('(b)  %s,  at %d m', opts.labels{2}, opts.zref);
else
  tb = sprintf('(b)  %s;  pairs: %s', opts.labels{2}, med(@(l) sprintf('%d\x2192%d', l.npa, l.npb)));
end
panels = { ...
  sprintf('(a)  %s,  at %d m', opts.labels{1}, opts.zref), tb, ...
  sprintf('(c)  Change (b - a) in combined sigma;  median |change|: %s', ...
    med(@(l) sprintf('%.2f', median(abs(l.z), 'omitnan'))))};
AX = gobjects(1, NP);
for p = 1:NP
  pos = [nx(ML) ny(ybot(p)) nx(aw) ny(ah)];
  ax = axes('parent', h, 'Position', pos); hold(ax, 'on'); AX(p) = ax;
  for i = 1:numel(L)
    [sx, sy] = block_segments(L(i).xr, L(i).yr);
    plot(ax, sx, sy, '-', 'Color', PAL.case, 'LineWidth', LW + 1.0);   % casing
    switch p
      case 1, v = L(i).ab; s = L(i).sb; cl = CLIM;
      case 2, v = L(i).ad; s = L(i).sd; cl = CLIM;
      case 3, v = L(i).z;  s = [];      cl = DCLIM;
    end
    [vv, ss] = block_values(v, s);
    ribbon_line(ax, sx, sy, vv, ss, dmap, cl, NSIG, PAL.div_mid, LW);
    text(ax, min(sx) - 0.05, sy(1), L(i).lab, 'FontSize', 7, 'Color', PAL.ink, ...
      'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
  end
  for q = 1:numel(AP)
    % ApRES stations, on every panel so the eye can carry a position
    % between them; the label is on a white chip because the stations sit
    % in the 120-210 m gaps between legs, where bare text lands on a ribbon
    plot(ax, AP(q).xr, AP(q).yr, 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', PAL.ink, ...
      'MarkerEdgeColor', 'w', 'LineWidth', 0.6);
    text(ax, AP(q).xr + 0.035, AP(q).yr, AP(q).name, 'FontSize', 6, 'Color', PAL.ink, ...
      'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
      'BackgroundColor', 'w', 'Margin', 0.5);
  end
  xlim(ax, xl); ylim(ax, yl); daspect(ax, [1 1 1]);
  set(ax, 'Position', pos, 'FontSize', 7, 'Box', 'on', 'XColor', PAL.ink, ...
    'YColor', PAL.ink, 'TickDir', 'out', 'TickLength', [0.004 0.004], 'Layer', 'top', ...
    'XTick', ceil(2*xl(1))/2:0.5:floor(xl(2)), 'YTick', -0.4:0.2:0.4);
  grid(ax, 'on'); set(ax, 'GridAlpha', 0.10);
  ylabel(ax, 'Across (km)', 'FontSize', 7);
  if p == NP
    xlabel(ax, 'Distance along the survey (km), zero at the mean grounding-end hinge', 'FontSize', 7.5);
  else
    set(ax, 'XTickLabel', []);
  end
  annotation(h, 'textbox', [nx(ML) ny(ybot(p) + ah) nx(aw) ny(TT)], 'String', panels{p}, ...
    'EdgeColor', 'none', 'FontSize', 7.5, 'Color', PAL.ink, 'Interpreter', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Margin', 0);
  colormap(ax, dmap);
  if p == 3, clim(ax, [-DCLIM DCLIM]); else, clim(ax, [-CLIM CLIM]); end
end
north_arrow(AX(1), u, nrm, org, xl, yl, PAL.ink);
cbx = nx(ML + aw + 2.5); cbw = nx(2.2);
cb = colorbar(AX(2), 'eastoutside');
set(cb, 'Position', [cbx ny(ybot(2)) cbw ny(ybot(1) + ah - ybot(2))], 'FontSize', 7, ...
  'Color', PAL.ink, 'TickDirection', 'out');
set(cb.Label, 'String', 'mm per m of tide', 'FontSize', 7);
cap = sprintf(['Each segment is one block at its own position: EPSG:3031, rotated so the survey ' ...
   'runs left to right, true scale. (a, b) reach full colour at %d sigma and wash toward ' ...
   'neutral below it.'], NSIG);
if NP == 3
  cb3 = colorbar(AX(3), 'eastoutside');
  set(cb3, 'Position', [cbx ny(ybot(3)) cbw ny(ah)], 'FontSize', 7, 'Color', PAL.ink, ...
    'TickDirection', 'out', 'Ticks', linspace(-DCLIM, DCLIM, 5));
  set(cb3.Label, 'String', 'sigma', 'FontSize', 7);
  cap = [cap ' (c) is not washed, so neutral means unchanged.'];
end
for p = 1:NP, set(AX(p), 'Position', [nx(ML) ny(ybot(p)) nx(aw) ny(ah)]); end
annotation(h, 'textbox', [nx(ML) 0 nx(aw + 6) ny(CAP)], 'String', cap, ...
  'EdgeColor', 'none', 'FontSize', 6.8, 'Color', PAL.soft, 'Margin', 0, ...
  'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

out_fn = fullfile(opts.out_dir, opts.out_name);
print(h, out_fn, '-dpng', sprintf('-r%d', GRL.dpi)); close(h);
fprintf('Wrote %s\n', out_fn);
for i = 1:numel(L)
  if NP < 3, break; end
  fprintf('%s: pairs %d -> %d, median |change| %.2f sigma, max %.2f\n', L(i).lab, ...
    L(i).npa, L(i).npb, median(abs(L(i).z), 'omitnan'), max(abs(L(i).z)));
end
end

%% ========================================================================
function [lat, lon, lat0, lon0] = block_latlon(O, F)
%BLOCK_LATLON Block centres, and the line's seaward-frame zero, in lat/lon.
%   Uses the positions tidal_stack.m stores; for a product saved before it
%   stored them, inverts x_sea back to the main-pass along-track coordinate
%   with the flexure fit's own flip, and reads the fit's full track there.
if isempty(F), error('flexure_fit_cats.mat is needed to place line %s', O.name); end
i = find(strcmp({F.lines.name}, O.name), 1);
assert(~isempty(i), '%s is not in flexure_fit_cats.mat', O.name);
T = F.lines(i);
at = T.along_full(:); ok = isfinite(at) & isfinite(T.lat_full(:)) & isfinite(T.lon_full(:));
[at, iu] = unique(at(ok)); la = T.lat_full(ok); lo = T.lon_full(ok); la = la(iu); lo = lo(iu);
lat0 = interp1(at, la, T.x_flip_ref, 'linear', 'extrap');
lon0 = interp1(at, lo, T.x_flip_ref, 'linear', 'extrap');
if isfield(O, 'lat') && ~isempty(O.lat)
  lat = O.lat(:); lon = O.lon(:); return;
end
xb = T.x_flip_ref + O.x_sea(:) / T.x_flip_sign;
assert(all(xb >= at(1) - 1 & xb <= at(end) + 1), '%s: blocks fall off the fit track', O.name);
lat = interp1(at, la, xb, 'linear'); lon = interp1(at, lo, xb, 'linear');
end

%% ========================================================================
function [sx, sy] = block_segments(x, y)
%BLOCK_SEGMENTS One segment per block, from the midpoint with its landward
%   neighbour to the midpoint with its seaward one (half a spacing beyond
%   the end blocks), NaN-separated so ribbon_line colours each by its own
%   block and never bridges two.
n = numel(x); x = x(:); y = y(:);
mx = (x(1:end-1) + x(2:end))/2; my = (y(1:end-1) + y(2:end))/2;
ex0 = [2*x(1) - mx(1); mx]; ey0 = [2*y(1) - my(1); my];
ex1 = [mx; 2*x(end) - mx(end)]; ey1 = [my; 2*y(end) - my(end)];
sx = reshape([ex0 ex1 nan(n,1)].', [], 1); sy = reshape([ey0 ey1 nan(n,1)].', [], 1);
end

%% ========================================================================
function [vv, ss] = block_values(v, s)
%BLOCK_VALUES Each block's value repeated at both ends of its segment.
n = numel(v);
vv = reshape([v(:) v(:) nan(n,1)].', [], 1);
if isempty(s), ss = []; else, ss = reshape([s(:) s(:) nan(n,1)].', [], 1); end
end

%% ========================================================================
function north_arrow(ax, u, nrm, org, xl, yl, ink)
%NORTH_ARROW True north at the survey, drawn in the rotated frame.
%   In polar stereographic, true north at any point is straight toward the
%   pole, i.e. along -(x, y) from that point. Expressed on the survey axes
%   (u along, nrm across), that direction is the arrow.
g = -org / norm(org);
d = [g.'*u; g.'*nrm]; d = d / norm(d);
len = 0.24;                                % km, fits the reserved seaward room
c = [xl(2) - 0.20; mean(yl) - 0.04];       % arrow centre
p0 = c - 0.5*len*d; p1 = c + 0.5*len*d;
th = atan2(d(2), d(1)); hl = 0.07;
plot(ax, [p0(1) p1(1)], [p0(2) p1(2)], '-', 'Color', ink, 'LineWidth', 0.9);
patch(ax, [p1(1), p1(1) - hl*cos(th - 0.4), p1(1) - hl*cos(th + 0.4)], ...
  [p1(2), p1(2) - hl*sin(th - 0.4), p1(2) - hl*sin(th + 0.4)], ink, 'EdgeColor', 'none');
text(ax, c(1), c(2) + 0.10, 'N', 'FontSize', 7, 'FontWeight', 'bold', 'Color', ink, ...
  'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
