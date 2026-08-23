function OUT = basal_crevasse_map(opts)
%BASAL_CREVASSE_MAP Basal crevasses in the tracked basal topography.
%
%   WHY. The local-E* profile softens seaward of the grounding zone, and
%   one candidate mechanism is basal crevassing: rigidity goes as
%   effective h^3, a basal crevasse incises upward into the column, and
%   the bed tracker picks the deepest coherent return so the weakened
%   section still reads as full thickness - the weakness lands in E*
%   instead. WHAT THE MAP FINDS, in the driver's own x_sea frame (origin
%   at the centre of the first valid 500 m admittance block, see below):
%   isolated keels at x_sea ~2.04 km (6.6 m, all four segments) and
%   ~2.69-2.74 km (5.9-7.6 m, three to four segments, crossing adjacent
%   corridors), plus smaller confirmed features near 2.96 and 3.50 km,
%   and dense confirmed relief at ~0.31-0.90 km on the rising bed of the
%   grounding zone, where flexure works the base but where slope-break
%   topography is an equally available reading. Against the soft-E* band
%   at x_sea 2.25-3.25 km the co-location is PARTIAL: the ~2.7 km keel
%   lies inside the band, while the deepest keel at ~2.04 km sits at the
%   band's landward edge, about 0.2 km outside it. So the crevasses are
%   consistent with driving the seaward softening without the alignment
%   being one-to-one. The visible keel height understates the fracture:
%   a 7 m open keel in a 290 m column cuts D by only ~7%, so the
%   patch-scale softening requires the weakened zone to extend well above
%   what the pick can see, which is ordinary for basal crevasses.
%
%   HOW. Full-resolution (~15 m) basal ELEVATION from the tracked bed of
%   the four main-pass segments: z_b = surface elevation - firn-corrected
%   thickness. Each segment walks the same line four times, so the
%   segment is split into continuous traverses at turnarounds, each
%   traverse is high-passed against a 1 km running median, and upward
%   anomalies wider than one sample and taller than a FIXED 3 m floor
%   become candidates (see the threshold comment in the body: a sigma
%   multiplier self-masks crevasse fields). A candidate only COUNTS when
%   it repeats in a second segment within 40 m on the map: the segments
%   are different days, so pick noise and off-nadir clutter do not repeat
%   while geometry does. Everything else is reported as unconfirmed
%   texture.
%
%   THE ALONG-TRACK FRAME. x_sea here is the SAME coordinate the driver's
%   gps_profile assigns to the elasticity products: seaward distance on
%   GL3's main-pass axis with the origin at the CENTRE OF THE FIRST VALID
%   500 m ADMITTANCE BLOCK, replicated below from the multipass ref_z
%   exactly as elastic_modulus.m computes it. Anchoring at the track
%   endpoint instead would slide everything ~0.25 km seaward of the frame
%   the soft-E* band was read in, and the co-location claim would be made
%   across two different rulers.
%
%   The soft-E* band drawn on the profile panel is x_sea 2.25-3.25 km,
%   the mid-line dip of elasticity_map.m - marked as context, and only
%   marginally resolved there, which the docstring of that script owns.
%
%   opts, all optional:
%     .out_dir   where the png goes
%     .min_amp   incision height floor [m] (default 3)
%     .hp_win    high-pass window [m] (default 1000)
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "addpath('<code>/scripts/figures'); basal_crevasse_map"

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));                       % +vdef

if ~isfield(opts,'out_dir') || isempty(opts.out_dir)
  opts.out_dir = '/kucresis/scratch/hoffmana_sta/vvel/figures_flexure';
end
if ~isfield(opts,'min_amp') || isempty(opts.min_amp), opts.min_amp = 3;    end
if ~isfield(opts,'hp_win')  || isempty(opts.hp_win),  opts.hp_win  = 1000; end

LB = '/kucresis/scratch/hoffmana_sta/vvel/opr_out/accum/2022_Antarctica_Ground/CSARP_layer_bed';
LA = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_layer';
MP = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
GIS = '/kucresis/scratch/hoffmana_sta/vvel/gis';
SEGS  = {'20221212_03','20221210_05','20221211_09','20221211_02'};
SOFT  = [2.25 3.25];                 % km, the marginal mid-line E* dip

Pfirn = vdef.firnColumn(vdef.defaultParams());
ps = projcrs(3031);

%% Reference along-track frame: GL3's main pass, in the DRIVER's x_sea
% Replicates elastic_modulus.m gps_profile block for block: 200-sample
% (500 m) blocks of the main pass, a block VALID when at least MIN_PASS
% passes carry a finite block mean and the tide regression is solvable,
% orientation from the latitudes of the first and last valid blocks, and
% the origin at the centre of the first valid block on the grounding
% (north) end. Everything drawn or printed in x_sea below shares the
% elasticity figures' ruler because it IS their ruler.
BLOCK = 200; MIN_PASS = 5;
D = load(fullfile(MP,'EAGER_2022_GL3_multipass03.mat'),'pass','param_multipass');
mi = D.param_multipass.multipass.baseline_master_idx;
rlat = D.pass(mi).lat(:); rlon = D.pass(mi).lon(:);
ralong = D.pass(mi).along_track(:);
Nx = numel(D.pass(mi).ref_z); Np = numel(D.pass);
Z = nan(Nx, Np); tide = nan(1, Np);
for k = 1:Np
  z = D.pass(k).ref_z(:);
  if numel(z) ~= Nx, continue; end
  Z(:,k)  = z;
  tide(k) = mean(z, 'omitnan');
end
clear D;
nbk = floor(Nx/BLOCK);
okb = false(nbk,1); xbv = nan(nbk,1); latv = nan(nbk,1);
for b = 1:nbk
  idx = (b-1)*BLOCK+1 : b*BLOCK;
  zb  = mean(Z(idx,:), 1, 'omitnan');
  okp = isfinite(zb) & isfinite(tide);
  if nnz(okp) < MIN_PASS, continue; end
  X = [ones(nnz(okp),1), tide(okp).'];
  if rcond(X.'*X) < 1e-12, continue; end
  okb(b)  = true;
  xbv(b)  = mean(ralong(idx), 'omitnan');
  latv(b) = mean(rlat(idx),   'omitnan');
end
clear Z;
assert(any(okb), 'no valid admittance block on the GL3 reference line');
fb = find(okb,1,'first'); lb = find(okb,1,'last');
if latv(fb) > latv(lb), sgn = +1; ref = xbv(fb);
else,                   sgn = -1; ref = xbv(lb);
end
rx_sea = sgn * (ralong - ref);
fprintf(['x_sea frame: origin at the centre of GL3''s first valid %d-sample ' ...
         'block (%.0f m from the track end), as elastic_modulus gps_profile\n'], ...
        BLOCK, min(abs(ref - ralong(1)), abs(ralong(end) - ref)));
[rxk, ryk] = ps_km(ps, rlon, rlat);

%% Per segment: basal elevation, traverse splitting, incision candidates
T = struct('seg',{},'xk',{},'yk',{},'anom',{},'zb',{},'xsea',{});
C = struct('seg',{},'xk',{},'yk',{},'xsea',{},'amp',{},'width',{});
for s = 1:numel(SEGS)
  seg = SEGS{s};
  B = load(fullfile(LB,seg,sprintf('Data_%s_001.mat',seg)));
  Qb = load(fullfile(LB,seg,sprintf('layer_%s.mat',seg)));
  nb = cellfun(@char,Qb.lyr_name,'uni',0);
  jb = find(strcmp(nb,'bottom'), 1);
  rb = []; if ~isempty(jb), rb = find(B.id == Qb.lyr_id(jb), 1); end
  if isempty(rb)
    fprintf('%s: layer "bottom" not in the tracked product (have: %s) - segment skipped\n', ...
      seg, strjoin(nb, ', '));
    continue;
  end
  tb = B.twtt(rb, :).';

  A = load(fullfile(LA,seg,sprintf('Data_%s_001.mat',seg)));
  Qa = load(fullfile(LA,seg,sprintf('layer_%s.mat',seg)));
  na = cellfun(@char,Qa.lyr_name,'uni',0);
  ja = find(strcmp(na,'surface'), 1);
  ra = []; if ~isempty(ja), ra = find(A.id == Qa.lyr_id(ja), 1); end
  if isempty(ra)
    fprintf('%s: layer "surface" not in %s (have: %s) - segment skipped\n', ...
      seg, LA, strjoin(na, ', '));
    continue;
  end
  ts = interp1(A.gps_time(:), A.twtt(ra,:).', ...
               B.gps_time(:), 'linear', 'extrap');
  ev = interp1(A.gps_time(:), A.elev(:), B.gps_time(:), 'linear', 'extrap');

  h  = vdef.depthFromTwtt(Pfirn, tb - ts);
  zb = ev - h;                                  % basal elevation
  [xk, yk] = ps_km(ps, B.lon(:), B.lat(:));

  % MASK before anything else, in three layers, because the first run of
  % this script "confirmed" seven incisions of 30-170 m and every one was
  % an artefact at the line ends:
  %   - the turnaround arcs sit up to 200 m off the line, and the
  %     projection clamps them all onto the reference track's ENDPOINTS,
  %     so four days of arc garbage stacked into repeatable-looking
  %     clusters at the two track ends;
  %   - the tracked bed fails SHALLOW by 100+ m at the turnarounds and
  %     toward the grounding end (h down to ~55 m against a 240-295 m
  %     column) - those are pick failures, not geometry;
  %   - end garbage inflated each traverse's MAD, so the detection
  %     threshold mid-line was tens of metres and nothing real could show.
  % Repetition across segments does NOT protect against any of this: the
  % arcs and the pick failures repeat too, because the geometry repeats.
  [xsall, doff, bref] = proj_ref(xk, yk, rxk, ryk, rx_sea);
  % Heading of the sled at each sample, against the line's own bearing at
  % the projection. The four leg corridors sit 120-213 m off the GL3
  % reference, so an offset cut tight enough to kill the turnaround arcs
  % also kills three of the four corridors; bearing separates them
  % cleanly instead, because a leg runs along the line in either
  % direction while an arc sweeps through every heading.
  dx = diff(xk)*1e3; dy = diff(yk)*1e3;
  hd = atan2([dy; dy(end)], [dx; dx(end)]);
  dang = abs(atan2(sin(hd - bref), cos(hd - bref)));
  along_line = min(dang, pi - dang) < deg2rad(25);
  hmed = median(h(isfinite(h)));
  hmad = 1.4826*mad_(h(isfinite(h)));
  good = isfinite(zb) & doff < 260 & along_line ...    % a leg, not an arc
       & abs(h - hmed) < 5*max(hmad, 1) ...            % pick is sane
       & xsall > min(rx_sea) + 150 ...                 % off the endpoints
       & xsall < max(rx_sea) - 150;
  zb(~good) = NaN;

  % Split into continuous traverses at the gaps the mask just opened (and
  % at any real spacing jump). With the arcs removed this yields the four
  % walks of the line, each analysed on its own.
  stp = hypot(dx, dy);
  brk = [false; stp > 60] | [false; diff(~good) ~= 0];
  id  = cumsum(brk);

  for t = unique(id).'
    sel = find(id == t & isfinite(zb) & isfinite(xk));
    if numel(sel) < 40, continue; end
    zt = zb(sel);
    win = max(3, 2*floor(opts.hp_win/mean(stp(min(sel(1:end-1),numel(stp))))/2)+1);
    base = movmedian(zt, win, 'omitnan');
    an   = zt - base;
    % A FIXED amplitude floor, deliberately with NO sigma multiplier. A
    % per-traverse robust sigma SELF-MASKS a crevasse field: mid-line the
    % traverse's own keels inflate sigma to 1.6-2.5 m, a 4-sigma bar to
    % 6-10 m, and the 5-7.6 m incisions at x_sea ~2.04 and ~2.70 km -
    % present in every segment that crosses them, repeating to <30 m -
    % fall exactly under it. Pick noise between same-corridor repeats is well under a
    % metre (the four segments agree on each keel's amplitude to ~0.5 m),
    % so the floor guards against ripple and the CROSS-SEGMENT gate below
    % is the real false-positive control.
    thr  = opts.min_amp;

    xs = proj_ref(xk(sel), yk(sel), rxk, ryk, rx_sea);
    T(end+1) = struct('seg',seg,'xk',xk(sel),'yk',yk(sel),'anom',an, ...
                      'zb',zt,'xsea',xs); %#ok<AGROW>

    % Runs of upward anomaly: base locally HIGHER, column locally thinner
    up = an > thr;
    d  = diff([0; up; 0]);
    i0 = find(d == 1); i1 = find(d == -1) - 1;
    for r = 1:numel(i0)
      w = (i1(r)-i0(r)+1) * mean(stp(min(sel(1:end-1),numel(stp))));
      if w < 20 || w > 600, continue; end       % crevasse-scale widths only
      [amp, ip] = max(an(i0(r):i1(r)));
      ip = i0(r) + ip - 1;
      C(end+1) = struct('seg',seg,'xk',xk(sel(ip)),'yk',yk(sel(ip)), ...
        'xsea',xs(ip),'amp',amp,'width',w); %#ok<AGROW>
    end
  end
  fprintf('%s: %d traverses kept, %d candidate incisions\n', seg, ...
    nnz(arrayfun(@(q) strcmp(q.seg,seg), T)), ...
    nnz(arrayfun(@(q) strcmp(q.seg,seg), C)));
end
assert(~isempty(T), 'no traverses survived - check the tracked bed product');

%% Cross-segment confirmation
% A real incision is geometry: it repeats on another day. Greedy cluster
% by map distance; confirmed = members from >= 2 distinct segments.
used = false(1, numel(C));
K = struct('xk',{},'yk',{},'xsea',{},'amp',{},'width',{},'nseg',{},'segs',{});
for i = 1:numel(C)
  if used(i), continue; end
  d = hypot([C.xk]-C(i).xk, [C.yk]-C(i).yk)*1e3;
  m = find(d < 40 & ~used);
  used(m) = true;
  segs_m = unique({C(m).seg});
  K(end+1) = struct('xk',mean([C(m).xk]),'yk',mean([C(m).yk]), ...
    'xsea',mean([C(m).xsea]),'amp',max([C(m).amp]), ...
    'width',mean([C(m).width]),'nseg',numel(segs_m), ...
    'segs',{segs_m}); %#ok<AGROW>
end
conf = [K.nseg] >= 2;
fprintf('\n%d candidate clusters, %d CONFIRMED in >=2 segments:\n', ...
  numel(K), nnz(conf));
fprintf('(x_sea in the elastic_modulus gps_profile frame, the elasticity figures'' axis)\n');
fprintf('%9s %8s %8s %6s %28s\n','x_sea km','amp m','width m','nseg','segments');
[~, order] = sort([K.xsea]);
for i = order
  if ~conf(i), continue; end
  fprintf('%9.2f %8.1f %8.0f %6d %28s\n', K(i).xsea/1e3, K(i).amp, ...
    K(i).width, K(i).nseg, strjoin(K(i).segs,' '));
end

%% Draw
PAL.ink = [0 0 0]; PAL.ink_soft = [0.45 0.45 0.45];
PAL.div_neg = [0.698 0.094 0.169]; PAL.div_mid = [0.941 0.937 0.925];
PAL.div_pos = [0.165 0.471 0.839];
ndiv = 256; half = ndiv/2;
dmap = [interp1([0 1],[PAL.div_pos; PAL.div_mid], linspace(0,1,half)); ...
        interp1([0 1],[PAL.div_mid; PAL.div_neg], linspace(0,1,ndiv-half))];
CLIMZ = 4;                                       % m of basal relief

hf = figure('Visible','off','Position',[100 100 1240 560],'Color','w');
set(0,'CurrentFigure',hf);

% (a) map of basal relief anomaly, confirmed incisions circled
axm = axes('parent',hf,'Position',[0.055 0.11 0.42 0.83]);
hold(axm,'on');
for q = 1:numel(T)
  ci = max(1, min(ndiv, round((T(q).anom+CLIMZ)/(2*CLIMZ)*(ndiv-1))+1));
  scatter(axm, T(q).xk, T(q).yk, 6, dmap(ci,:), 'filled');
end
for i = find(conf)
  plot(axm, K(i).xk, K(i).yk, 'o', 'MarkerSize', 11, ...
    'MarkerEdgeColor', PAL.ink, 'LineWidth', 1.2);
end
axis(axm,'equal');
xl = xlim(axm); ylm = ylim(axm);
xlim(axm, xl + 0.08*diff(xl)*[-1 1]); ylim(axm, ylm + 0.08*diff(ylm)*[-1 1]);
overlay_gl(axm, GIS);
grid(axm,'on');
set(axm,'GridAlpha',0.15,'XColor',PAL.ink,'YColor',PAL.ink,'Box','off');
xlabel(axm,'Polar stereographic x (km, EPSG:3031)','Color',PAL.ink);
ylabel(axm,'Polar stereographic y (km, EPSG:3031)','Color',PAL.ink);
colormap(axm, dmap); caxis(axm, [-CLIMZ CLIMZ]);
cb = colorbar(axm);
set(get(cb,'ylabel'),'string','Basal relief anomaly (m, + = up into ice)','Color',PAL.ink);
set(cb,'XColor',PAL.ink,'YColor',PAL.ink);
if rema_underlay(axm, fullfile(GIS,'rema','17_33_10m_v2.0_browse.tif'))
  grid(axm,'off');
end

% (b) basal topography along the line, incisions and the soft-E* band
axp = axes('parent',hf,'Position',[0.55 0.42 0.42 0.52]);
hold(axp,'on');
yl_all = [];
for q = 1:numel(T)
  [xs, is] = sort(T(q).xsea/1e3);
  plot(axp, xs, T(q).zb(is), '-', 'Color', [0.55 0.55 0.55 0.35], ...
    'LineWidth', 0.5);
  yl_all = [yl_all; T(q).zb(:)]; %#ok<AGROW>
end
ylp = [min(yl_all)-2 max(yl_all)+2];
patch(axp, [SOFT fliplr(SOFT)], ylp([1 1 2 2]), [0.97 0.92 0.90], ...
  'EdgeColor','none');
for q = 1:numel(T)
  [xs, is] = sort(T(q).xsea/1e3);
  plot(axp, xs, T(q).zb(is), '-', 'Color', [0.45 0.45 0.45 0.5], ...
    'LineWidth', 0.5);
end
for i = find(conf)
  plot(axp, K(i).xsea/1e3, ylp(2)-1.0, 'v', 'MarkerSize', 7, ...
    'MarkerFaceColor', PAL.div_neg, 'MarkerEdgeColor', 'none');
end
XLP = [min(rx_sea) max(rx_sea)]/1e3 + [-0.2 0.2];
ylim(axp, ylp); xlim(axp, XLP);
grid(axp,'on'); box(axp,'on');
set(axp,'XTickLabel',[],'XColor',PAL.ink,'YColor',PAL.ink);
ylabel(axp,'Basal elevation (m)','Color',PAL.ink);
title(axp,'(b)  Basal topography; triangles = confirmed incisions','Color',PAL.ink);

% (c) incision amplitude against distance from the grounding zone
axh = axes('parent',hf,'Position',[0.55 0.11 0.42 0.24]);
hold(axh,'on');
patch(axh, [SOFT fliplr(SOFT)], [0 0 1 1]*20, [0.97 0.92 0.90], ...
  'EdgeColor','none');
for i = find(conf)
  plot(axh, K(i).xsea/1e3, K(i).amp, 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', PAL.div_neg, 'MarkerEdgeColor','none');
end
for i = find(~conf)
  plot(axh, K(i).xsea/1e3, K(i).amp, 'o', 'MarkerSize', 4, ...
    'MarkerEdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.7);
end
xlim(axh, XLP);
amax = 1; if ~isempty(K), amax = max([K.amp]); end
ylim(axh, [0 amax*1.15]);
grid(axh,'on'); box(axh,'on');
set(axh,'XColor',PAL.ink,'YColor',PAL.ink);
xlabel(axh,'Seaward distance (km)','Color',PAL.ink);
ylabel(axh,'Incision height (m)','Color',PAL.ink);
title(axh,'(c)  Filled = confirmed in >=2 segments; open = single-segment','Color',PAL.ink);

if ~exist(opts.out_dir,'dir'), mkdir(opts.out_dir); end
out_fn = fullfile(opts.out_dir,'EAGER_2022_basal_crevasses.png');
print(hf, out_fn, '-dpng', '-r140');
fprintf('\nWrote %s\n', out_fn);

OUT = struct('clusters', K, 'confirmed', conf, 'traverses', numel(T));
end

%% ========================================================================
function m = mad_(v)
m = median(abs(v - median(v,'omitnan')),'omitnan');
end

%% ========================================================================
function [xs, doff, bref] = proj_ref(xk, yk, rxk, ryk, rx_sea)
%PROJ_REF x_sea of map points, by nearest-segment projection onto the
%   GL3 reference track, plus the perpendicular distance in metres and
%   the line's own bearing at the projection - the pair that separates
%   on-line samples from turnaround arcs. Same construction as the site
%   projections elsewhere in this figure set; never the bounding-box
%   diagonal.
xs = nan(size(xk)); doff = nan(size(xk)); bref = nan(size(xk));
% decimate the 1928-pt reference to ~10 m segments for speed
step = 4;
rx = rxk(1:step:end); ry = ryk(1:step:end); rv = rx_sea(1:step:end);
for n = 1:numel(xk)
  d2 = (rx - xk(n)).^2 + (ry - yk(n)).^2;
  [~, j] = min(d2);
  j0 = max(1, j-1); j1 = min(numel(rx), j+1);
  vx = rx(j1)-rx(j0); vy = ry(j1)-ry(j0);
  L2 = vx^2 + vy^2;
  if L2 <= 0
    xs(n) = rv(j); doff(n) = sqrt(d2(j))*1e3;
    continue;
  end
  t = ((xk(n)-rx(j0))*vx + (yk(n)-ry(j0))*vy) / L2;
  t = min(max(t,0),1);
  xs(n)   = rv(j0) + t*(rv(j1)-rv(j0));
  doff(n) = hypot(xk(n) - (rx(j0)+t*vx), yk(n) - (ry(j0)+t*vy))*1e3;
  bref(n) = atan2(vy, vx);
end
end

%% ========================================================================
% The three map helpers below are copied verbatim from admittance_map.m,
% following the project's existing convention of duplicating
% rema_underlay per figure script rather than sharing a path-dependent
% helper (see also tidal_deformation.m, eastwind_survey_movie.m).

%% ========================================================================
function [xk, yk] = ps_km(ps, lon, lat)
%PS_KM Project lon/lat to EPSG:3031 Antarctic Polar Stereographic, in km.
[x, y] = projfwd(ps, lat, lon);
xk = x/1e3; yk = y/1e3;
end

%% ========================================================================
function ok = rema_underlay(ax, tif)
%REMA_UNDERLAY Draw the REMA v2 hillshade under an EPSG:3031 axes in km.
%   The axes are already in the raster's own projection, so the block
%   covering the current view is read and dropped straight in - no
%   resampling, no rotation, no reprojection. It is lifted into a light
%   gray range so the data drawn on top stays dominant, and pushed to the
%   bottom of the draw order. The block is padded 8% beyond the limits so
%   a small later limit adjustment does not expose white strips. Returns
%   false (with a message) when the tile is missing or unreadable, and
%   the figure then renders exactly as it did without imagery.
%
%   The browse tif is a COG whose overview IFDs make geotiffinfo error
%   ("multiple images ... sizes are different"), so the georeference comes
%   from georasterinfo and the pixels from imread on IFD 1.
ok = false;
if ~exist(tif,'file')
  fprintf('REMA hillshade not found (%s) - no imagery underlay.\n', tif);
  return;
end
try
  R3 = georasterinfo(tif).RasterReference;
  xl = xlim(ax); yl = ylim(ax);                     % km, EPSG:3031
  xg = xl + 0.08*diff(xl)*[-1 1];
  yg = yl + 0.08*diff(yl)*[-1 1];
  px = R3.CellExtentInWorldX;
  py = R3.CellExtentInWorldY;                       % may differ from px
  c0 = max(1, floor((xg(1)*1e3 - R3.XWorldLimits(1))/px));
  c1 = min(R3.RasterSize(2), ceil((xg(2)*1e3 - R3.XWorldLimits(1))/px));
  r0 = max(1, floor((R3.YWorldLimits(2) - yg(2)*1e3)/py));
  r1 = min(R3.RasterSize(1), ceil((R3.YWorldLimits(2) - yg(1)*1e3)/py));
  if c1 <= c0 || r1 <= r0
    fprintf('REMA tile does not cover this view - no imagery underlay.\n');
    return;
  end
  A = double(imread(tif, 'Index', 1, 'PixelRegion', {[r0 r1],[c0 c1]}));
  A(A == 0) = NaN;                                  % nodata
  vv = sort(A(isfinite(A)));
  if isempty(vv)
    fprintf('REMA tile is empty over this view - no imagery underlay.\n');
    return;
  end
  lo = vv(max(1,round(0.02*numel(vv)))); hi = vv(round(0.98*numel(vv)));
  g = min(max((A - lo)/max(hi-lo, 1), 0), 1);
  g = 0.58 + 0.40*g;                                % recessive light grays
  g(~isfinite(A)) = 1;                              % nodata as paper white
  % cell centres of the block actually read, in km; y descends with row
  ximg = (R3.XWorldLimits(1) + ([c0 c1] - 0.5)*px)/1e3;
  yimg = (R3.YWorldLimits(2) - ([r0 r1] - 0.5)*py)/1e3;
  hImg = image(ax, 'XData', ximg, 'YData', yimg, 'CData', repmat(g,[1 1 3]));
  uistack(hImg, 'bottom');
  xlim(ax, xl); ylim(ax, yl);
  ok = true;
  fprintf('REMA hillshade underlay: %d x %d px, native EPSG:3031\n', ...
    r1-r0+1, c1-c0+1);
catch ME
  fprintf('REMA underlay failed (%s) - continuing without imagery.\n', ME.message);
end
end

%% ========================================================================
function drawn = overlay_gl(ax, gis_dir)
drawn = false;
fn = fullfile(gis_dir,'GroundingLine_Antarctica_v02.shp');
if ~exist(fn,'file'), return; end
try
  S = shaperead(fn);
  X = []; Y = [];
  for k = 1:numel(S)
    X = [X; S(k).X(:); NaN]; Y = [Y; S(k).Y(:); NaN]; %#ok<AGROW>
  end
  gx = X/1e3; gy = Y/1e3;      % the shapefile is already EPSG:3031 metres
  xl = xlim(ax); yl = ylim(ax); pad = 3;
  bad = gx < xl(1)-pad | gx > xl(2)+pad | gy < yl(1)-pad | gy > yl(2)+pad;
  gx(bad) = NaN; gy(bad) = NaN;
  if all(isnan(gx)), return; end
  plot(ax, gx, gy, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2);
  xlim(ax, xl); ylim(ax, yl);
  drawn = true;
catch
end

end
