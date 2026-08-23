%EASTWIND_SURVEY_MOVIE Animate the EAGER 2022 surveys against the CATS2008 tide.
%   Two fixed panels, one above the other, both created once and only
%   updated per frame so their positions never shift:
%
%     top     mean CATS2008 tide over the survey array across the whole
%             window, with a marker riding the curve at the current time and
%             the timestamp in the top right corner
%     bottom  the radar system position, moving as time advances, over the
%             three days of surveying
%
%   NO PANEL TITLES OR IN-PLOT LABELS: the axis labels carry it, and any
%   naming belongs in the slide or caption, same convention as the other
%   figures in this repo. The only text in a frame is the axis labels,
%   the tick labels and the running timestamp - the timestamp stays
%   because it is the movie's clock, not a caption (the dated x axis and
%   the riding marker give the same information more coarsely).
%
%   The tide is the MEAN over the array: CATS2008 is predicted at a set of
%   points spread across the survey footprint and averaged, rather than at a
%   single centroid. At this site the model is ~4 km resolution against a
%   ~5 km array, so the spread is small, but averaging is what was asked for
%   and it costs nothing.
%
%   Uses the MATLAB Tide Model Driver already on the server rather than
%   pyTMD, which is not installed. Same CATS2008 model either way.
%
%   Run on the server, where the products, the TMD driver, and the CATS2008
%   model live:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../eastwind_survey_movie.m')"
%
%   The output is Motion JPEG AVI rather than MPEG-4 because the Linux
%   MATLAB on the server has no MPEG-4 VideoWriter profile and no ffmpeg is
%   available; VLC and QuickTime both play it. Copy the result into figs/
%   locally (gitignored - a 200+ MB artefact, not a repo file).

TMD_DIR   = '/kucresis/scratch/hoffmana_sta/scripts/TMD';
GIS_DIR   = '/kucresis/scratch/hoffmana_sta/vvel/gis';
GL_SHP    = fullfile(GIS_DIR,'GroundingLine_Antarctica_v02.shp');
% LIMA 240 m, used ONLY for the continent inset ("where in Antarctica").
% The map panel and the neighbourhood inset are both REMA.
LIMA_TIF  = fullfile(GIS_DIR,'lima','tiff_90pct','00000-20080319-092059124.tif');
% REMA v2 mosaic hillshade (10 m browse, tile 17_33, EPSG:3031) under the
% map panel and the neighbourhood inset.
REMA_TIF  = fullfile(GIS_DIR,'rema','17_33_10m_v2.0_browse.tif');
CATS      = fullfile(TMD_DIR,'usapdc_601772','CATS2008_v2023.nc');
MP_DIR    = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
OUT_DIR   = '/kucresis/scratch/hoffmana_sta/vvel/figures';
PASS_NAMES = {'EAGER_2022','EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};

STEP_MIN  = 3;      % [min] animation time step
FPS       = 24;     % 40 -> 16 (2.5x slower) -> 24 (1.5x faster again)
N_TIDE_PT = 12;     % points across the array that the mean tide averages over
ZOOM_OUT  = 1.15;   % how far to pull the map back from the survey extent
LIMA_DEC  = 20;     % decimation for the continent inset

addpath(TMD_DIR);
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

%% Collect every distinct pass that went into the multipass products
% The five products are the four out-and-back legs of one line, and leg 1
% appears in two of them, so the same physical traverse would otherwise be
% drawn twice. Dedup on segment plus start time.
P = struct('seg',{},'t0',{},'t1',{},'gps',{},'lat',{},'lon',{});
seen = {};
for i = 1:numel(PASS_NAMES)
  fn = fullfile(MP_DIR, [PASS_NAMES{i} '_multipass03.mat']);
  if ~exist(fn,'file'), continue; end
  L = load(fn,'pass');
  for k = 1:numel(L.pass)
    p = L.pass(k);
    if ~isfield(p,'gps_time') || isempty(p.gps_time), continue; end
    seg = p.param_pass.day_seg;
    key = sprintf('%s_%.0f', seg, round(mean(p.gps_time)));
    if any(strcmp(seen, key)), continue; end
    seen{end+1} = key; %#ok<SAGROW>
    g = p.gps_time(:).';
    P(end+1) = struct('seg', seg, 't0', min(g), 't1', max(g), ...
      'gps', g, 'lat', p.lat(:).', 'lon', p.lon(:).'); %#ok<SAGROW>
  end
  clear L;
end
assert(~isempty(P), 'no passes found');
fprintf('%d distinct traverses from %d products\n', numel(P), numel(PASS_NAMES));

t0 = min([P.t0]); t1 = max([P.t1]);
fprintf('window %s to %s UTC (%.2f days)\n', ...
  datestr(epoch2dn(t0),'yyyy-mm-dd HH:MM'), datestr(epoch2dn(t1),'yyyy-mm-dd HH:MM'), ...
  (t1-t0)/86400);

%% Local map frame, equal scale
alllat = [P.lat]; alllon = [P.lon];
% EPSG:3031 Antarctic Polar Stereographic, in km - the standard Antarctic
% frame, and the one the REMA tile and the MEaSUREs grounding line
% already ship in. NOTE at lon ~168 E the 3031 grid
% runs ~168 deg from local north, so north points roughly DOWN on the map
% panel.
lat0 = mean(alllat,'omitnan'); lon0 = mean(alllon,'omitnan');
ps = projcrs(3031);
[sx_km, sy_km] = ps_km(ps, alllon, alllat);      % every survey fix, in km

%% Mean CATS2008 tide over the array
pad = 0.15*86400;                       % show a little context either side
tt  = (t0-pad) : STEP_MIN*60 : (t1+pad);
tdn = epoch2dn(tt);

idx = round(linspace(1, numel(alllat), N_TIDE_PT));
plat = alllat(idx); plon = alllon(idx);
fprintf('predicting CATS2008 at %d points across the array...\n', numel(plat));
Z = nan(numel(plat), numel(tt));
for q = 1:numel(plat)
  Z(q,:) = tmd_predict(CATS, plat(q), plon(q), tdn, 'h');
end
tide = mean(Z, 1, 'omitnan');
spread = max(max(Z,[],1) - min(Z,[],1));
fprintf('mean tide range %.3f m; max spread across the array %.3f m\n', ...
  max(tide)-min(tide), spread);

%% Context layers: grounding line
% The MEaSUREs grounding line ships in EPSG:3031 metres, and
% the main panel is now in that same grid, so both layers are used as-is - no
% lat/lon round trip, no resampling.
[gx_km, gy_km] = deal([], []);
try
  S = shaperead(GL_SHP);
  GX = []; GY = [];
  for k = 1:numel(S)
    GX = [GX; S(k).X(:); NaN]; GY = [GY; S(k).Y(:); NaN]; %#ok<AGROW>
  end
  gx_km = GX/1e3; gy_km = GY/1e3;    % already EPSG:3031 metres
  near = ~(gx_km < min(sx_km)-20 | gx_km > max(sx_km)+20 | ...
           gy_km < min(sy_km)-20 | gy_km > max(sy_km)+20);
  gx_km(~near) = NaN; gy_km(~near) = NaN;
  fprintf('grounding line: %d vertices, %d near the survey\n', ...
    sum(isfinite(GX)), nnz(near & isfinite(gx_km)));
catch ME
  fprintf('grounding line unavailable (%s)\n', ME.message);
end

%% Figure: build every graphics object ONCE, then only update data
ink = [0 0 0]; ink_soft = [0.45 0.45 0.45];   % ink: all axis text
accent = [0.165 0.471 0.839];   % #2a78d6
hot    = [0.922 0.408 0.204];   % #eb6834

h = figure('Visible','off','Position',[100 100 900 900],'Color','w');
axstyle = {'GridAlpha',0.15,'XColor',ink,'YColor',ink,'Box','off'};

% Panel positions are literal and never touched again
% The survey is a narrow NE-SW strip, so under equal aspect a full-width map
% axes is mostly blank. The map gets a tall narrow box matching the strip and
% the context inset sits beside it, which keeps the two panels one above the
% other while using the space.
axT = axes('parent',h,'Position',[0.10 0.72 0.86 0.21]);
% The survey is a tall narrow NE-SW strip, so the map panel is shaped to
% match it and centred. A full-width panel under equal aspect would show
% mostly empty ice either side rather than the survey.
axM_pos = [0.34 0.07 0.32 0.57];
axM = axes('parent',h,'Position',axM_pos);

% --- top: tide
plot(axT, tdn, tide, '-', 'Color', ink_soft, 'LineWidth', 1.5); hold(axT,'on');
% mark when the radar was actually running
for k = 1:numel(P)
  plot(axT, epoch2dn([P(k).t0 P(k).t1]), [0 0]+min(tide)-0.06, '-', ...
    'Color', accent, 'LineWidth', 3);
end
hTideMark = plot(axT, tdn(1), tide(1), 'o', 'MarkerSize', 11, ...
  'MarkerFaceColor', hot, 'MarkerEdgeColor','w', 'LineWidth', 1.5);
grid(axT,'on'); set(axT, axstyle{:});
xlim(axT, [tdn(1) tdn(end)]);
ylim(axT, [min(tide)-0.12, max(tide)+0.30]);   % headroom for the clock
datetick(axT,'x','mmm dd HH:MM','keeplimits');
ylabel(axT,'Mean tide over the array (m)','Color',ink);
% timestamp, top right of the top panel
hClock = text(axT, 0.985, 0.95, '', 'Units','normalized', ...
  'HorizontalAlignment','right', 'VerticalAlignment','top', ...
  'FontName','Courier', 'FontSize', 13, 'FontWeight','bold', 'Color', ink, ...
  'BackgroundColor', [1 1 1], 'Margin', 3);

% --- bottom: map
hold(axM,'on');
hTrk = gobjects(1, numel(P));
for k = 1:numel(P)
  [tx, ty] = ps_km(ps, P(k).lon, P(k).lat);
  hTrk(k) = plot(axM, tx, ty, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.5);
end
if ~isempty(gx_km) && any(isfinite(gx_km))
  plot(axM, gx_km, gy_km, '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 2);
end
hTrail = plot(axM, NaN, NaN, '-', 'Color', accent, 'LineWidth', 2.5);
hRadar = plot(axM, NaN, NaN, 'o', 'MarkerSize', 13, ...
  'MarkerFaceColor', hot, 'MarkerEdgeColor','w', 'LineWidth', 2);
% Equal scale, with the shown extent matched to the panel's own aspect so a
% narrow NE-SW strip does not leave the panel half empty. ZOOM_OUT pulls back
% a little from the survey so the grounding line has room.
grid(axM,'on'); set(axM, axstyle{:}); axis(axM,'equal');
mx = (min(sx_km) + max(sx_km))/2;
my = (min(sy_km) + max(sy_km))/2;
panel_aspect = (axM_pos(3)*900) / (axM_pos(4)*900);
% grow whichever axis is short until the view matches the panel aspect, so
% the whole survey fits however the strip happens to lie in the 3031 grid
half_y = max(ZOOM_OUT*(max(sy_km)-min(sy_km))/2, ...
             ZOOM_OUT*(max(sx_km)-min(sx_km))/2 / panel_aspect);
half_x = half_y * panel_aspect;
xlim(axM, mx + [-half_x half_x]);
ylim(axM, my + [-half_y half_y]);
% Printed so it can be checked against survey_locator.m, whose panel (d)
% must be this same view for the true-scale box drawn in its panel (c) to
% mean anything.
fprintf(['map panel extent: x %.2f..%.2f km, y %.2f..%.2f km ' ...
         '(%.2f x %.2f km)\n'], mx-half_x, mx+half_x, my-half_y, my+half_y, ...
        2*half_x, 2*half_y);
xlabel(axM, 'Polar stereographic x (km, EPSG:3031)','Color',ink);
ylabel(axM, 'Polar stereographic y (km, EPSG:3031)','Color',ink);

% REMA hillshade under the map, static, drawn once. Tracks flip to white:
% light gray disappears on the hillshade.
if rema_underlay(axM, REMA_TIF)
  set(hTrk, 'Color', 'w', 'LineWidth', 0.7);
  grid(axM, 'off');
end

% ---- locator insets: continent -> region -> neighbourhood -------------
% The same chain survey_locator.m draws as panels (a)-(c), stacked in the
% empty margin to the LEFT of the map panel. Read top to bottom they answer
% "where in Antarctica", then narrow twice into the map beside them.
%
% DELIBERATELY OUTSIDE THE MAP AXES. They used to be drawn inside it, which
% put static furniture on top of the very scene the radar marker and its
% trail move through - the animated content had to be routed around them,
% and any change to the view risked the tracks running underneath. Sitting
% in the figure margin they cannot occlude the moving scene or be occluded
% by it, whatever the map extent does. The margin is dead space otherwise,
% and the insets get more pixels there (135 px rather than 95) so every
% box reads better too.
%
% EVERY BOX IS TRUE SCALE, and it is the intermediate steps that make that
% possible. The old single continent inset had to inflate its box about a
% hundredfold - 5 km on 5482 km is under one pixel - because it tried to
% make the whole jump at once. Split into steps, the SMALLEST box here is
% still ~7% of its inset, so nothing needs inflating and every box means
% exactly what it looks like.
REGION_KM = 400;   % must match survey_locator.m REGION_KM
LOCAL_KM  = 30;    % must match survey_locator.m LOCAL_KM
INS_X = 0.10;   INS_W = 0.15;          % left margin, aligned with the tide
INS_Y = [0.45 0.28 0.11];              % top to bottom, clear of both panels
assert(INS_X + INS_W < axM_pos(1), ...
  'locator insets overlap the map axes: they must stay outside the animated scene');

icx = (min(sx_km)+max(sx_km))/2;  icy = (min(sy_km)+max(sy_km))/2;
REGION = [icx-REGION_KM/2, icx+REGION_KM/2, icy-REGION_KM/2, icy+REGION_KM/2];
LOCAL  = [icx-LOCAL_KM/2,  icx+LOCAL_KM/2,  icy-LOCAL_KM/2,  icy+LOCAL_KM/2];
AOIBOX = [mx-half_x, mx+half_x, my-half_y, my+half_y];

% (a) Antarctica. ALL of it must be in frame: the inset axes is square but
% LIMA spans 5482 x 4657 km, so the window is built explicitly square on
% the longer side and centred rather than left to aspect resolution. The
% padding falls outside the raster, so the axes background is set to
% LIMA's own nodata black and the join is invisible.
try
  gi = geotiffinfo(LIMA_TIF);
  lpx = abs(gi.PixelScale(1));
  LIM = imread(LIMA_TIF, 'PixelRegion', ...
    {[1 LIMA_DEC gi.Height],[1 LIMA_DEC gi.Width]});
  lxk = (gi.BoundingBox(1,1) + [0 gi.Width-1]*lpx)/1e3;
  lyk = (gi.BoundingBox(2,2) - [0 gi.Height-1]*lpx)/1e3;
  axA = axes('parent',h,'Position',[INS_X INS_Y(1) INS_W INS_W]);
  image(axA, lxk, lyk, LIM); set(axA,'YDir','normal'); hold(axA,'on');
  inset_box(axA, REGION, hot);
  lxs = sort(lxk); lys = sort(lyk);
  % 3% margin: sized on the longer side alone the continent lands EXACTLY
  % on the frame, which both looks pinched and makes the coverage test a
  % floating-point coin toss (it failed on the x edges at 0% margin).
  halfL = max(diff(lxs), diff(lys))/2 * 1.03;
  xlim(axA, mean(lxs) + halfL*[-1 1]); ylim(axA, mean(lys) + halfL*[-1 1]);
  style_inset(axA, ink);
  % Explicit black backdrop rather than the axes Color: the padding above
  % and below the raster rendered WHITE through 'Color','k', so the join
  % with LIMA's own nodata black was a visible seam. A patch behind the
  % image does not depend on how the axes background is resolved.
  xl_ = xlim(axA); yl_ = ylim(axA);
  pb = patch(axA, 'XData', xl_([1 2 2 1]), 'YData', yl_([1 1 2 2]), ...
    'FaceColor','k', 'EdgeColor','none');
  uistack(pb, 'bottom');
  tol = 1e-6*halfL;
  covers = lxs(1) >= xl_(1)-tol && lxs(2) <= xl_(2)+tol && ...
           lys(1) >= yl_(1)-tol && lys(2) <= yl_(2)+tol;
  fprintf(['inset (a) Antarctica: %.0f km square window, whole continent ' ...
           'inside: %d, box %.1f%% of frame\n'], 2*halfL, covers, ...
          100*REGION_KM/(2*halfL));
  if ~covers, warning('continent inset does not contain the full LIMA extent'); end

  % (b) the region, from the same LIMA
  [IM2, rx, ry] = lima_window(LIMA_TIF, gi, REGION, 1);
  axB = axes('parent',h,'Position',[INS_X INS_Y(2) INS_W INS_W]);
  image(axB, rx, ry, IM2); set(axB,'YDir','normal'); hold(axB,'on');
  inset_box(axB, LOCAL, hot);
  xlim(axB, REGION(1:2)); ylim(axB, REGION(3:4));
  style_inset(axB, ink);
  fprintf('inset (b) region: %.0f km, box %.1f%% of frame\n', ...
    REGION_KM, 100*LOCAL_KM/REGION_KM);
catch ME
  fprintf('LIMA insets unavailable (%s)\n', ME.message);
end

% (c) the neighbourhood, REMA, carrying the map panel's own extent
axC = axes('parent',h,'Position',[INS_X INS_Y(3) INS_W INS_W]);
hold(axC,'on');
xlim(axC, LOCAL(1:2)); ylim(axC, LOCAL(3:4));
rema_underlay(axC, REMA_TIF);
if ~isempty(gx_km) && any(isfinite(gx_km))
  plot(axC, gx_km, gy_km, '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 1.2);
end
inset_box(axC, AOIBOX, hot);
xlim(axC, LOCAL(1:2)); ylim(axC, LOCAL(3:4));
style_inset(axC, ink);
fprintf('inset (c) neighbourhood: %.0f km, box %.1f%% x %.1f%% of frame\n', ...
  LOCAL_KM, 100*2*half_x/LOCAL_KM, 100*2*half_y/LOCAL_KM);

%% Render
% PREVIEW_ONLY writes ONE frame as a png and stops, for checking layout and
% the inset without waiting on a 200 MB encode. Define it before running:
%   matlab -batch "PREVIEW_ONLY=true; run('.../eastwind_survey_movie.m')"
if ~exist('PREVIEW_ONLY','var'), PREVIEW_ONLY = false; end

active = false(1,numel(tt));
for f = 1:numel(tt)
  active(f) = any(tt(f) >= [P.t0] & tt(f) <= [P.t1]);
end

af = find(active);
prev_fn = fullfile(OUT_DIR,'EAGER_2022_survey_preview.png');

if PREVIEW_ONLY
  % A frame with the radar actually RUNNING, so the trail, the moving marker
  % and the status line are all exercised - the idle state would hide most
  % of what a preview is for
  if isempty(af), frames = 1; else, frames = af(max(1,round(numel(af)/2))); end
  fprintf('PREVIEW_ONLY: rendering frame %d of %d to %s\n', frames, numel(tt), prev_fn);
else
  out_fn = fullfile(OUT_DIR,'EAGER_2022_survey_movie.avi');
  v = VideoWriter(out_fn, 'Motion JPEG AVI'); % no MPEG-4 on this Linux MATLAB
  v.FrameRate = FPS; v.Quality = 92;
  open(v);
  frames = 1:numel(tt);
  fprintf('rendering %d frames to %s\n', numel(tt), out_fn);
  fprintf('%d of %d frames have the radar running (%.0f%%); examples: %s\n', ...
    numel(af), numel(tt), 100*numel(af)/numel(tt), ...
    strjoin(arrayfun(@(v) sprintf('%d',v), af(round(linspace(1,numel(af),5))), ...
    'UniformOutput', false), ', '));
end

for f = frames
  now_t = tt(f);
  set(hTideMark, 'XData', tdn(f), 'YData', tide(f));
  set(hClock, 'String', datestr(epoch2dn(now_t),'yyyy-mm-dd HH:MM UTC'));

  % Which traverse, if any, is running right now
  act = find(now_t >= [P.t0] & now_t <= [P.t1], 1);
  if isempty(act)
    set(hRadar,'XData',NaN,'YData',NaN);
    set(hTrail,'XData',NaN,'YData',NaN);
  else
    p = P(act);
    la = interp1(p.gps, p.lat, now_t, 'linear');
    lo = interp1(p.gps, p.lon, now_t, 'linear');
    [rx, ry] = ps_km(ps, lo, la);
    set(hRadar,'XData',rx,'YData',ry);
    sel = p.gps <= now_t;
    [trx, try_] = ps_km(ps, p.lon(sel), p.lat(sel));
    set(hTrail,'XData',trx,'YData',try_);
  end

  if PREVIEW_ONLY
    print(h, prev_fn, '-dpng', '-r96');
  else
    % print -RGBImage rather than getframe: reliable on a headless display
    writeVideo(v, print(h,'-RGBImage','-r96'));
    if mod(f,200)==0, fprintf('  frame %d/%d\n', f, numel(tt)); end
  end
end
close(h);
if PREVIEW_ONLY
  fprintf('wrote %s\n', prev_fn);
else
  close(v);
  d = dir(out_fn);
  fprintf('wrote %s (%.1f MB, %.0f s at %d fps)\n', out_fn, d.bytes/1e6, numel(tt)/FPS, FPS);
end

%% ========================================================================
function [IM, xk, yk] = lima_window(tif, gi, ext_km, dec)
%LIMA_WINDOW Read LIMA whole (ext_km empty) or over an EPSG:3031 km window.
px = abs(gi.PixelScale(1));
x0 = gi.BoundingBox(1,1); y1 = gi.BoundingBox(2,2);
if isempty(ext_km)
  r = [1 dec gi.Height]; c = [1 dec gi.Width];
  xk = (x0 + [0 gi.Width-1]*px)/1e3;
  yk = (y1 - [0 gi.Height-1]*px)/1e3;
else
  c0 = max(1, floor((ext_km(1)*1e3 - x0)/px));
  c1 = min(gi.Width,  ceil((ext_km(2)*1e3 - x0)/px));
  r0 = max(1, floor((y1 - ext_km(4)*1e3)/px));
  r1 = min(gi.Height, ceil((y1 - ext_km(3)*1e3)/px));
  r = [r0 dec r1]; c = [c0 dec c1];
  xk = (x0 + ([c0 c1]-0.5)*px)/1e3;
  yk = (y1 - ([r0 r1]-0.5)*px)/1e3;
end
IM = imread(tif, 'PixelRegion', {r, c});
end

%% ========================================================================
function inset_box(ax, ext, col)
% True-scale extent box. White underlay first so it reads on both bright
% imagery and dark hillshade.
x = ext([1 2 2 1 1]); y = ext([3 3 4 4 3]);
plot(ax, x, y, '-', 'Color','w', 'LineWidth', 2.6);
plot(ax, x, y, '-', 'Color',col, 'LineWidth', 1.4);
end

%% ========================================================================
function style_inset(ax, ink)
% Common inset styling: true scale, no ticks, thin frame.
set(ax,'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[],'Box','on', ...
  'XColor',ink,'YColor',ink,'LineWidth',1);
end

%% ========================================================================
function dn = epoch2dn(t)
% OPR gps_time is seconds since 1970-01-01
dn = datenum(1970,1,1) + t/86400;
end

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
