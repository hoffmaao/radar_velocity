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
LIMA_TIF  = fullfile(GIS_DIR,'lima','tiff_90pct','00000-20080319-092059124.tif');
CATS      = fullfile(TMD_DIR,'usapdc_601772','CATS2008_v2023.nc');
MP_DIR    = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
OUT_DIR   = '/kucresis/scratch/hoffmana_sta/vvel/figures';
PASS_NAMES = {'EAGER_2022','EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};

STEP_MIN  = 3;      % [min] animation time step
FPS       = 40;
N_TIDE_PT = 12;     % points across the array that the mean tide averages over
LIMA_DEC  = 20;     % decimation when reading LIMA for the continent inset
ZOOM_OUT  = 1.15;   % how far to pull the map back from the survey extent

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
lat0 = mean(alllat,'omitnan'); lon0 = mean(alllon,'omitnan');
xkm = @(lon) (lon-lon0)*111320*cosd(lat0)/1e3;
ykm = @(lat) (lat-lat0)*110540/1e3;

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

%% Context layers: grounding line and LIMA imagery
% The MEaSUREs grounding line ships in EPSG:3031 metres, and so does LIMA, so
% the inset works natively in that grid and only the main panel needs the
% conversion into the local tangent frame.
[gx_km, gy_km] = deal([], []);
try
  S = shaperead(GL_SHP);
  GX = []; GY = [];
  for k = 1:numel(S)
    GX = [GX; S(k).X(:); NaN]; GY = [GY; S(k).Y(:); NaN]; %#ok<AGROW>
  end
  [glat, glon] = projinv(projcrs(3031), GX, GY);
  gx_km = xkm(glon); gy_km = ykm(glat);
  near = ~(gx_km < -20 | gx_km > 20 | gy_km < -20 | gy_km > 20);
  gx_km(~near) = NaN; gy_km(~near) = NaN;
  fprintf('grounding line: %d vertices, %d near the survey\n', ...
    sum(isfinite(GX)), nnz(near & isfinite(gx_km)));
catch ME
  fprintf('grounding line unavailable (%s)\n', ME.message);
end

% Survey extent in EPSG:3031, and the whole continent from LIMA, decimated,
% for a small locator inset. The point of the inset is "where in Antarctica",
% so it is continent scale rather than a regional view.
[aoi_x, aoi_y] = projfwd(projcrs(3031), alllat, alllon);
aoi_x = aoi_x/1e3; aoi_y = aoi_y/1e3;   % km, EPSG:3031
[cx, cy] = projfwd(projcrs(3031), lat0, lon0);
lima_ok = false;
try
  ginfo = geotiffinfo(LIMA_TIF);
  px = abs(ginfo.PixelScale(1));
  x_min = ginfo.BoundingBox(1,1); y_max = ginfo.BoundingBox(2,2);
  IM = imread(LIMA_TIF, 'PixelRegion', ...
    {[1 LIMA_DEC ginfo.Height],[1 LIMA_DEC ginfo.Width]});
  lima_x = (x_min + [0 ginfo.Width-1]*px)/1e3;    % km, EPSG:3031
  lima_y = (y_max - [0 ginfo.Height-1]*px)/1e3;
  lima_ok = true;
  fprintf('LIMA continent inset %d x %d px (decimated %dx)\n', ...
    size(IM,1), size(IM,2), LIMA_DEC);
catch ME
  fprintf('LIMA unavailable (%s)\n', ME.message);
end

%% Figure: build every graphics object ONCE, then only update data
ink = [0.20 0.20 0.20]; ink_soft = [0.45 0.45 0.45];
accent = [0.165 0.471 0.839];   % #2a78d6
hot    = [0.922 0.408 0.204];   % #eb6834

h = figure('Visible','off','Position',[100 100 900 900],'Color','w');
axstyle = {'GridAlpha',0.15,'XColor',ink_soft,'YColor',ink_soft,'Box','off'};

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
title(axT,'CATS2008 tide and EAGER 2022 survey timing','Color',ink);
% timestamp, top right of the top panel
hClock = text(axT, 0.985, 0.95, '', 'Units','normalized', ...
  'HorizontalAlignment','right', 'VerticalAlignment','top', ...
  'FontName','Courier', 'FontSize', 13, 'FontWeight','bold', 'Color', ink, ...
  'BackgroundColor', [1 1 1], 'Margin', 3);

% --- bottom: map
hold(axM,'on');
for k = 1:numel(P)
  plot(axM, xkm(P(k).lon), ykm(P(k).lat), '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.5);
end
if ~isempty(gx_km) && any(isfinite(gx_km))
  plot(axM, gx_km, gy_km, '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 2);
  text(axM, 0.04, 0.985, 'grounding line (MEaSUREs)', 'Units','normalized', ...
    'VerticalAlignment','top', 'FontSize', 9, 'Color', [0.10 0.10 0.10]);
end
hTrail = plot(axM, NaN, NaN, '-', 'Color', accent, 'LineWidth', 2.5);
hRadar = plot(axM, NaN, NaN, 'o', 'MarkerSize', 13, ...
  'MarkerFaceColor', hot, 'MarkerEdgeColor','w', 'LineWidth', 2);
% Equal scale, with the shown extent matched to the panel's own aspect so a
% narrow NE-SW strip does not leave the panel half empty. ZOOM_OUT pulls back
% a little from the survey so the grounding line has room.
grid(axM,'on'); set(axM, axstyle{:}); axis(axM,'equal');
mx = mean([min(xkm(alllon)) max(xkm(alllon))]);
my = mean([min(ykm(alllat)) max(ykm(alllat))]);
half_y = ZOOM_OUT * max(ykm(alllat) - my);
panel_aspect = (axM_pos(3)*900) / (axM_pos(4)*900);
half_x = half_y * panel_aspect;
xlim(axM, mx + [-half_x half_x]);
ylim(axM, my + [-half_y half_y]);
xlabel(axM, sprintf('East of %.4f deg (km)', lon0),'Color',ink);
ylabel(axM, sprintf('North of %.4f deg (km)', lat0),'Color',ink);
hStatus = title(axM,'','Color',ink,'Interpreter','none');

% LIMA context inset, inside the bottom panel so the two panel positions
% the user specified are untouched. Static: it never updates per frame.
if lima_ok
  % Small locator tucked into the bottom left of the map panel, so the figure
  % stays two panels rather than three.
  % Bottom-right corner of the map panel: the survey occupies the left and
  % centre, so the locator never sits over the tracks or the moving marker
  axI = axes('parent',h,'Position',[axM_pos(1)+axM_pos(3)-0.113, axM_pos(2)+0.008, 0.105, 0.105]);
  image(axI, lima_x, lima_y, IM);
  set(axI,'YDir','normal'); hold(axI,'on'); axis(axI,'equal');
  % AOI bounding box rather than a point marker.
  %
  % HONESTY ABOUT THE SCALE: the survey is about 5 km across and this inset
  % spans the whole continent, roughly 5500 km, so a true-scale box would be
  % about a thousandth of the panel - well under one pixel, and invisible.
  % The box is therefore grown to a minimum on-screen size. It marks WHERE
  % the survey is, not how big it is; the main panel carries the real
  % extent, and the box is deliberately square so nobody reads its shape as
  % the survey's footprint.
  % The inset spans ~5500 km in ~130 px, so about 42 km per pixel: a 200 km
  % box came out 5 px and read as a dot. 350 km is ~8 px, which reads as a
  % square outline, and is still a small fraction of the continent.
  MIN_BOX_KM = 350;
  bx0 = mean([min(aoi_x) max(aoi_x)]); by0 = mean([min(aoi_y) max(aoi_y)]);
  bhw = max((max(aoi_x)-min(aoi_x))/2, MIN_BOX_KM/2);
  bhh = max((max(aoi_y)-min(aoi_y))/2, MIN_BOX_KM/2);
  bhw = max(bhw, bhh); bhh = bhw;
  box_x = bx0 + [-bhw  bhw bhw -bhw -bhw];
  box_y = by0 + [-bhh -bhh bhh  bhh -bhh];
  % White underlay first: LIMA is bright in places and dark in others, so a
  % single-colour outline is not reliably legible on its own
  plot(axI, box_x, box_y, '-', 'Color', 'w', 'LineWidth', 3.0);
  plot(axI, box_x, box_y, '-', 'Color', hot, 'LineWidth', 1.6);
  xlim(axI, sort(lima_x)); ylim(axI, sort(lima_y));
  set(axI,'XTick',[],'YTick',[],'Box','on','XColor',ink,'YColor',ink,'LineWidth',1);
  text(axI, 0.5, -0.06, 'LIMA', 'Units','normalized', 'HorizontalAlignment','center', ...
    'VerticalAlignment','top', 'FontSize', 8, 'Color', ink);
end

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
    set(hStatus,'String','radar idle');
  else
    p = P(act);
    la = interp1(p.gps, p.lat, now_t, 'linear');
    lo = interp1(p.gps, p.lon, now_t, 'linear');
    set(hRadar,'XData',xkm(lo),'YData',ykm(la));
    sel = p.gps <= now_t;
    set(hTrail,'XData',xkm(p.lon(sel)),'YData',ykm(p.lat(sel)));
    set(hStatus,'String', sprintf('%s  -  %.0f%% along', p.seg, ...
      100*nnz(sel)/numel(p.gps)));
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
function dn = epoch2dn(t)
% OPR gps_time is seconds since 1970-01-01
dn = datenum(1970,1,1) + t/86400;
end
