%SURVEY_LOCATOR Honest nested zoom from Antarctica to the EAGER survey.
%
%   WHY THIS EXISTS. No single locator box can answer "where is this?" to
%   scale: the survey is about 5 km across against a 5500 km continent, so
%   a true-scale box on a continent panel would be under one pixel. The
%   answer has to be built in STEPS, and this figure is where that
%   true-scale nested chain is defined. The movie's locator insets draw
%   the same chain - its comments cross-reference this file - after its
%   old single continent inset, which inflated its box about a hundredfold
%   just to be visible, was replaced by these steps.
%
%   Every box here is TRUE SCALE; the readability an inflated box would
%   buy is bought instead with intermediate panels and connector lines.
%   The zoom factors are ~14x, ~13x, then ~5x, which is what it takes to
%   get from a continent to a 5 km line without a step where the box
%   vanishes.
%
%     (a) Antarctica, LIMA 240 m.        Box = (b), 400 km  -> 7% of frame
%     (b) Ross Island region, LIMA.      Box = (c),  30 km  -> 8%
%     (c) The survey's neighbourhood, REMA 10 m hillshade.
%                                        Box = (d), 3.3x5.9 km -> 20%
%     (d) The survey itself, over REMA, with the passes and the MEaSUREs
%         grounding line.
%
%   Four panels rather than three because the honest box sizes demand it.
%   Going straight from the continent to a 100 km frame leaves a 1.8% box,
%   and the REMA tile is a poor middle step anyway: the survey sits at its
%   top edge, so most of that frame is featureless shelf. Two LIMA steps
%   then two REMA steps keeps every box between 7% and 20% of its frame -
%   visible on its own, before the connector lines help at all.
%
%   PANEL (d) IS THE MOVIE'S MAP PANEL, to the kilometre. Its extent is
%   computed with the same formula and the same constants the movie uses
%   (ZOOM_OUT, and the movie's map-panel aspect), so the box drawn in (c)
%   really is the frame the tracks are drawn in. Both scripts print the
%   extent they computed, so a divergence is visible rather than silent -
%   see the CONSISTENCY block below, which is the thing to check if these
%   two figures ever stop agreeing.
%
%   Run on the server, where LIMA, the REMA tile and the products live:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../survey_locator.m')"

GIS_DIR   = '/kucresis/scratch/hoffmana_sta/vvel/gis';
GL_SHP    = fullfile(GIS_DIR,'GroundingLine_Antarctica_v02.shp');
LIMA_TIF  = fullfile(GIS_DIR,'lima','tiff_90pct','00000-20080319-092059124.tif');
REMA_TIF  = fullfile(GIS_DIR,'rema','17_33_10m_v2.0_browse.tif');
MP_DIR    = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef
addpath(fileparts(mfilename('fullpath')));            % grl_figure, map helpers
OUT_DIR   = vdef.figureDir();
PASS_NAMES = vdef.surveyLines();   % EAGER_2022 duplicates GL1

% ---- CONSISTENCY WITH eastwind_survey_movie.m -------------------------
% These three MUST match the movie or panel (d) stops being the movie's
% map panel and the box in (c) becomes a lie of exactly the kind this
% figure exists to avoid. Both scripts print the extent they derive; if
% they ever disagree, this block is why.
ZOOM_OUT     = 1.15;                    % movie: ZOOM_OUT
MOVIE_AXM    = [0.34 0.07 0.32 0.57];   % movie: axM_pos
MOVIE_FIG_PX = 900;                     % movie: figure Position 3rd/4th
% -----------------------------------------------------------------------

LIMA_DEC  = 12;   % continent panel decimation (240 m -> ~2.9 km/px)
REGION_KM = 400;  % (b) half-continent to regional
LOCAL_KM  = 30;   % (c) regional to neighbourhood. 30 not 32: the REMA tile
                  % stops at y = -1300 km and the survey sits near that
                  % edge, so a larger centred window would run off the tile.

ink = [0 0 0]; ink_soft = [0.45 0.45 0.45];
hot = [0.922 0.408 0.204];              % the box colour, as in the movie

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end
ps = projcrs(3031);

%% Every distinct traverse, exactly as the movie collects them
P = struct('lat',{},'lon',{});
seen = {};
for i = 1:numel(PASS_NAMES)
  fn = fullfile(MP_DIR, [PASS_NAMES{i} '_multipass03.mat']);
  if ~exist(fn,'file'), continue; end
  L = load(fn,'pass');
  for k = 1:numel(L.pass)
    p = L.pass(k);
    if ~isfield(p,'gps_time') || isempty(p.gps_time), continue; end
    key = sprintf('%s_%.0f', p.param_pass.day_seg, round(mean(p.gps_time)));
    if any(strcmp(seen, key)), continue; end
    seen{end+1} = key; %#ok<SAGROW>
    P(end+1) = struct('lat', p.lat(:).', 'lon', p.lon(:).'); %#ok<SAGROW>
  end
  clear L;
end
assert(~isempty(P), 'no passes found');
fprintf('%d distinct traverses\n', numel(P));

alllat = [P.lat]; alllon = [P.lon];
[sx_km, sy_km] = ps_km(ps, alllon, alllat);

%% The movie's map extent, by the movie's own formula
mx = (min(sx_km) + max(sx_km))/2;
my = (min(sy_km) + max(sy_km))/2;
panel_aspect = (MOVIE_AXM(3)*MOVIE_FIG_PX) / (MOVIE_AXM(4)*MOVIE_FIG_PX);
half_y = max(ZOOM_OUT*(max(sy_km)-min(sy_km))/2, ...
             ZOOM_OUT*(max(sx_km)-min(sx_km))/2 / panel_aspect);
half_x = half_y * panel_aspect;
AOI = [mx-half_x, mx+half_x, my-half_y, my+half_y];
fprintf(['survey view (must equal the movie''s map panel): ' ...
         'x %.2f..%.2f km, y %.2f..%.2f km  (%.2f x %.2f km)\n'], ...
        AOI(1), AOI(2), AOI(3), AOI(4), AOI(2)-AOI(1), AOI(4)-AOI(3));

%% Grounding line, already EPSG:3031 metres
gx = []; gy = [];
try
  S = shaperead(GL_SHP);
  for k = 1:numel(S)
    gx = [gx; S(k).X(:); NaN]; gy = [gy; S(k).Y(:); NaN]; %#ok<AGROW>
  end
  gx = gx/1e3; gy = gy/1e3;
catch ME
  fprintf('grounding line unavailable (%s)\n', ME.message);
end

%% Figure: four panels, each sized to its own data aspect so nothing is
%% letterboxed and every box stays true scale
[h, GRL] = grl_figure(170, 59.5);
axA = axes('parent',h,'Position',[0.015 0.20 0.250 0.64]);   % continent
axB = axes('parent',h,'Position',[0.300 0.20 0.210 0.64]);   % region
axC = axes('parent',h,'Position',[0.545 0.20 0.210 0.64]);   % neighbourhood
axD = axes('parent',h,'Position',[0.815 0.12 0.150 0.80]);   % the survey

cx = (min(sx_km)+max(sx_km))/2; cy = (min(sy_km)+max(sy_km))/2;
REGION = [cx-REGION_KM/2, cx+REGION_KM/2, cy-REGION_KM/2, cy+REGION_KM/2];
LOCAL  = [cx-LOCAL_KM/2,  cx+LOCAL_KM/2,  cy-LOCAL_KM/2,  cy+LOCAL_KM/2];

%% (a) Antarctica, and (b) the region, both from LIMA
gi = geotiffinfo(LIMA_TIF);
[IM, lx, ly] = lima_window(LIMA_TIF, gi, [], LIMA_DEC);
image(axA, lx, ly, IM); set(axA,'YDir','normal'); hold(axA,'on');
fprintf('LIMA continent %d x %d px (decimated %dx)\n', size(IM,1), size(IM,2), LIMA_DEC);

[IM2, rx, ry] = lima_window(LIMA_TIF, gi, REGION, 1);
image(axB, rx, ry, IM2); set(axB,'YDir','normal'); hold(axB,'on');
fprintf('LIMA region %d x %d px over %.0f km\n', size(IM2,1), size(IM2,2), REGION_KM);

%% (c) the neighbourhood and (d) the survey, both REMA
axes(axC); hold(axC,'on'); %#ok<LAXES>
xlim(axC, LOCAL(1:2)); ylim(axC, LOCAL(3:4));
rema_underlay(axC, REMA_TIF);
if ~isempty(gx), plot(axC, gx, gy, '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 1.4); end
xlim(axC, LOCAL(1:2)); ylim(axC, LOCAL(3:4));

axes(axD); hold(axD,'on'); %#ok<LAXES>
hTrk = gobjects(1,numel(P));
for k = 1:numel(P)
  [tx_, ty_] = ps_km(ps, P(k).lon, P(k).lat);
  hTrk(k) = plot(axD, tx_, ty_, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.8);
end
if ~isempty(gx)
  plot(axD, gx, gy, '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 2);
end
xlim(axD, AOI(1:2)); ylim(axD, AOI(3:4));
if rema_underlay(axD, REMA_TIF)
  set(hTrk, 'Color', 'w', 'LineWidth', 1.0);
end

%% True-scale boxes, each marking the NEXT panel
box_on(axA, REGION, hot);
box_on(axB, LOCAL,  hot);
box_on(axC, AOI,    hot);

%% Axes styling
for ax = [axA axB axC axD]
  set(ax, 'DataAspectRatio',[1 1 1], 'Box','on', ...
    'XColor',ink,'YColor',ink,'LineWidth',0.8, 'FontSize', 8);
end
xlim(axA, lx); ylim(axA, sort(ly));
set(axA,'XTick',[],'YTick',[]);      % continent: a locator, not a chart
xlabel(axD,'Polar stereographic x (km, EPSG:3031)','Color',ink,'FontSize', 8);
ylabel(axD,'Polar stereographic y (km, EPSG:3031)','Color',ink,'FontSize', 8);

%% Connector lines, box corners to the next panel's corners
connect(h, axA, REGION, axB, hot);
connect(h, axB, LOCAL,  axC, hot);
connect(h, axC, AOI,    axD, hot);

out = fullfile(OUT_DIR,'EAGER_2022_survey_locator.png');
print(h, out, '-dpng', sprintf('-r%d', GRL.dpi));
close(h);
fprintf('wrote %s\n', out);

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
function box_on(ax, ext, col)
% True-scale extent box: white underlay so it reads on both bright and
% dark imagery, then the colour on top.
x = ext([1 2 2 1 1]); y = ext([3 3 4 4 3]);
plot(ax, x, y, '-', 'Color','w', 'LineWidth', 3.0);
plot(ax, x, y, '-', 'Color',col, 'LineWidth', 1.6);
end

%% ========================================================================
function connect(h, ax_from, ext, ax_to, col)
% Leader lines from a box in ax_from to the whole of ax_to, in figure
% coordinates. Drawn from the box's right-hand corners to the target
% panel's left-hand corners, which is the standard nested-locator idiom
% and reads correctly for panels laid out left to right. Both ends use the
% axes' DRAWN plot box, not its Position rectangle: DataAspectRatio
% [1 1 1] letterboxes the plot box inside Position, and leaders anchored
% to Position would visibly miss the true-scale corners this figure
% exists to register.
p1 = data2fig(ax_from, ext(2), ext(4));    % box top right
p2 = data2fig(ax_from, ext(2), ext(3));    % box bottom right
q  = plotbox(ax_to);
for pr = {[p1; q(1), q(2)+q(4)], [p2; q(1), q(2)]}
  a = pr{1};
  annotation(h,'line', [a(1,1) a(2,1)], [a(1,2) a(2,2)], ...
    'Color', col, 'LineWidth', 0.8, 'LineStyle',':');
end
end

%% ========================================================================
function p = data2fig(ax, x, y)
% Data coordinates to normalised figure coordinates, through the axes'
% drawn plot box (see plotbox).
pos = plotbox(ax); xl = xlim(ax); yl = ylim(ax);
p = [pos(1) + (x-xl(1))/diff(xl)*pos(3), pos(2) + (y-yl(1))/diff(yl)*pos(4)];
end

%% ========================================================================
function pb = plotbox(ax)
% The rectangle the axes actually draw in, in normalised figure units.
% With DataAspectRatio [1 1 1] and manual limits, MATLAB shrinks the plot
% box to the data aspect and centres it inside Position - the letterboxed
% margins are what leaders must not treat as axes.
pos = get(ax,'Position');
fp  = get(ancestor(ax,'figure'),'Position');            % pixels
w = pos(3)*fp(3); ht = pos(4)*fp(4);
xl = xlim(ax); yl = ylim(ax);
ar = abs(diff(yl)/diff(xl));                            % units are equal km
if w*ar <= ht
  hh = w*ar;                                            % letterboxed top/bottom
  pb = [pos(1), pos(2) + (ht-hh)/2/fp(4), pos(3), hh/fp(4)];
else
  ww = ht/ar;                                           % letterboxed left/right
  pb = [pos(1) + (w-ww)/2/fp(3), pos(2), ww/fp(3), pos(4)];
end
end
