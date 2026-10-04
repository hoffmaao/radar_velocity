%TIDAL_RESPONSE_MAP Map view of the tidal response, at the density the data supports.
%
%   ONE PANEL. This figure used to carry a map beside an along-track
%   "pattern test" scatter, but the two encoded the SAME per-block
%   admittance in different coordinates - the map added geography and the
%   grounding line, the scatter added nothing the map did not already
%   have. The scatter is gone. What it did carry that the map could not -
%   the radar-free GPS-curvature prediction - lives in
%   scripts/diagnostics/gps_flexure.m, and the numbers it reported are
%   still printed to stdout here.
%
%   WHAT IS DRAWN. Per-leg tidal response - mm of column change per metre
%   of tide over the top REF_DEPTH metres - as a continuous colour-coded
%   ribbon along each leg's own track, over a REMA v2 10 m hillshade in
%   EPSG:3031 Antarctic Polar Stereographic km, with the MEaSUREs
%   grounding line and the ApRES stations. Diverging colour about zero.
%
%   THE ApRES STATIONS ARE FILLED FROM THE SAME COLOUR SCALE as the
%   ribbons, so the two instruments can be compared by eye. Their values
%   come from scripts/diagnostics/apres_site_admittance.py, which
%   regresses each station's own strain rate on the CATS2008 tide RATE.
%   Two things about that are worth knowing before reading the map:
%
%   1. DEPTH-MATCHED, AND THAT CHANGED THE ANSWER. The station used to
%      carry the ApRES scalar, which is the slope of ApRES's own fit
%      interval - 101-224 m as reported, 115-238 m once its constant
%      er = 3.18 range is mapped onto the radar's firn scale. That is a
%      DEEPER column than the ribbon, and since strain varies strongly
%      with depth here the station rendered light red on a blue ribbon
%      and read as an instrument conflict. It was a depth mismatch. Over
%      38-100 m, the shallowest the ApRES record reaches, the SAME data
%      gives +3.17 +/- 0.34 against the ribbon's +3.21 +/- 0.75 - they
%      agree to 0.04 mm. The station now carries that matched value.
%      ApRES cannot see the firn, so 38-100 m against 0-100 m is as close
%      as the match gets; the colour bar states both.
%
%   2. ONLY STATIONS WITH A USABLE TIDAL SIGNAL APPEAR. Of the four sites
%      with pair_results, one survives: GA04 (n=276, R2 0.75). GA01 is
%      dropped for a tracked bed that drifts ~82 m, GA05 for having 7
%      pairs inside the short-dt window the rate method needs, GA10 for
%      having none. Those are instrument records that cannot answer the
%      question, not places where the ice does nothing.
%
%   ON THE TIDE MODEL. CATS2008 returns NaN at GA01, GA04 and GA05 - they
%   fall inside its grounded mask, and only GA10 is wet where it stands,
%   which is consistent with the first three being at the grounding-line
%   end. Each station is therefore served by the nearest WET cell, 2.0 km
%   away for those three. It makes almost no difference which cell is
%   used: across the four, the predicted series differ by 0.0019 m at most
%   against a 0.546 m amplitude, and pairwise r is 0.999998 or better.
%   Per-station tide is the defensible choice, not a correction.
%
%   ORIENTATION: at lon ~168 E the 3031 grid runs ~168 deg from local
%   north, so north points roughly DOWN in this view. That is what a polar
%   stereographic map looks like at this longitude, and it is the frame
%   the REMA tile and the grounding line already ship in, so neither
%   overlay needs reprojecting.
%
%   SMOOTHED AT WIN_M. The ribbons are a moving window along each leg -
%   inverse-variance weights times a tricube kernel in along-track distance
%   - because a map is read synoptically and the raw 125 m blocks alternate
%   too fast to show a spatial pattern at this scale. line_sections.m uses
%   the same width on its response curve. The two
%   use the SAME kernel width, so a colour means the same number smoothed
%   the same way in both. Change one and the other must follow.
%
%   THE DENSITY, AND WHERE IT COMES FROM. The standing 500 m network build
%   gives 10-11 blocks per leg. This figure reads the 125 m build
%   (CSARP_vvel_net_fine) instead, which gives 39-41 - four times the
%   along-track sampling for about 1.3x the per-block sigma, so pooling a
%   few of them back together beats one coarse block outright. Blocks are
%   drawn AS THEMSELVES - see the no-smoothing note above. STEP_M only
%   subdivides a block for rendering; it is not a resolution, and the
%   independent count is the block count.
%
%   WHY PER-LEG AND NOT ONE POOLED RIBBON. Pooling all legs would collapse
%   four spatially distinct measurements onto one curve, and a map whose
%   four tracks all showed the identical pooled value would imply an
%   agreement it had not tested. Each leg is therefore stacked against
%   itself, and the legs agreeing (or not) is something the reader can see
%   directly. The POOLED profile is still computed - it is what the
%   ApRES-site comparison on stdout is evaluated from.
%
%   THE FOUR RIBBONS ARE ALSO THE ERROR BAR, and that is the main reason
%   to prefer them here. Per-leg sigma runs 1.85 mm/m (GL3) to 3.33 (GL2,
%   the line with the tens-of-metres cross-track baselines), against a
%   colour scale of +/-4 - so a single leg drawn alone would render its
%   own noise at nearly full saturation and a reader could not tell that
%   from signal. Across four legs the distinction is visible without any
%   extra encoding: a feature carried by ALL FOUR at the same PLACE is
%   real, and one leg's excursion that its neighbours do not share is not.
%   Read the map that way. The mid-line thinning and the reversal to
%   thickening at the grounding-line approach pass that test; the
%   isolated patches mid-line, which sit on different legs at different
%   places, do not.
%
%   Run on the server (needs CSARP_vvel_net_fine):
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../tidal_response_map.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef
addpath(fileparts(mfilename('fullpath')));            % grl_figure
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'diagnostics'));  % pass_tide

root    = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
mp_dir  = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
net_dir = fullfile(root,'CSARP_vvel_net_fine');      % 125 m blocks
alt_dir = fullfile(root,'CSARP_vvel_net');           % 500 m, per-line fallback
gis_dir = '/kucresis/scratch/hoffmana_sta/vvel/gis';
out_dir = vdef.figureDir();
% REMA v2 mosaic hillshade (10 m browse, tile 17_33, EPSG:3031) as the map
% background; fetched from the PGC open-data S3 bucket into gis/rema
rema_tif = fullfile(gis_dir,'rema','17_33_10m_v2.0_browse.tif');
% ApRES site positions, EPSG:3031 metres (name in column 3). The file
% lists exactly the four sites that have pair_results - GA01, GA04, GA05,
% GA10 - so it also identifies which of the twelve GA sites are ApRES
% rather than GNSS. A copy lives in the repo at data/gis/.
apres_xy = fullfile(gis_dir,'eastwind_2022_2023_apres_xy.txt');
% ApRES tidal admittance per station, `site mm_per_m sigma R2 n`, written by
% scripts/diagnostics/apres_site_admittance.py. Only stations that pass its
% inclusion rule appear in it; a station absent from the file recorded no
% usable tidal signal and is NOT drawn. A copy lives in the repo at data/gis/.
apres_adm = fullfile(gis_dir,'eastwind_2022_apres_admittance.txt');

% EAGER_2022 is deliberately absent: it is the uncalibrated duplicate of
% GL1 (same physical leg, stale-calibration build, non-reproducible from
% the archive), so plotting it would double-plot the same ice.
PASS_NAMES = {'EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};
REF_DEPTH = 100; MAX_BASELINE = 10;
WIN_M = 500;     % moving-window kernel FULL width, for the drawn ribbons
                 % AND the pooled cross-leg profile the stdout ApRES
                 % comparison reads. At 125 m blocks this holds 4 per leg.
                 %
                 % line_sections.m uses the SAME width, so the two figures
                 % share both the colour scale and the smoothing. Keep them
                 % in step: changing one alone makes identical colours mean
                 % differently-averaged numbers.
STEP_M = 25;     % ribbon sampling step; rendering only, not resolution
CLIM = 4;        % mm/m colour range; values beyond are clipped to the ends
NSIG = 2;        % sigma for full colour saturation; below it, washed out
SHOW_APRES = true;   % stations on; the table and the
                     % projection code stay, so this is a one-word
                     % change to bring them back

PAL.ink = [0 0 0]; PAL.ink_soft = [0.45 0.45 0.45];
PAL.div_neg = [0.698 0.094 0.169]; PAL.div_mid = [0.941 0.937 0.925];
PAL.div_pos = [0.165 0.471 0.839];

%% Per-block network admittance for every line
R = [];
for n = 1:numel(PASS_NAMES)
  S = one_line(PASS_NAMES{n}, net_dir, mp_dir, REF_DEPTH, MAX_BASELINE);
  if isempty(S)
    % A line absent from the fine build still belongs on the map; fall
    % back to the 500 m one and say so, rather than silently dropping a
    % quarter of the survey.
    S = one_line(PASS_NAMES{n}, alt_dir, mp_dir, REF_DEPTH, MAX_BASELINE);
    if isempty(S)
      fprintf('*** %s: no usable product in EITHER build - not drawn\n', ...
        PASS_NAMES{n});
      continue;
    end
    fprintf(['*** %s: absent from the 125 m build, using the 500 m one - ' ...
             'its ribbon is smoother than the others for that reason\n'], ...
      PASS_NAMES{n});
  end
  if isempty(R), R = S; else, R(end+1) = S; end %#ok<AGROW>
  fprintf('%-16s %3d blocks at %5.1f m, adm %+6.2f..%+6.2f, median sigma %.2f\n', ...
    S.name, numel(S.adm), S.step, min(S.adm), max(S.adm), ...
    median(S.adm_std(isfinite(S.adm_std))));
end
assert(numel(R) >= 3, 'need at least three lines');

%% ApRES stations: positions, and the admittance that colours them
% A station is drawn ONLY if it has a usable tidal admittance. Sites that
% recorded too few short pairs for the rate method, or whose tracked bed
% wanders, carry no measurement to show and would otherwise sit on the map
% as an unexplained marker inviting comparison with the ribbon under it.
if SHOW_APRES
  AP = read_apres_xy(apres_xy);
  ADM = read_apres_adm(apres_adm);
else   % reachable by flipping SHOW_APRES above, so the analyzer's dead-branch note is expected
  ap_fields = {'name','x','y','along','off','adm','adm_se','r2','n'}; %#ok<UNRCH>
  AP = reshape(cell2struct(cell(numel(ap_fields), 0), ap_fields, 1), 0, 0);
  ADM = [];
  fprintf('ApRES stations suppressed (SHOW_APRES = false)\n');
end
keep = false(1, numel(AP));
for q = 1:numel(AP)
  [AP(q).along, AP(q).off] = apres_along(R, ps_crs(), AP(q).x, AP(q).y);
  j = find(strcmp({ADM.site}, AP(q).name), 1);
  if isempty(j)
    fprintf('ApRES %-5s along %.2f km - no usable tidal admittance, not drawn\n', ...
      AP(q).name, AP(q).along/1e3);
    continue;
  end
  AP(q).adm = ADM(j).mm; AP(q).adm_se = ADM(j).se;
  AP(q).r2 = ADM(j).r2; AP(q).n = ADM(j).n;
  keep(q) = true;
  fprintf(['ApRES %-5s along %.2f km, %.2f km off: %+.2f +/- %.2f mm/m ' ...
           '(R2 %.2f, n=%d)\n'], AP(q).name, AP(q).along/1e3, ...
    AP(q).off/1e3, AP(q).adm, AP(q).adm_se, AP(q).r2, AP(q).n);
end
AP = AP(keep);
fprintf('%d ApRES station(s) drawn\n', numel(AP));

%% Per-leg ribbon profile, and the pooled one behind the stdout numbers
for i = 1:numel(R)

  % min_n scales with what the line's own spacing can deliver, so a leg
  % that fell back to the 500 m build does not silently produce an empty
  % ribbon against a floor set for 125 m blocks.
  min_n = max(2, min(3, floor(WIN_M/max(R(i).step,eps))));
  [R(i).xg, R(i).vg, R(i).sg, R(i).ng] = ...
    kernel_profile(R(i).along, R(i).adm, R(i).adm_std, WIN_M, STEP_M, min_n);
  % ribbon geometry: the leg's own full-resolution track, sampled at xg
  R(i).glat = interp1(R(i).track_at, R(i).track_lat, R(i).xg, 'linear', NaN);
  R(i).glon = interp1(R(i).track_at, R(i).track_lon, R(i).xg, 'linear', NaN);
  fprintf('%-16s ribbon %3d samples, %d..%d blocks per %.0f m window\n', ...
    R(i).name, nnz(isfinite(R(i).vg)), min(R(i).ng(R(i).ng>0)), ...
    max(R(i).ng), WIN_M);
end

% Pooled: every leg's blocks in one moving window. Not drawn - the map
% shows the legs separately on purpose - but it is the profile the
% ApRES-site comparison below is read from.
pa = []; pv = []; psd = [];
for i = 1:numel(R)
  pa = [pa; R(i).along(:)]; pv = [pv; R(i).adm(:)]; psd = [psd; R(i).adm_std(:)]; %#ok<AGROW>
end
[pxg, pvg, psg, png_] = kernel_profile(pa, pv, psd, WIN_M, STEP_M, 3);
fprintf(['\npooled moving-window profile (%d blocks from %d legs, %.0f m ' ...
         'kernel):\n'], nnz(isfinite(pv)), numel(R), WIN_M);
for xq = 0:500:5000
  [~, k] = min(abs(pxg - xq));
  if isfinite(pvg(k)) && abs(pxg(k)-xq) < STEP_M
    fprintf('  %.2f km  %+6.2f +/- %.2f  (%d blocks)\n', ...
      pxg(k)/1e3, pvg(k), psg(k), png_(k));
  end
end
fprintf('\nradar (pooled, %.0f m kernel) at each ApRES site:\n', WIN_M);
for q = 1:numel(AP)
  if ~isfinite(AP(q).along), continue; end
  [~, k] = min(abs(pxg - AP(q).along));
  if abs(pxg(k) - AP(q).along) > WIN_M/2 || ~isfinite(pvg(k))
    fprintf('  %-5s %.2f km along: outside the profile\n', AP(q).name, ...
      AP(q).along/1e3);
  else
    fprintf('  %-5s %.2f km along (%.2f km off): %+.2f +/- %.2f mm/m (%d blocks)\n', ...
      AP(q).name, AP(q).along/1e3, AP(q).off/1e3, pvg(k), psg(k), png_(k));
  end
end

%% Figure - one portrait map panel
% Portrait because the survey is a ~5 km by ~0.5 km sliver and the map is
% equal-aspect; and large, because the legs are only 120-213 m apart, so
% the px-per-km the panel HEIGHT buys is what decides whether four ribbons
% read as four or as one.
% The panel is narrow as well as tall. Under `axis equal` the visible x
% range is (panel aspect) x (y range), and the y range is fixed at ~6 km by
% the line plus the outlying GA01; so the ONLY way to stop the 0.7 km-wide
% survey floating in an expanse of empty hillshade is to make the panel
% narrower. Height is left long because px-per-km - hence whether four
% ribbons 120-213 m apart read as four - comes from the height.
[h, GRL] = grl_figure(85, 132.8);
set(0,'CurrentFigure',h);
% Margins sized for 8 pt text on an 85 mm canvas, which is the whole
% difficulty of a one-column map: a y tick label like -1313.5 is ~11 mm,
% an eighth of the figure width, and the bottom has to stack tick labels,
% an axis label, the colour bar, ITS tick labels and ITS label. At the old
% pixel proportions the x label ran straight through the bar.
MAP_POS = [0.255 0.205 0.695 0.755];

ps = projcrs(3031);
axm = axes('parent',h,'Position',MAP_POS);
set(0,'CurrentFigure',h); hold(axm,'on');
ndiv = 256; half = ndiv/2;
dmap = [interp1([0 1],[PAL.div_neg; PAL.div_mid], linspace(0,1,half)); ...
        interp1([0 1],[PAL.div_mid; PAL.div_pos], linspace(0,1,ndiv-half))];

% Tracks first, as a thin spine under each ribbon: it keeps a leg visible
% where its profile has no samples (the ends, mostly) instead of letting
% the leg simply stop.
hTrk = gobjects(1, numel(R));
mx = []; my = [];
for i = 1:numel(R)
  [xk, yk] = ps_km(ps, R(i).track_lon, R(i).track_lat);
  hTrk(i) = plot(axm, xk, yk, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
  mx = [mx; xk(:)]; my = [my; yk(:)]; %#ok<AGROW>
end
nsig_tot = 0; nseg_tot = 0;
for i = 1:numel(R)
  [gx, gy] = ps_km(ps, R(i).glon, R(i).glat);
  n = ribbon_line(axm, gx, gy, R(i).vg, R(i).sg, dmap, CLIM, NSIG, ...
    PAL.div_mid, 5);
  nsig_tot = nsig_tot + n; nseg_tot = nseg_tot + max(numel(gx)-1,0);
end
% Ribbon samples overlap, so the raw count overstates the evidence: the
% independent count along a leg is its length divided by the kernel width.
n_indep = 0;
for i = 1:numel(R)
  n_indep = n_indep + (max(R(i).along)-min(R(i).along))/WIN_M;
end
fprintf(['\nsignificance: %d of %d drawn segments reach %.0f sigma (%.0f%%); ' ...
         'about %.0f INDEPENDENT samples across the four legs, since ' ...
         'neighbours share a %.0f m kernel\n'], ...
  nsig_tot, nseg_tot, NSIG, 100*nsig_tot/max(nseg_tot,1), n_indep, WIN_M);
% ApRES stations, FILLED BY THEIR OWN TIDAL ADMITTANCE on the same scale
% as the ribbons, so instrument against instrument is a colour match and
% not a lookup. Larger than the ribbon is wide, with a heavy black edge,
% because it is a point measurement and should not be mistaken for a
% length of ribbon. Read the depth caveat on the colour bar before
% comparing a station against the ribbon it sits on.
for q = 1:numel(AP)
  ci = max(1, min(ndiv, round((AP(q).adm+CLIM)/(2*CLIM)*(ndiv-1))+1));
  % same significance rule as the ribbons, so a station and the ice under
  % it are washed out by the same criterion
  wq = 1;
  if isfinite(AP(q).adm_se) && AP(q).adm_se > 0
    wq = min(1, abs(AP(q).adm)/(NSIG*AP(q).adm_se));
  end
  fc = PAL.div_mid + wq*(dmap(ci,:) - PAL.div_mid);
  plot(axm, AP(q).x/1e3, AP(q).y/1e3, 'o', 'MarkerSize', 11, ...
    'MarkerFaceColor', fc, 'MarkerEdgeColor', PAL.ink, 'LineWidth', 1.6);
  mx(end+1) = AP(q).x/1e3; my(end+1) = AP(q).y/1e3; %#ok<AGROW>
end

% Limits BEFORE the grounding line and the hillshade, both of which clip to
% them. The box is padded out to the panel's own pixel aspect so that
% `daspect [1 1 1]` fits it exactly, leaving no dead map.
fp = get(h,'Position'); mp = get(axm,'Position');
asp = (mp(3)*fp(3)) / (mp(4)*fp(4));
padkm = 0.25;
wx = (max(mx)-min(mx)) + 2*padkm;
wy = (max(my)-min(my)) + 2*padkm;
if wx/wy < asp, wx = asp*wy; else, wy = wx/asp; end
xlim(axm, (min(mx)+max(mx))/2 + wx/2*[-1 1]);
ylim(axm, (min(my)+max(my))/2 + wy/2*[-1 1]);
daspect(axm, [1 1 1]);
pxkm = (mp(4)*fp(4)/2.54*GRL.dpi)/wy;   % figure units are cm now, not px
fprintf('\nmap: %.2f x %.2f km, %.0f px per km (leg spacing 0.12-0.21 km = %.0f-%.0f px)\n', ...
  wx, wy, pxkm, 0.12*pxkm, 0.21*pxkm);

gl_drawn = overlay_gl(axm, gis_dir);
place_site_labels(axm, AP, PAL.ink);
grid(axm,'on');
set(axm,'GridAlpha',0.15,'XColor',PAL.ink,'YColor',PAL.ink,'Box','off');
xlabel(axm, 'Easting (km)','Color',PAL.ink);
ylabel(axm, 'Northing (km)','Color',PAL.ink);
% no panel title: lettering and captions are added in the slide or
% manuscript, same convention as the other figures in this repo
if gl_drawn, fprintf('map: grounding line drawn in black\n'); end

% Colour bar laid HORIZONTALLY under the map: a vertical one costs width,
% which is the scarce dimension for a portrait map. Creating it shrinks
% the parent axes, so the position is put back - the limits above were
% padded to that exact rectangle.
colormap(axm, dmap); caxis(axm, [-CLIM CLIM]);
cb = colorbar(axm, 'southoutside');
set(axm, 'Position', MAP_POS);
set(cb, 'Position', [MAP_POS(1) 0.078 MAP_POS(3) 0.016], ...
  'XColor', PAL.ink, 'YColor', PAL.ink);
% cb.Label, not get(cb,'xlabel'): Label tracks the bar's orientation
% Both intervals are named because they are not identical - but they now
% very nearly are. The stations carry a DEPTH-MATCHED value (38-100 m)
% rather than the ApRES scalar over its own 101-224 m fit interval, which
% is a deeper column and made a station read as though it contradicted the
% ribbon it sits on. ApRES cannot see the firn, so 38 m is as shallow as
% the match goes.
% Two words and the units. What the colour scale means in detail - the
% 0-REF_DEPTH column, and full saturation only at NSIG sigma - belongs in
% the caption, not on an 85 mm axis.
set(cb.Label, 'String', 'Tidal response (mm/m)', 'Color', PAL.ink);

% Imagery under everything, drawn last so the axis limits are final. The
% survey tracks flip to white: light gray disappears on the hillshade.
if rema_underlay(axm, rema_tif)
  set(hTrk, 'Color', 'w', 'LineWidth', 0.7);
  grid(axm, 'off');
end

print(h, fullfile(out_dir,'EAGER_2022_tidal_response_map.png'), '-dpng', sprintf('-r%d', GRL.dpi));
close(h);
fprintf('Wrote %s\n', fullfile(out_dir,'EAGER_2022_tidal_response_map.png'));

%% ========================================================================
function S = one_line(pn, net_dir, mp_dir, REF_DEPTH, MAX_BASELINE)
%ONE_LINE Per-block tidal admittance for one repeat-pass line.
%   Also returns the leg's FULL-RESOLUTION track (2.5 m along-track, from
%   the multipass main pass, which is the geometry the per-block lat/lon
%   were averaged from). The map ribbon is drawn on that track rather than
%   on an interpolation through the block centres.
S = [];
f = dir(fullfile(net_dir, [pn '_vvel_*.mat']));
if isempty(f), return; end
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

P = []; D = []; W = []; along = []; lat = []; lon = []; Nblk = 0; imain = 1;
for q = 1:numel(f)
  tok = regexp(f(q).name, ['^' regexptranslate('escape',pn) '_vvel_(\d+)_(\d+)\.mat$'], ...
    'tokens','once');
  o = load(fullfile(net_dir, f(q).name));
  if ~vdef.pairAligned(o), continue; end
  if max(abs(o.baseline_y)) > MAX_BASELINE, continue; end
  if Nblk == 0
    Nblk = numel(o.S1); along = o.Along_track(:);
    lat = o.Latitude(:); lon = o.Longitude(:);
    if isfield(o,'baseline_main_idx') && ~isempty(o.baseline_main_idx)
      imain = o.baseline_main_idx;
    end
  end
  sv = nan(Nblk,1);
  for b = 1:Nblk
    d = o.depth_blk(:,b); ok = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(ok) || max(d(ok)) < REF_DEPTH, continue; end
    sv(b) = interp1(d(ok), o.dh_blk(ok,b), REF_DEPTH,'linear',NaN)/REF_DEPTH;
  end
  if all(~isfinite(sv)), continue; end
  if isempty(D), D = sv; else, D(:,end+1) = sv; end %#ok<AGROW>
  P(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
  W(end+1) = max(mean(o.coh_blk(:),'omitnan'),1e-3); %#ok<AGROW>
end
if size(P,1) < 10, return; end

N = vdef.invertNetwork(P, D, struct('n_sigma',3,'weights',W,'n_epoch',Np));
tday = (ptime - min(ptime))/86400;
% The tide is CATS2008 at each pass mid-time, with the pass gate, from the
% same helper the flexure inversion uses (scripts/diagnostics/pass_tide.m)
% - NOT the line-mean GPS height, which carries a per-pass height error
% that attenuated this admittance by 15-25% and let a 1.1 m bad pass
% through (see vdef.surfaceAdmittance).
tide = pass_tide(pn, mp_dir);
A = vdef.fitTideAdmittance(N.x, tday, tide);

if imain < 1 || imain > Np, imain = 1; end
S = struct('name',pn, 'along',along, 'lat',lat, 'lon',lon, ...
  'adm',     1e3*REF_DEPTH*A.admittance(:), ...
  'adm_std', 1e3*REF_DEPTH*A.admittance_std(:), ...
  'step',    median(diff(along)), ...
  'track_at',  L.pass(imain).along_track(:), ...
  'track_lat', L.pass(imain).lat(:), ...
  'track_lon', L.pass(imain).lon(:), ...
  'xg',[], 'vg',[], 'sg',[], 'ng',[], 'glat',[], 'glon',[]);
end

%% ========================================================================
function [xg, vg, sg, ng] = kernel_profile(along, val, sigma, W, step, min_n)
%KERNEL_PROFILE Locally weighted inverse-variance mean along track.
%   Blocks inside the window are weighted by 1/sigma^2 TIMES a tricube
%   kernel in along-track distance, so they fade in and out rather than
%   stepping in - a boxcar over ~40 blocks per leg puts the window's own
%   edges into the profile as if they were structure.
%
%   W is the FULL kernel width and is the honest along-track resolution.
%   `step` only sets how finely the result is sampled for drawing;
%   neighbouring samples share most of their blocks and are strongly
%   correlated, so `step` must never be quoted as a resolution.
%
%   sg is the standard error of the WEIGHTED MEAN, sqrt(sum(w^2 s^2))/sum(w).
%   It assumes the blocks are independent, which is the same assumption
%   the disjoint-bin stack it replaces was making.
ok = isfinite(along) & isfinite(val) & isfinite(sigma) & sigma > 0;
along = along(ok); val = val(ok); sigma = sigma(ok);
if isempty(along)
  xg = []; vg = []; sg = []; ng = []; return;
end
xg = (min(along) : step : max(along)).';
vg = nan(size(xg)); sg = nan(size(xg)); ng = zeros(size(xg));
for k = 1:numel(xg)
  u  = abs(along - xg(k)) / (W/2);
  in = u < 1;
  if nnz(in) < min_n, continue; end
  kern = (1 - u(in).^3).^3;                 % tricube
  w  = kern ./ sigma(in).^2;
  sw = sum(w);
  if sw <= 0, continue; end
  vg(k) = sum(w .* val(in)) / sw;
  sg(k) = sqrt(sum((w.^2) .* (sigma(in).^2))) / sw;
  ng(k) = nnz(in);
end
end

%% ========================================================================
function ps = ps_crs()
%PS_CRS The map projection, as a function so helpers need no plumbing.
persistent p
if isempty(p), p = projcrs(3031); end
ps = p;
end

%% ========================================================================
function ADM = read_apres_adm(fn)
%READ_APRES_ADM Per-station ApRES tidal admittance, `site mm sigma R2 n`.
%   Comment lines start with #. A station absent from the file has no
%   usable admittance and is not drawn; an absent FILE means no station is
%   drawn at all, which is the safe failure - a bare marker on a
%   value-coloured map invites a comparison it cannot support.
ADM = struct('site',{},'mm',{},'se',{},'r2',{},'n',{});
if ~exist(fn,'file')
  fprintf('ApRES admittance table not found (%s) - no stations drawn.\n', fn);
  return;
end
try
  fid = fopen(fn,'r');
  C = textscan(fid, '%s %f %f %f %f', 'CommentStyle', '#');
  fclose(fid);
  for k = 1:numel(C{1})
    ADM(end+1) = struct('site',C{1}{k}, 'mm',C{2}(k), 'se',C{3}(k), ...
      'r2',C{4}(k), 'n',C{5}(k)); %#ok<AGROW>
  end
  fprintf('ApRES admittance read for: %s\n', strjoin({ADM.site}, ', '));
catch ME
  fprintf('Could not read %s (%s) - no stations drawn.\n', fn, ME.message);
end
end

%% ========================================================================
function [along_m, off_m] = apres_along(R, ps, X, Y)
%APRES_ALONG Along-track position of an EPSG:3031 point (metres).
%   The point is PROJECTED onto the leg track, not snapped to the nearest
%   sample: it returns the along-track coordinate in the radar's own frame
%   plus the perpendicular distance, which matters here because the sites
%   sit a few hundred metres OFF the line rather than on it. Runs on the
%   FULL-RESOLUTION track, so the old caution about approximating the line
%   by its bounding-box diagonal (wrong by up to 0.3 km along, 0.5 km
%   across) no longer applies.
along_m = NaN; off_m = NaN; best = inf;
for i = 1:numel(R)
  [bx, by] = ps_km(ps, R(i).track_lon, R(i).track_lat);
  bx = bx*1e3; by = by*1e3;
  [~, k] = min(hypot(bx - X, by - Y));
  for kk = [k-1, k+1]
    if kk < 1 || kk > numel(bx), continue; end
    vx = bx(kk)-bx(k); vy = by(kk)-by(k);
    L2 = vx^2 + vy^2;
    if L2 == 0, continue; end
    t = max(0, min(1, ((X-bx(k))*vx + (Y-by(k))*vy)/L2));
    d = hypot(bx(k) + t*vx - X, by(k) + t*vy - Y);
    if d < best
      best = d; off_m = d;
      along_m = R(i).track_at(k) + t*(R(i).track_at(kk) - R(i).track_at(k));
    end
  end
end
end

%% ========================================================================
function [xg, vg, sg, ng] = block_profile(along, val, sigma, step)
%BLOCK_PROFILE The blocks themselves, sampled piecewise-constant for drawing.
%   Returns a fine grid carrying each block's own value and sigma over the
%   block's own extent - nearest-neighbour, so a sample belongs to whichever
%   block it falls in and nothing is interpolated across a block boundary.
%   ng is 1 everywhere a block exists: no averaging is happening, and the
%   name is kept only so callers need not special-case it.
ok = isfinite(along) & isfinite(val) & isfinite(sigma) & sigma > 0;
along = along(ok); val = val(ok); sigma = sigma(ok);
if isempty(along)
  xg = []; vg = []; sg = []; ng = []; return;
end
xg = (min(along) : step : max(along)).';
vg = interp1(along, val,   xg, 'nearest', NaN);
sg = interp1(along, sigma, xg, 'nearest', NaN);
ng = double(isfinite(vg));
end
