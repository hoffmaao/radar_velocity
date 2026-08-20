%ADMITTANCE_MAP Map view of the tidal response, against the GPS-predicted pattern.
%
%   Two panels:
%
%   (a) MAP. Per-block tidal response (mm of column change per metre of
%       tide, top 100 m, network inversion) for the four calibrated lines,
%       in EPSG:3031 Antarctic Polar Stereographic km with the MEaSUREs
%       grounding line, over a REMA v2 10 m hillshade. Diverging colour
%       about zero. This is where a spatial pattern would show. Note that
%       at lon ~168 E the 3031 grid runs ~168 deg from local north, so
%       north points roughly DOWN in this view.
%
%   (b) THE PATTERN TEST. All calibrated lines' per-block responses
%       against along-track position, with their inverse-variance stack -
%       and over it, the response PREDICTED FROM GPS ALONE: plate bending
%       says the tidal strain is set by the curvature of the deflection
%       profile, eps_zz = nu/(1-nu) * (H/2) * a''(x) * (1 - z/H), and
%       a(x) is measured by the GPS with no radar involved. If the ApRES
%       and radar magnitudes are right, the stacked radar profile should
%       sit on the GPS-curvature prediction.
%
%   THE ApRES SITE POSITIONS WERE RECOVERED 19 Aug 2026 and the four
%   sites with pair_results (GA01, GA04, GA05, GA10) are now drawn on the
%   map, with GA04 - the only site whose tidal admittance is trustworthy,
%   GA01's tracked bed drifts 81 m - placed at its OWN along-track
%   position in panel (b) instead of the band across the axis that the
%   unknown position used to force. GA04 projects to 4.70 km along, 0.07 km
%   OFF the line, i.e. in the grounding-line approach rather than on the
%   flat floating section the earlier write-up assumed. Those numbers come
%   from apres_along() below, which projects each site onto the per-block
%   track; approximating the line by its bounding-box diagonal instead was
%   wrong by up to 0.3 km along and 0.5 km across.
%
%   Run on the server (needs CSARP_vvel_net):
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../admittance_map.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef

root    = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
mp_dir  = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
net_dir = fullfile(root,'CSARP_vvel_net');
gis_dir = '/kucresis/scratch/hoffmana_sta/vvel/gis';
out_dir = '/kucresis/scratch/hoffmana_sta/vvel/figures';
% REMA v2 mosaic hillshade (10 m browse, tile 17_33, EPSG:3031) as the map
% background; fetched from the PGC open-data S3 bucket into gis/rema
rema_tif = fullfile(gis_dir,'rema','17_33_10m_v2.0_browse.tif');
% ApRES site positions, EPSG:3031 metres (name in column 3). The file
% lists exactly the four sites that have pair_results - GA01, GA04, GA05,
% GA10 - so it also identifies which of the twelve GA sites are ApRES
% rather than GNSS. A copy lives in the repo at data/gis/.
apres_xy = fullfile(gis_dir,'eastwind_2022_2023_apres_xy.txt');

% EAGER_2022 is deliberately absent: it is the uncalibrated duplicate of
% GL1 (same physical leg, stale-calibration build, non-reproducible from
% the archive), so plotting it would double-plot the same ice.
PASS_NAMES = {'EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};
CALIB = [true true true true];
% Palette/marker index into PAL.cat/PAL.cat_mk, fixed to product identity
% (GL1 orange square, GL2 aqua triangle, GL3 yellow diamond, GL4 magenta
% inverted-triangle) - do not re-index by position.
PAL_IDX = [2 3 4 5];
REF_DEPTH = 100; MAX_BASELINE = 10; BLOCK = 200;
H_ICE = 300; NU = 0.33;
APRES_MM = -1.24; APRES_SE = 0.04;    % rate method, apres_rate_check.py

PAL.cat = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; ...
           0.929 0.631 0.000; 0.910 0.482 0.643];
PAL.cat_mk = {'o','s','^','d','v'};
% ink is pure black: axis labels, tick labels and axis lines all read black
PAL.ink = [0 0 0]; PAL.ink_soft = [0.45 0.45 0.45];
PAL.div_neg = [0.698 0.094 0.169]; PAL.div_mid = [0.941 0.937 0.925];
PAL.div_pos = [0.165 0.471 0.839];
PAL.expect = [0.35 0.35 0.35];

%% Per-block network admittance + GPS a(x) for every line
R = [];
for n = 1:numel(PASS_NAMES)
  S = one_line(PASS_NAMES{n}, net_dir, mp_dir, REF_DEPTH, MAX_BASELINE, BLOCK);
  if isempty(S), continue; end
  S.calib = CALIB(n);
  % palette colour/marker stay tied to the product identity even when an
  % earlier product was skipped
  S.pi = PAL_IDX(n);
  if isempty(R), R = S; else, R(end+1) = S; end %#ok<AGROW>
  fprintf('%-16s %2d blocks, adm %.2f..%.2f mm/m\n', S.name, numel(S.adm), ...
    min(S.adm), max(S.adm));
end
assert(numel(R) >= 4, 'need the lines');

%% ApRES site positions
% Located 19 Aug 2026; before that the site position was unrecorded and
% the ApRES value had to be drawn as a band across the whole axis.
AP = read_apres_xy(apres_xy);
for q = 1:numel(AP)
  [AP(q).along, AP(q).off] = apres_along(R, ps_crs(), AP(q).x, AP(q).y);
  fprintf('ApRES %-5s along %.2f km, %.2f km off the line\n', ...
    AP(q).name, AP(q).along/1e3, AP(q).off/1e3);
end

%% GPS-curvature prediction along track
% mean a(x) over the calibrated lines, quartic fit, analytic curvature
ax_all = []; xx_all = [];
for i = 1:numel(R)
  if ~R(i).calib, continue; end
  ax_all = [ax_all; R(i).ax_gps(:)]; %#ok<AGROW>
  xx_all = [xx_all; R(i).along(:)];  %#ok<AGROW>
end
ok = isfinite(ax_all) & isfinite(xx_all);
x0 = mean(xx_all(ok));
% Curvature from BOTH a cubic and a quartic fit, drawn as a range: the
% second derivative of a polynomial is least constrained at the profile
% ends, so a single confident curve there would overstate the prediction.
% The display is also trimmed half a kilometre at each end for the same
% reason.
xg = linspace(min(xx_all(ok))+500, max(xx_all(ok))-500, 200);
pred = nan(2, numel(xg));
for dg = 3:4
  pp = polyfit(xx_all(ok)-x0, ax_all(ok), dg);
  app = polyval(polyder(polyder(pp)), xg-x0);    % a''(x) [1/m^2]
  pred(dg-2,:) = 1e3 * (NU/(1-NU)) * (H_ICE/2) .* app ...
    * (1 - REF_DEPTH/H_ICE) * REF_DEPTH;
end
pred_lo = min(pred,[],1); pred_hi = max(pred,[],1);
pred_mm = mean(pred,1);
fprintf('GPS-curvature prediction: %.2f .. %.2f mm/m along the line\n', ...
  min(pred_mm), max(pred_mm));

%% Figure
h = figure('Visible','off','Position',[100 100 1180 560],'Color','w');
set(0,'CurrentFigure',h);
axst = {'GridAlpha',0.15,'XColor',PAL.ink,'YColor',PAL.ink,'Box','off'};

% EPSG:3031 Antarctic Polar Stereographic, in km. This is the standard
% Antarctic frame, and it is the frame the REMA tile and the MEaSUREs
% grounding line already ship in, so neither overlay needs reprojecting -
% the hillshade drops straight in with no resampling.
%
% ORIENTATION: at lon ~168 E the 3031 grid is rotated ~168 deg from local
% north, so north points roughly DOWN in this view. That is what a polar
% stereographic map looks like at this longitude; the previous local
% tangent frame was north-up but is not a projection anyone else's data
% is in.
ps = projcrs(3031);

% (a) map
axm = axes('parent',h,'Position',[0.06 0.12 0.40 0.78]);
set(0,'CurrentFigure',h); hold(axm,'on');
CLIM = 4;   % mm/m colour range; values beyond are clipped to the ends
ndiv = 256; half = ndiv/2;
dmap = [interp1([0 1],[PAL.div_neg; PAL.div_mid], linspace(0,1,half)); ...
        interp1([0 1],[PAL.div_mid; PAL.div_pos], linspace(0,1,ndiv-half))];
hTrk = gobjects(1, numel(R));
for i = 1:numel(R)
  [xk, yk] = ps_km(ps, R(i).lon, R(i).lat);
  hTrk(i) = plot(axm, xk, yk, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
  for b = 1:numel(R(i).adm)
    if ~isfinite(R(i).adm(b)), continue; end
    ci = max(1, min(ndiv, round((R(i).adm(b)+CLIM)/(2*CLIM)*(ndiv-1))+1));
    mk = PAL.cat_mk{R(i).pi};
    ec = PAL.ink_soft; if ~R(i).calib, ec = [0.75 0.75 0.75]; end
    plot(axm, xk(b), yk(b), mk, 'MarkerSize', 10, ...
      'MarkerFaceColor', dmap(ci,:), 'MarkerEdgeColor', ec, 'LineWidth', 0.6);
  end
end
% ApRES sites: white fill, black edge, so they read against both the
% hillshade and the value-coloured radar markers
for q = 1:numel(AP)
  plot(axm, AP(q).x/1e3, AP(q).y/1e3, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor','w', 'MarkerEdgeColor', PAL.ink, 'LineWidth', 1.4);
  text(axm, AP(q).x/1e3 + 0.12, AP(q).y/1e3, AP(q).name, 'FontSize', 8, ...
    'Color', PAL.ink, 'VerticalAlignment','middle');
end
gl_drawn = overlay_gl(axm, gis_dir);
grid(axm,'on'); set(axm, axst{:}); axis(axm,'equal');
xlabel(axm, 'Polar stereographic x (km, EPSG:3031)','Color',PAL.ink);
ylabel(axm, 'Polar stereographic y (km, EPSG:3031)','Color',PAL.ink);
% no panel title: the (a)/(b) lettering and captions are added separately
% in the slide or manuscript, same convention as the method schematic
if gl_drawn, fprintf('map: grounding line drawn in black\n'); end
colormap(axm, dmap); caxis(axm, [-CLIM CLIM]);
cb = colorbar(axm);
set(get(cb,'ylabel'),'string','mm per m of tide (top 100 m)','Color',PAL.ink);
set(cb,'XColor',PAL.ink,'YColor',PAL.ink);

% Imagery under everything, drawn last so the axis limits are final. The
% survey tracks flip to white: light gray disappears on the hillshade.
if rema_underlay(axm, rema_tif)
  set(hTrk, 'Color', 'w', 'LineWidth', 0.7);
  grid(axm, 'off');
end

% (b) pattern test
axp = axes('parent',h,'Position',[0.58 0.12 0.385 0.78]);
set(0,'CurrentFigure',h); hold(axp,'on');
xmax = 0;
for i = 1:numel(R), xmax = max(xmax, max(R(i).along)/1e3); end
plot(axp, [0 xmax], [0 0], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
hl = []; lb = {};
for i = 1:numel(R)
  if ~R(i).calib, continue; end
  hh = plot(axp, R(i).along/1e3, R(i).adm, PAL.cat_mk{R(i).pi}, 'MarkerSize', 5, ...
    'MarkerFaceColor', PAL.cat(R(i).pi,:), 'MarkerEdgeColor','w', 'LineWidth', 0.6);
  hl(end+1) = hh; lb{end+1} = strrep(R(i).name,'EAGER_2022_',''); %#ok<AGROW>
end
% inverse-variance stack in 500 m bins across the calibrated lines
edges = 0:0.5:ceil(xmax*2)/2;
xc = edges(1:end-1) + 0.25;
sm = nan(size(xc)); ss = nan(size(xc));
for k = 1:numel(xc)
  v = []; w = [];
  for i = 1:numel(R)
    if ~R(i).calib, continue; end
    al = R(i).along(:)/1e3; av = R(i).adm(:); as_ = R(i).adm_std(:);
    % everything as columns: a row/column mix here broadcasts to a matrix
    inb = al >= edges(k) & al < edges(k+1) & isfinite(av) & isfinite(as_) & as_ > 0;
    v = [v av(inb).']; w = [w 1./as_(inb).'.^2]; %#ok<AGROW>
  end
  if numel(v) >= 2
    sm(k) = sum(w.*v)/sum(w); ss(k) = sqrt(1/sum(w));
  end
end
okb = isfinite(sm);
fprintf('\nstacked response by 500 m bin (mm per m of tide, top 100 m):\n');
for k = 1:numel(xc)
  if okb(k), fprintf('  %.2f km  %+6.2f +/- %.2f\n', xc(k), sm(k), ss(k)); end
end
% The POSITION-MATCHED comparison: the radar stack interpolated to each
% ApRES site's own along-track position. Before the sites were located
% (19 Aug 2026) the only available comparison was ApRES against the radar
% LINE MEAN, which is a different quantity wherever the profile varies.
% A site past the last stacked bin centre cannot be interpolated to; report
% the bin that contains it instead of a silent NaN, and say which was used.
for q = 1:numel(AP)
  if ~isfinite(AP(q).along), continue; end
  [v, e, how] = at_site(xc, sm, ss, AP(q).along/1e3, 0.5);
  fprintf('radar stack at %-5s (%.2f km along, %.2f km off): %+.2f +/- %.2f mm/m  [%s]\n', ...
    AP(q).name, AP(q).along/1e3, AP(q).off/1e3, v, e, how);
end
he = errbars(axp, xc(okb), sm(okb), ss(okb), PAL.ink);
hs = plot(axp, xc(okb), sm(okb), 'o', 'MarkerSize', 9, 'MarkerFaceColor', PAL.ink, ...
  'MarkerEdgeColor','w', 'LineWidth', 1);
fill(axp, [xg fliplr(xg)]/1e3, [pred_lo fliplr(pred_hi)], PAL.expect, ...
  'FaceAlpha', 0.15, 'EdgeColor','none');
hp = plot(axp, xg/1e3, pred_mm, '--', 'Color', PAL.expect, 'LineWidth', 2.4);
% ApRES at its OWN along-track position (GA04 is the only site with a
% trustworthy tidal admittance - GA01's tracked bed drifts 81 m, so its
% +9.3 is contrast only). Site coordinates recovered 19 Aug 2026.
ha = [];
iga4 = find(strcmp({AP.name},'GA04'), 1);
if ~isempty(iga4) && isfinite(AP(iga4).along)
  xa = AP(iga4).along/1e3;
  plot(axp, [xa xa], APRES_MM + APRES_SE*[-1 1], '-', 'Color', PAL.ink, ...
    'LineWidth', 1.4);
  ha = plot(axp, xa, APRES_MM, 'o', 'MarkerSize', 9, 'MarkerFaceColor','w', ...
    'MarkerEdgeColor', PAL.ink, 'LineWidth', 1.6);
  fprintf('ApRES GA04 drawn at %.2f km along (%.2f km off the line)\n', ...
    xa, AP(iga4).off/1e3);
end
grid(axp,'on'); set(axp, axst{:}); xlim(axp,[0 xmax]);
xlabel(axp,'Along track (km)','Color',PAL.ink);
ylabel(axp,'mm per m of tide (top 100 m)','Color',PAL.ink);
% no panel title and no in-plot ApRES text: the ApRES value is the point
% and error bar at its own along-track position, and the legend names it;
% lettering and captions are added separately
leg_h = [hl hs hp]; leg_l = [lb {'stack (4 lines)','GPS a''''(x) prediction'}];
if ~isempty(ha)
  leg_h(end+1) = ha; leg_l{end+1} = 'ApRES GA04';
end
lg = legend(axp, leg_h, leg_l, ...
  'Location','southoutside','Orientation','horizontal','Interpreter','none');
set(lg,'TextColor',PAL.ink,'Box','off','FontSize',8);

print(h, fullfile(out_dir,'EAGER_2022_admittance_map.png'), '-dpng','-r120');
close(h);
fprintf('Wrote %s\n', fullfile(out_dir,'EAGER_2022_admittance_map.png'));

%% ========================================================================
function S = one_line(pn, net_dir, mp_dir, REF_DEPTH, MAX_BASELINE, BLOCK)
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

P = []; D = []; W = []; along = []; lat = []; lon = []; Nblk = 0;
for q = 1:numel(f)
  tok = regexp(f(q).name, ['^' regexptranslate('escape',pn) '_vvel_(\d+)_(\d+)\.mat$'], ...
    'tokens','once');
  o = load(fullfile(net_dir, f(q).name));
  if isfield(o,'coalign_applied') && ~o.coalign_applied, continue; end
  if max(abs(o.baseline_y)) > MAX_BASELINE, continue; end
  if Nblk == 0
    Nblk = numel(o.S1); along = o.Along_track(:);
    lat = o.Latitude(:); lon = o.Longitude(:);
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
tday = (ptime - min(ptime))/86400; tide = elev - mean(elev);
A = vdef.fitTideAdmittance(N.x, tday, tide);

% GPS a(x): per-block regression of ref_z on its line mean
Nx = numel(L.pass(1).ref_z);
Z = nan(Nx, Np); tid = nan(1,Np);
for k = 1:Np
  z = L.pass(k).ref_z(:);
  if numel(z) == Nx, Z(:,k) = z; tid(k) = mean(z,'omitnan'); end
end
ax_gps = nan(Nblk,1);
for b = 1:min(floor(Nx/BLOCK), Nblk)
  idx = (b-1)*BLOCK+1 : b*BLOCK;
  zb = mean(Z(idx,:),1,'omitnan');
  okz = isfinite(zb) & isfinite(tid);
  if nnz(okz) < 5, continue; end
  pz = polyfit(tid(okz), zb(okz), 1); ax_gps(b) = pz(1);
end
clear L Z;

S = struct('name',pn,'along',along,'lat',lat,'lon',lon, ...
  'adm', 1e3*REF_DEPTH*A.admittance(:), ...
  'adm_std', 1e3*REF_DEPTH*A.admittance_std(:), ...
  'ax_gps', ax_gps, 'calib', true);
end

%% ========================================================================
function [v, e, how] = at_site(xc, sm, ss, x0, w)
% Stacked value at along-track position x0, and a string saying where it
% came from. Interpolation is only defined between the first and last bin
% centre that carry a stack; a site outside that range - the far-end sites
% especially - falls back to the bin that CONTAINS it rather than to NaN.
ok = isfinite(sm) & isfinite(ss);
v = NaN; e = NaN; how = 'no stacked bin at this position';
if ~any(ok), return; end
xo = xc(ok); so = sm(ok); eo = ss(ok);
if numel(xo) >= 2 && x0 >= min(xo) && x0 <= max(xo)
  v = interp1(xo, so, x0, 'linear');
  e = interp1(xo, eo, x0, 'linear');
  how = 'interpolated between bin centres';
  return;
end
[dmin, k] = min(abs(xo - x0));
if dmin <= w/2
  v = so(k); e = eo(k);
  how = sprintf('%.2f-%.2f km bin', xo(k)-w/2, xo(k)+w/2);
end
end

%% ========================================================================
function hh = errbars(ax, x, y, s, col)
hh = [];
for k = 1:numel(x)
  if ~isfinite(s(k)), continue; end
  hh = plot(ax, [x(k) x(k)], [y(k)-s(k) y(k)+s(k)], '-', 'Color', col, 'LineWidth', 1.4);
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
function AP = read_apres_xy(fn)
%READ_APRES_XY Site name and EPSG:3031 metres from the whitespace table
%   `x  y  name`. Returns an empty struct (with a note) when the file is
%   absent, so the figure still builds without the site coordinates.
AP = struct('name',{},'x',{},'y',{},'along',{},'off',{});
if ~exist(fn,'file')
  fprintf('ApRES site coordinates not found (%s) - sites not drawn.\n', fn);
  return;
end
try
  fid = fopen(fn,'r');
  C = textscan(fid, '%f %f %s');
  fclose(fid);
  for k = 1:numel(C{3})
    AP(end+1) = struct('name', C{3}{k}, 'x', C{1}(k), 'y', C{2}(k), ...
      'along', NaN, 'off', NaN); %#ok<AGROW>
  end
  fprintf('ApRES sites read: %s\n', strjoin({AP.name}, ', '));
catch ME
  fprintf('Could not read %s (%s) - sites not drawn.\n', fn, ME.message);
end
end

%% ========================================================================
function [along_m, off_m] = apres_along(R, ps, X, Y)
%APRES_ALONG Along-track position of an EPSG:3031 point (metres).
%   The point is PROJECTED onto the block track of the calibrated lines,
%   not snapped to the nearest block: blocks are ~500 m apart, so snapping
%   would quantise the position to a quarter of the flexure scale. Returns
%   the along-track coordinate in the radar's own frame plus the
%   perpendicular distance, which matters here because the sites sit a few
%   hundred metres OFF the line rather than on it.
along_m = NaN; off_m = NaN; best = inf;
for i = 1:numel(R)
  if ~R(i).calib, continue; end
  [bx, by] = ps_km(ps, R(i).lon, R(i).lat);
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
      along_m = R(i).along(k) + t*(R(i).along(kk) - R(i).along(k));
    end
  end
end
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
