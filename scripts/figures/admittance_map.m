%ADMITTANCE_MAP Map view of the tidal response, against the GPS-predicted pattern.
%
%   Two panels:
%
%   (a) MAP. Per-block tidal response (mm of column change per metre of
%       tide, top 100 m, network inversion) for all five lines, on the
%       local tangent plane with the MEaSUREs grounding line. Diverging
%       colour about zero. This is where a spatial pattern would show.
%
%   (b) THE PATTERN TEST. All calibrated lines' per-block responses
%       against along-track position, with their inverse-variance stack -
%       and over it, the response PREDICTED FROM GPS ALONE: plate bending
%       says the tidal strain is set by the curvature of the deflection
%       profile, eps_zz = nu/(1-nu) * (H/2) * a''(x) * (1 - z/H), and
%       a(x) is measured by the GPS with no radar involved. If the ApRES
%       and radar magnitudes are right, the stacked radar profile should
%       sit on the GPS-curvature prediction. The ApRES value is drawn as
%       a band, not a point on the map - its site position was never
%       recorded (GPS was off for the whole deployment).
%
%   Run on the server (needs CSARP_vvel_net):
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../admittance_map.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef

root    = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
mp_dir  = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
net_dir = fullfile(root,'CSARP_vvel_net');
gis_dir = '/kucresis/scratch/hoffmana_sta/vvel/gis';
out_dir = '/kucresis/scratch/hoffmana_sta/vvel/figures';

PASS_NAMES = {'EAGER_2022','EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};
CALIB = [false true true true true];
REF_DEPTH = 100; MAX_BASELINE = 10; BLOCK = 200;
H_ICE = 300; NU = 0.33;
APRES_MM = -1.24; APRES_SE = 0.04;    % rate method, apres_rate_check.py

PAL.cat = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; ...
           0.929 0.631 0.000; 0.910 0.482 0.643];
PAL.cat_mk = {'o','s','^','d','v'};
PAL.ink = [0.20 0.20 0.20]; PAL.ink_soft = [0.45 0.45 0.45];
PAL.div_neg = [0.698 0.094 0.169]; PAL.div_mid = [0.941 0.937 0.925];
PAL.div_pos = [0.165 0.471 0.839]; PAL.band = [0.90 0.90 0.88];
PAL.expect = [0.35 0.35 0.35];

%% Per-block network admittance + GPS a(x) for every line
R = [];
for n = 1:numel(PASS_NAMES)
  S = one_line(PASS_NAMES{n}, net_dir, mp_dir, REF_DEPTH, MAX_BASELINE, BLOCK);
  if isempty(S), continue; end
  S.calib = CALIB(n);
  if isempty(R), R = S; else, R(end+1) = S; end %#ok<AGROW>
  fprintf('%-16s %2d blocks, adm %.2f..%.2f mm/m\n', S.name, numel(S.adm), ...
    min(S.adm), max(S.adm));
end
assert(numel(R) >= 4, 'need the lines');

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
axst = {'GridAlpha',0.15,'XColor',PAL.ink_soft,'YColor',PAL.ink_soft,'Box','off'};

% local tangent frame
lat0 = mean(cellfun(@(v) mean(v,'omitnan'), {R.lat}));
lon0 = mean(cellfun(@(v) mean(v,'omitnan'), {R.lon}));
xy = @(lon,lat) deal((lon-lon0)*111320*cosd(lat0)/1e3, (lat-lat0)*110540/1e3);

% (a) map
axm = axes('parent',h,'Position',[0.06 0.12 0.40 0.78]);
set(0,'CurrentFigure',h); hold(axm,'on');
CLIM = 4;   % mm/m colour range; values beyond are clipped to the ends
ndiv = 256; half = ndiv/2;
dmap = [interp1([0 1],[PAL.div_neg; PAL.div_mid], linspace(0,1,half)); ...
        interp1([0 1],[PAL.div_mid; PAL.div_pos], linspace(0,1,ndiv-half))];
for i = 1:numel(R)
  [xk, yk] = xy(R(i).lon, R(i).lat);
  plot(axm, xk, yk, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
  for b = 1:numel(R(i).adm)
    if ~isfinite(R(i).adm(b)), continue; end
    ci = max(1, min(ndiv, round((R(i).adm(b)+CLIM)/(2*CLIM)*(ndiv-1))+1));
    mk = PAL.cat_mk{i};
    ec = PAL.ink_soft; if ~R(i).calib, ec = [0.75 0.75 0.75]; end
    plot(axm, xk(b), yk(b), mk, 'MarkerSize', 10, ...
      'MarkerFaceColor', dmap(ci,:), 'MarkerEdgeColor', ec, 'LineWidth', 0.6);
  end
end
gl_drawn = overlay_gl(axm, gis_dir, xy);
grid(axm,'on'); set(axm, axst{:}); axis(axm,'equal');
xlabel(axm, sprintf('East of %.4f deg (km)', lon0),'Color',PAL.ink);
ylabel(axm, sprintf('North of %.4f deg (km)', lat0),'Color',PAL.ink);
ttl = '(a) Tidal response in map view';
if gl_drawn, ttl = [ttl ' (black: grounding line)']; end
title(axm, ttl, 'Color', PAL.ink);
colormap(axm, dmap); caxis(axm, [-CLIM CLIM]);
cb = colorbar(axm);
set(get(cb,'ylabel'),'string','mm per m of tide (top 100 m)','Color',PAL.ink);
set(cb,'XColor',PAL.ink_soft,'YColor',PAL.ink_soft);

% (b) pattern test
axp = axes('parent',h,'Position',[0.58 0.12 0.385 0.78]);
set(0,'CurrentFigure',h); hold(axp,'on');
xmax = 0;
for i = 1:numel(R), xmax = max(xmax, max(R(i).along)/1e3); end
% ApRES band (site position unknown, so a band across the axis)
fill(axp, [0 xmax xmax 0], APRES_MM + APRES_SE*[-1 -1 1 1], PAL.band, ...
  'EdgeColor','none');
plot(axp, [0 xmax], [0 0], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
hl = []; lb = {};
for i = 1:numel(R)
  if ~R(i).calib, continue; end
  hh = plot(axp, R(i).along/1e3, R(i).adm, PAL.cat_mk{i}, 'MarkerSize', 5, ...
    'MarkerFaceColor', PAL.cat(i,:), 'MarkerEdgeColor','w', 'LineWidth', 0.6);
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
he = errbars(axp, xc(okb), sm(okb), ss(okb), PAL.ink);
hs = plot(axp, xc(okb), sm(okb), 'o', 'MarkerSize', 9, 'MarkerFaceColor', PAL.ink, ...
  'MarkerEdgeColor','w', 'LineWidth', 1);
fill(axp, [xg fliplr(xg)]/1e3, [pred_lo fliplr(pred_hi)], PAL.expect, ...
  'FaceAlpha', 0.15, 'EdgeColor','none');
hp = plot(axp, xg/1e3, pred_mm, '--', 'Color', PAL.expect, 'LineWidth', 2.4);
grid(axp,'on'); set(axp, axst{:}); xlim(axp,[0 xmax]);
xlabel(axp,'Along track (km)','Color',PAL.ink);
ylabel(axp,'mm per m of tide (top 100 m)','Color',PAL.ink);
title(axp,'(b) Stacked response vs the GPS-curvature prediction','Color',PAL.ink);
lg = legend(axp, [hl hs hp], [lb {'stack (4 lines)','GPS a''''(x) prediction'}], ...
  'Location','southoutside','Orientation','horizontal','Interpreter','none');
set(lg,'TextColor',PAL.ink,'Box','off','FontSize',8);
text(axp, 0.02*xmax, APRES_MM, ' ApRES (site position unrecorded)', ...
  'Color', PAL.ink_soft, 'FontSize', 8.5, 'VerticalAlignment','bottom');

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

N = vdef.invertNetwork(P, D, struct('n_sigma',3,'weights',W));
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
function hh = errbars(ax, x, y, s, col)
hh = [];
for k = 1:numel(x)
  if ~isfinite(s(k)), continue; end
  hh = plot(ax, [x(k) x(k)], [y(k)-s(k) y(k)+s(k)], '-', 'Color', col, 'LineWidth', 1.4);
end
end

%% ========================================================================
function drawn = overlay_gl(ax, gis_dir, xy)
drawn = false;
fn = fullfile(gis_dir,'GroundingLine_Antarctica_v02.shp');
if ~exist(fn,'file'), return; end
try
  S = shaperead(fn);
  X = []; Y = [];
  for k = 1:numel(S)
    X = [X; S(k).X(:); NaN]; Y = [Y; S(k).Y(:); NaN]; %#ok<AGROW>
  end
  [glat, glon] = projinv(projcrs(3031), X, Y);
  [gx, gy] = xy(glon, glat);
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
