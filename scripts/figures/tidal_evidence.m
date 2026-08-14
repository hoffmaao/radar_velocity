%TIDAL_EVIDENCE What the tidal analysis actually established, in one figure.
%
%   This is the retraction figure. The project's headline result was a
%   tidal flexure hinge - vertical column strain correlating with the tide
%   and reversing sign along track. It did not survive correct
%   coalignment. Prose in the README says so; this shows it.
%
%   FOUR PANELS, in the order the argument runs:
%
%   (a) THE CLAIM. Tidal response along track from the products built with
%       a SCALAR coalignment (CSARP_vvel_v2). All five lines swing across
%       zero: an apparent hinge, at a consistent place.
%
%   (b) THE RETRACTION. The same thing from the per-column coalignment
%       (_v3). Same axes, same products, same ice. The swing is gone and
%       everything sits inside the method's own floor.
%
%   (c) WHY IT WAS NEVER REAL. A scalar correction removes a line MEAN, so
%       what survives goes as (a(x) - 1), where a(x) is the local surface
%       tidal admittance from GPS normalised to line-mean 1. That residual
%       must change sign exactly where a(x) = 1 - a position fixed by the
%       arbitrary normalisation of the survey, not by the ice. Panel (c)
%       plots a(x) with that crossing marked and the panel-(a) hinges on
%       top. Four of the five land within about 0.5 km of the crossing.
%       The exception is EAGER_2022 at 1.50 km, which is also the only
%       wholly uncalibrated product (coregistration_time_shift AND
%       equalization both absent), so its residual is the largest and
%       least well described by the simple (a(x) - 1) form. Real flexure
%       would instead change sign at an INFLECTION of a(x), and a(x) is
%       concave-down throughout the surveyed window.
%
%   (d) THE NUMBER BEHIND IT. Correlation between the measured tidal
%       response and the residual-misalignment predictor, per product,
%       before and after. It collapses from -0.84..-0.97 (every product
%       p < 0.005) to scattered and insignificant.
%
%   WHAT SURVIVES is drawn on (b): the ApRES measurement, +3.8 mm per
%   metre of tide at 100 m from an instrument sitting on the ice. The
%   tidal strain is real. It is simply ~2x below what this method can
%   currently resolve, which is the honest conclusion.
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../tidal_evidence.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef

root     = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
mp_dir   = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
out_dir  = '/kucresis/scratch/hoffmana_sta/vvel/figures';
apres_fn = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
  'diagnostics', 'apres_GA04_tide_profile.csv');

PASS_NAMES = {'EAGER_2022','EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};
REF_DEPTH    = 100;
MAX_BASELINE = 10;
BLOCK        = 200;
ALPHA        = 1.03;    % fraction of the erroneous compensation that survives
c_light      = 299792458;

if ~exist(out_dir,'dir'), mkdir(out_dir); end

PAL.cat = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; ...
           0.929 0.631 0.000; 0.910 0.482 0.643];
PAL.cat_mk  = {'o','s','^','d','v'};
PAL.ink      = [0.20 0.20 0.20];
PAL.ink_soft = [0.45 0.45 0.45];
PAL.band     = [0.90 0.90 0.88];
PAL.expect   = [0.35 0.35 0.35];

%% Measure both generations
V2 = analyse_all(fullfile(root,'CSARP_vvel_v2'), mp_dir, PASS_NAMES, ...
  REF_DEPTH, MAX_BASELINE, BLOCK, ALPHA, c_light);
V3 = analyse_all(fullfile(root,'CSARP_vvel_v3'), mp_dir, PASS_NAMES, ...
  REF_DEPTH, MAX_BASELINE, BLOCK, ALPHA, c_light);
assert(~isempty(V2) && ~isempty(V3), 'need both _v2 and _v3 products');

% same-leg floor from v3, in mm
FLOOR_MM = same_leg_floor(V3, 'EAGER_2022', 'EAGER_2022_GL1');

apres_mm = NaN;
if exist(apres_fn,'file')
  A = importdata(apres_fn, ',', 1);
  apres_mm = interp1(A.data(:,1), A.data(:,2), REF_DEPTH, 'linear', NaN);
end

%% Report
fprintf('\n%-16s | %-28s | %-28s\n','product','SCALAR coalign (v2)','PER-COLUMN coalign (v3)');
fprintf('%-16s | %10s %8s %8s | %10s %8s %8s\n','', ...
  'hinge[km]','corr','p','hinge[km]','corr','p');
for i = 1:numel(V3)
  fprintf('%-16s | %10s %8.2f %8.3f | %10s %8.2f %8.3f\n', V3(i).name, ...
    numstr(V2(i).hinge_km,'%.2f'), V2(i).r_art, V2(i).p_art, ...
    numstr(V3(i).hinge_km,'%.2f'), V3(i).r_art, V3(i).p_art);
end
fprintf('\nfloor %.1f mm per m of tide; ApRES measured %+.1f mm at %d m\n', ...
  FLOOR_MM, apres_mm, REF_DEPTH);

%% Figure
h = figure('Visible','off','Position',[100 100 1120 940],'Color','w');
set(0,'CurrentFigure',h);
axst = {'GridAlpha',0.15,'XColor',PAL.ink_soft,'YColor',PAL.ink_soft,'Box','off'};

xmax = 0;
for i = 1:numel(V3), xmax = max(xmax, max(V3(i).along)/1e3); end
% products have different block counts, so collect rather than concatenate
mx = FLOOR_MM;
for i = 1:numel(V2), mx = max(mx, max(abs(V2(i).adm_mm))); end
for i = 1:numel(V3), mx = max(mx, max(abs(V3(i).adm_mm))); end
ylim_mm = [-1 1] * 1.15 * mx;

% (a) and (b): before and after, on shared axes
for pnl = 1:2
  if pnl == 1, S = V2; ttl = '(a) With a SCALAR coalignment: an apparent hinge';
  else,        S = V3; ttl = '(b) With per-column coalignment: it is gone';
  end
  ax = axes('parent',h,'Position',[0.075 0.565-(pnl-1)*0.0 0.40 0.36]);
  if pnl == 2, set(ax,'Position',[0.555 0.565 0.40 0.36]); end
  set(0,'CurrentFigure',h); hold(ax,'on');
  fill(ax, [0 xmax xmax 0], [-1 -1 1 1]*FLOOR_MM, PAL.band, 'EdgeColor','none');
  plot(ax, [0 xmax], [0 0], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
  if pnl == 2 && isfinite(apres_mm)
    plot(ax, [0 xmax], [1 1]*apres_mm, '-', 'Color', PAL.expect, 'LineWidth', 3);
  end
  hl = [];
  for i = 1:numel(S)
    hl(end+1) = plot(ax, S(i).along/1e3, S(i).adm_mm, '-', ...
      'Color', PAL.cat(i,:), 'LineWidth', 1.8); %#ok<AGROW>
    plot(ax, S(i).along/1e3, S(i).adm_mm, PAL.cat_mk{i}, 'MarkerSize', 6, ...
      'MarkerFaceColor', PAL.cat(i,:), 'MarkerEdgeColor','w','LineWidth',0.8);
    if isfinite(S(i).hinge_km)
      plot(ax, S(i).hinge_km, 0, 'p', 'MarkerSize', 15, 'MarkerFaceColor','w', ...
        'MarkerEdgeColor', PAL.ink, 'LineWidth', 1.2);
    end
  end
  grid(ax,'on'); set(ax, axst{:}); xlim(ax,[0 xmax]); ylim(ax, ylim_mm);
  xlabel(ax,'Along track (km)','Color',PAL.ink);
  if pnl == 1
    ylabel(ax, sprintf('Tidal response, top %d m (mm per m of tide)', REF_DEPTH), ...
      'Color', PAL.ink);
  end
  title(ax, ttl, 'Color', PAL.ink);
  if pnl == 2
    lg = legend(ax, hl, {S.name}, 'Location','southeast','Interpreter','none');
    set(lg,'TextColor',PAL.ink,'Box','off','FontSize',8);
    text(ax, 0.98*xmax, apres_mm, 'ApRES measured ', 'Color', PAL.expect, ...
      'FontSize', 9, 'HorizontalAlignment','right','VerticalAlignment','bottom');
  else
    text(ax, 0.02*xmax, -FLOOR_MM, ' method floor', 'Color', PAL.ink_soft, ...
      'FontSize', 9, 'VerticalAlignment','top');
  end
end

% (c) the mechanism
ax3 = axes('parent',h,'Position',[0.075 0.085 0.40 0.36]);
set(0,'CurrentFigure',h); hold(ax3,'on');
plot(ax3, [0 xmax], [1 1], '--', 'Color', PAL.expect, 'LineWidth', 2);
for i = 1:numel(V3)
  plot(ax3, V3(i).along/1e3, V3(i).ax_gps, '-', 'Color', PAL.cat(i,:), 'LineWidth', 1.8);
  plot(ax3, V3(i).along/1e3, V3(i).ax_gps, PAL.cat_mk{i}, 'MarkerSize', 5, ...
    'MarkerFaceColor', PAL.cat(i,:), 'MarkerEdgeColor','w','LineWidth',0.8);
end
yl = ylim(ax3);
for i = 1:numel(V2)
  if isfinite(V2(i).hinge_km)
    plot(ax3, [1 1]*V2(i).hinge_km, yl, ':', 'Color', PAL.ink, 'LineWidth', 1.2);
  end
end
ylim(ax3, yl);
grid(ax3,'on'); set(ax3, axst{:}); xlim(ax3,[0 xmax]);
xlabel(ax3,'Along track (km)','Color',PAL.ink);
ylabel(ax3,'GPS surface tidal admittance a(x)','Color',PAL.ink);
title(ax3,'(c) The scalar residual must flip sign where a(x) = 1','Color',PAL.ink);
text(ax3, 0.02*xmax, 1, ' a(x) = 1', 'Color', PAL.expect, 'FontSize', 9, ...
  'VerticalAlignment','bottom');
text(ax3, 0.98*xmax, yl(2), 'dotted: panel (a) hinges ', 'Color', PAL.ink, ...
  'FontSize', 9, 'HorizontalAlignment','right','VerticalAlignment','top');

% (d) the correlation collapse
ax4 = axes('parent',h,'Position',[0.555 0.085 0.40 0.36]);
set(0,'CurrentFigure',h); hold(ax4,'on');
plot(ax4, [0.5 2.5], [0 0], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
for i = 1:numel(V3)
  plot(ax4, [1 2], [V2(i).r_art V3(i).r_art], '-', 'Color', PAL.cat(i,:), 'LineWidth', 1.8);
  plot(ax4, 1, V2(i).r_art, PAL.cat_mk{i}, 'MarkerSize', 9, ...
    'MarkerFaceColor', PAL.cat(i,:), 'MarkerEdgeColor','w','LineWidth',1);
  plot(ax4, 2, V3(i).r_art, PAL.cat_mk{i}, 'MarkerSize', 9, ...
    'MarkerFaceColor', PAL.cat(i,:), 'MarkerEdgeColor','w','LineWidth',1);
end
grid(ax4,'on'); set(ax4, axst{:});
xlim(ax4,[0.6 2.4]); ylim(ax4,[-1.05 1.05]);
set(ax4,'XTick',[1 2],'XTickLabel',{'scalar (v2)','per-column (v3)'});
ylabel(ax4,'Correlation with the artefact predictor','Color',PAL.ink);
title(ax4,'(d) The response stops tracking the artefact','Color',PAL.ink);

print(h, fullfile(out_dir,'EAGER_2022_tidal_evidence.png'), '-dpng', '-r120');
close(h);
fprintf('\nWrote %s\n', fullfile(out_dir,'EAGER_2022_tidal_evidence.png'));

%% ========================================================================
function R = analyse_all(vdir, mp_dir, names, REF_DEPTH, MAX_BASELINE, BLOCK, ALPHA, c_light)
R = [];
for n = 1:numel(names)
  pn = names{n};
  f = dir(fullfile(vdir, [pn '_vvel_*.mat']));
  keep = ~cellfun('isempty', regexp({f.name}, ...
    ['^' regexptranslate('escape',pn) '_vvel_\d+_\d+\.mat$'], 'once'));
  f = f(keep);
  if isempty(f), continue; end

  L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass');
  Np = numel(L.pass); elev = nan(1,Np);
  for k = 1:Np, elev(k) = mean(L.pass(k).elev,'omitnan'); end
  Zref = L.pass; % keep for a(x)
  spy = 365.25*86400; np = numel(f);
  strain = []; resid = []; tide = nan(1,np); tday = nan(1,np); mb = nan(1,np);
  along = []; ref0 = NaN;
  for i = 1:np
    o = load(fullfile(vdir, f(i).name));
    if isnan(ref0), ref0 = o.pass_idx_ref; end
    if isempty(strain)
      Nblk = numel(o.S1);
      strain = nan(Nblk, np); resid = nan(Nblk, np);
      along = o.Along_track(:);
    end
    for b = 1:size(strain,1)
      d = o.depth_blk(:,b); ok = isfinite(d) & isfinite(o.dh_blk(:,b));
      if ~any(ok) || max(d(ok)) < REF_DEPTH, continue; end
      strain(b,i) = interp1(d(ok), o.dh_blk(ok,b), REF_DEPTH,'linear',NaN)/REF_DEPTH;
    end
    % artefact predictor: along-track residual misalignment, line mean removed
    zs = Zref(o.pass_idx_sec).ref_z(:); zr = Zref(ref0).ref_z(:);
    dd = ALPHA * -(zs - zr)/(c_light/2); dd = dd - mean(dd,'omitnan');
    nb = floor(numel(dd)/BLOCK);
    for b = 1:min(nb, size(resid,1))
      resid(b,i) = mean(dd((b-1)*BLOCK+1 : b*BLOCK),'omitnan');
    end
    tide(i) = elev(o.pass_idx_sec) - elev(ref0);
    tday(i) = mean(o.GPS_time + o.delta_t_blk*spy,'omitnan')/86400;
    mb(i)   = max(abs(o.baseline_y));
  end
  use = ~(isfinite(mb) & mb > MAX_BASELINE);
  if nnz(use) < 6, clear L; continue; end
  strain = strain(:,use); resid = resid(:,use);
  tide = tide(use); tday = tday(use) - min(tday(use));

  A = vdef.fitTideAdmittance(strain, tday, tide);
  adm_mm = 1e3 * A.admittance(:) * REF_DEPTH;
  rp = A.r_partial(:);

  % per-block artefact predictor slope against tide
  Nblk = size(resid,1); g = nan(Nblk,1);
  for b = 1:Nblk
    v = resid(b,:); ok = isfinite(v) & isfinite(tide);
    if nnz(ok) < 5, continue; end
    p = polyfit(tide(ok), v(ok), 1); g(b) = p(1)*1e9;
  end
  ok = isfinite(adm_mm) & isfinite(g);
  r_art = NaN; p_art = NaN;
  if nnz(ok) >= 5
    m = corrcoef(g(ok), adm_mm(ok)); r_art = m(1,2);
    nn = nnz(ok); tt = r_art*sqrt((nn-2)/max(1-r_art^2,eps));
    p_art = betainc((nn-2)/((nn-2)+tt^2), (nn-2)/2, 0.5);
  end

  % GPS surface tidal admittance a(x), normalised to line mean 1
  Nx = numel(Zref(ref0).ref_z);
  Z = nan(Nx, Np); tid = nan(1,Np);
  for k = 1:Np
    z = Zref(k).ref_z(:);
    if numel(z) == Nx, Z(:,k) = z; tid(k) = mean(z,'omitnan'); end
  end
  nb = floor(Nx/BLOCK); ax_gps = nan(Nblk,1);
  for b = 1:min(nb,Nblk)
    idx = (b-1)*BLOCK+1 : b*BLOCK;
    zb = mean(Z(idx,:),1,'omitnan');
    okz = isfinite(zb) & isfinite(tid);
    if nnz(okz) < 5, continue; end
    pz = polyfit(tid(okz), zb(okz), 1); ax_gps(b) = pz(1);
  end
  clear L Zref;

  hinge_km = change_point(rp, along);
  S = struct('name',pn,'along',along,'adm_mm',adm_mm,'rp',rp,'g',g, ...
    'r_art',r_art,'p_art',p_art,'ax_gps',ax_gps,'hinge_km',hinge_km);
  if isempty(R), R = S; else, R(end+1) = S; end %#ok<AGROW>
end
end

%% ========================================================================
function xc = change_point(rr, along)
% Strongest positive-to-negative step, the same detector the tidal analysis
% uses. Returns NaN when no step clears MIN_CONTRAST, so a line with no
% hinge reports none rather than inventing one.
MIN_CONTRAST = 0.5; MIN_SIDE = 3;
xc = NaN;
ok = find(isfinite(rr)); n = numel(ok);
if n < 2*MIN_SIDE, return; end
rv = rr(ok); xx = along(ok);
best = -inf; kb = [];
for k = MIN_SIDE:(n-MIN_SIDE)
  d = mean(rv(1:k)) - mean(rv(k+1:end));
  if d > best, best = d; kb = k; end
end
if isempty(kb) || best < MIN_CONTRAST, return; end
xc = 0.5*(xx(kb)+xx(kb+1))/1e3;
end

%% ========================================================================
function fl = same_leg_floor(R, nameA, nameB)
fl = NaN;
ia = find(strcmp({R.name}, nameA),1); ib = find(strcmp({R.name}, nameB),1);
if isempty(ia) || isempty(ib), return; end
A = R(ia); B = R(ib); TOL = 100; d = [];
for k = 1:numel(A.along)
  [gap,m] = min(abs(B.along - A.along(k)));
  if gap <= TOL && isfinite(A.adm_mm(k)) && isfinite(B.adm_mm(m))
    d(end+1) = A.adm_mm(k) - B.adm_mm(m); %#ok<AGROW>
  end
end
if ~isempty(d), fl = sqrt(mean(d.^2)); end
end

%% ========================================================================
function s = numstr(v, fmt)
if isfinite(v), s = sprintf(fmt, v); else, s = 'none'; end
end
