function OUT = elasticity_results(opts)
%ELASTICITY_RESULTS Figure: effective elasticity of the Windless Bight
%   grounding zone, as results rather than diagnostics.
%
%   Three panels:
%     (a) THE OBSERVATION AND THE FITS. Surface tidal admittance per line
%         with the fitted beams and the fitted clamp positions. Lines the
%         gate excludes - E_hi unbounded, OR the beam misfitting the
%         line's own shape data (chi-squared >= 3; GL2 on current data,
%         at 3.7 on its known ref_z anomaly, the strain having closed its
%         upper bound) - are drawn in grey: shown, not averaged.
%     (b) THE THICKNESS THAT WENT IN. The tracked bed h(x) per line, the
%         BedMachine pseudo-layer it replaced (grey dashed, one line), and
%         the ApRES bed depths at the two stable sites, converted through
%         the SAME firn column so like is compared with like. This panel is
%         why E* moved from 5.9 to 3.8 GPa: BedMachine is ~45 m thin here.
%     (c) THE MODULUS IN CONTEXT. Forest of E* per line and the
%         three-line mean, the same lines under BedMachine thickness (the
%         data-product bias, kept visible), and the published estimates:
%         Vaughan (1995), Sayag and Worster (2013) - whose two values are
%         the same data under two bed models, the cautionary tale for this
%         parameter - Elgart, Minchew and Meyer (2025), and laboratory ice.
%         The grey band is Elgart's site-to-site range on the Ross shelf.
%
%   Fits come from scripts/diagnostics/elastic_modulus.m (pass a saved OUT
%   via opts.fit to redraw without refitting), so the inversion has one
%   implementation and this figure cannot drift from it.
%
%   opts, all optional:
%     .fit      OUT struct from elastic_modulus - skips the refit
%     .out_dir  where the png goes
%     .mp_dir, .products, .block   passed through when refitting
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "addpath('<code>/scripts/figures'); elasticity_results"

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root);                                   % +vdef
addpath(fullfile(root,'scripts','diagnostics')); % the inversion driver

if ~isfield(opts,'out_dir') || isempty(opts.out_dir)
  opts.out_dir = '/kucresis/scratch/hoffmana_sta/vvel/figures_flexure';
end
if ~isfield(opts,'fit'), opts.fit = []; end

if isempty(opts.fit)
  dopts = struct();
  for f = {'mp_dir','products','block'}
    if isfield(opts, f{1}), dopts.(f{1}) = opts.(f{1}); end
  end
  D = elastic_modulus(dopts);
else
  D = opts.fit;
end
OUT = D;
L = D.lines; R = D.fits;
% Case rows in the driver's H_CASES: 1 = tracked bed with the englacial
% strain joined, 2 = the same thickness a(x)-only, 3 = BedMachine.
CASE = 1; AXCASE = 2; BMCASE = 3;
nL = numel(L);

% Constrained = the data bound the modulus from above AND the beam
% actually describes the line's flexure data (shape chi-squared below 3).
% Both clauses are needed and the second became load-bearing when the
% strain joined the fit: the strain amplitude closes GL2's upper bound, so
% "unbounded" alone would silently promote into the mean a line whose a(x)
% carries a physically impossible negative block and misfits at 3.7.
% Criterion, not name, so the figure tracks the data across rebuilds.
con = false(1,nL);
for i = 1:nL
  Ri = R{CASE,i};
  if isempty(Ri), continue; end
  x2a = Ri.chi2red;
  if Ri.has_strain && isfinite(Ri.chi2_shape), x2a = Ri.chi2_shape; end
  con(i) = isfinite(Ri.E_hi) && isfinite(x2a) && x2a < 3;
end

% ApRES bed depths at the two stable sites, converted through the SAME
% firn column as the radar traveltimes (the sites' own constant-n depths
% are 271.8 and 250.5 m; er_ice = 3.18, no firn term). GA01/GA05 excluded
% for drifting bed picks, as everywhere else in this project.
SITES = struct( ...
  'name', {'GA10','GA04'}, ...
  'lat',  {-77.666664390990462, -77.652087106641048}, ...
  'lon',  { 168.145574208061106,  168.166278345287481}, ...
  'h',    { 286.4, 265.0});

cols = [0.11 0.42 0.69; 0.89 0.47 0.10; 0.20 0.60 0.25; 0.75 0.20 0.30];
grey = [0.62 0.62 0.62];
linecol = @(i) cols(mod(i-1,size(cols,1))+1,:);

hf = figure('Visible','off','Position',[100 100 1280 420]);
set(0,'CurrentFigure',hf);

%% ---- (a) observation and fits
ax = axes('parent',hf,'Position',[0.055 0.13 0.27 0.78]);
hold(ax,'on'); xlo = inf; xhi = -inf; hleg = []; lleg = {};
for i = 1:nL
  cc = linecol(i); if ~con(i), cc = grey; end
  ok = isfinite(L(i).a);
  for b = find(ok(:).')
    if ~isfinite(L(i).a_std(b)), continue; end
    plot(ax, [1 1]*L(i).x_sea(b)/1e3, L(i).a(b)+[-1 1]*L(i).a_std(b), ...
      '-', 'Color', cc, 'LineWidth', 0.9);
  end
  hleg(end+1) = plot(ax, L(i).x_sea(ok)/1e3, L(i).a(ok), 'o', ...
    'Color', cc, 'MarkerSize', 3.5, 'MarkerFaceColor', cc); %#ok<AGROW>
  lleg{end+1} = shortname(L(i).name); %#ok<AGROW>
  Ri = R{CASE,i};
  if isempty(Ri), continue; end
  plot(ax, Ri.x_grid/1e3, Ri.w_grid, '-', 'Color', cc, 'LineWidth', 1.3);
  plot(ax, Ri.x0/1e3, 0, '^', 'Color', cc, 'MarkerSize', 6, 'MarkerFaceColor', cc);
  xlo = min([xlo; Ri.x0/1e3 - 0.35]);
  xhi = max([xhi; max(L(i).x_sea(ok))/1e3 + 0.35]);
end
if isfinite(xlo) && xhi > xlo, xlim(ax,[xlo xhi]); end
grid(ax,'on'); box(ax,'on');
xlabel(ax,'Seaward distance (km)'); ylabel(ax,'a(x) (line-mean units)');
title(ax,'(a)  Tidal flexure and fitted beams');
legend(ax, hleg, lleg, 'Location','SouthEast');

%% ---- (b) thickness
ax = axes('parent',hf,'Position',[0.385 0.13 0.25 0.78]);
hold(ax,'on');
bm_drawn = false;
for i = 1:nL
  if ~con(i), continue; end
  hs = R{CASE,i}.h_spec;
  if isnumeric(hs) && size(hs,2) == 2
    plot(ax, hs(:,1)/1e3, hs(:,2), 'o-', 'Color', linecol(i), ...
      'MarkerSize', 3.5, 'LineWidth', 1.3, 'MarkerFaceColor', linecol(i));
  elseif isnumeric(hs) && isscalar(hs)
    % The driver fell back to its scalar segment median for this line; a
    % constrained line must not vanish from the thickness panel while it
    % stands in the E* forest, so the constant it actually used is drawn.
    fprintf(['panel (b): %s used the scalar fallback thickness %.0f m - ' ...
             'drawn as a flat dashed line\n'], shortname(L(i).name), hs);
    if isfinite(xlo) && xhi > xlo
      plot(ax, [max(xlo,-0.5) xhi], [hs hs], '--', 'Color', linecol(i), ...
        'LineWidth', 1.1);
    end
  end
  if ~bm_drawn && ~isempty(R{BMCASE,i})
    hb = R{BMCASE,i}.h_spec;
    if isnumeric(hb) && size(hb,2) == 2
      plot(ax, hb(:,1)/1e3, hb(:,2), '--', 'Color', grey, 'LineWidth', 1.3);
      bm_drawn = true;
    end
  end
end
% Sites, projected onto each constrained line's own blocks
for q = 1:numel(SITES)
  for i = 1:nL
    if ~con(i), continue; end
    xs = site_xsea(L(i), SITES(q).lat, SITES(q).lon);
    if ~isfinite(xs), continue; end
    plot(ax, xs/1e3, SITES(q).h, 's', 'Color', 'k', ...
      'MarkerSize', 7, 'LineWidth', 1.4);
  end
end
grid(ax,'on'); box(ax,'on');
if isfinite(xlo) && xhi > xlo, xlim(ax,[max(xlo,-0.5) xhi]); end
xlabel(ax,'Seaward distance (km)'); ylabel(ax,'Ice thickness (m)');
title(ax,'(b)  Tracked bed, BedMachine, ApRES');

%% ---- (c) the modulus in context
ax = axes('parent',hf,'Position',[0.70 0.13 0.285 0.78]);
hold(ax,'on');
XL = [0.4 22];

% Elgart's site-to-site range on the Ross shelf, behind everything
patch(ax, [0.6 9 9 0.6], [0 0 100 100], [0.92 0.92 0.92], 'EdgeColor','none');

rows = {};   % {label, y} accumulated as they are drawn
y = 0;

% per-line, tracked bed
for i = 1:nL
  if ~con(i), continue; end
  y = y + 1; Ri = R{CASE,i};
  plot(ax, [Ri.E_lo Ri.E_hi]/1e9, [y y], '-', 'Color', linecol(i), 'LineWidth', 1.6);
  plot(ax, Ri.E/1e9, y, 'o', 'Color', linecol(i), 'MarkerSize', 6, ...
    'MarkerFaceColor', linecol(i));
  rows(end+1,:) = {shortname(L(i).name), y}; %#ok<AGROW>
end

% three-line mean
Ec = nan(1,nL);
for i = find(con), Ec(i) = R{CASE,i}.E; end
Ec = Ec(isfinite(Ec));
y = y + 1;
plot(ax, (mean(Ec)+[-1 1]*std(Ec))/1e9, [y y], 'k-', 'LineWidth', 2.2);
plot(ax, mean(Ec)/1e9, y, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
rows(end+1,:) = {sprintf('mean of %d', numel(Ec)), y};

% excluded lines: whatever interval they have, in grey, labelled by why
for i = 1:nL
  if con(i) || isempty(R{CASE,i}), continue; end
  y = y + 1; Ri = R{CASE,i};
  if isfinite(Ri.E_hi)
    plot(ax, [Ri.E_lo Ri.E_hi]/1e9, [y y], '--', 'Color', grey, 'LineWidth', 1.3);
    lbl = sprintf('%s (misfits)', shortname(L(i).name));
  else
    plot(ax, [Ri.E_lo/1e9 XL(2)*0.92], [y y], '--', 'Color', grey, 'LineWidth', 1.3);
    lbl = sprintf('%s (unbounded)', shortname(L(i).name));
  end
  plot(ax, Ri.E/1e9, y, 'o', 'Color', grey, 'MarkerSize', 6);
  rows(end+1,:) = {lbl, y};
end

% the same constrained lines with the strain left out, then under
% BedMachine thickness: the two rows that say what the englacial data and
% the bed product each contribute to the headline number
Ex = nan(1,nL); Eb = nan(1,nL);
for i = find(con)
  if size(R,1) >= AXCASE && ~isempty(R{AXCASE,i}), Ex(i) = R{AXCASE,i}.E; end
  if size(R,1) >= BMCASE && ~isempty(R{BMCASE,i}), Eb(i) = R{BMCASE,i}.E; end
end
Ex = Ex(isfinite(Ex)); Eb = Eb(isfinite(Eb));
if numel(Ex) >= 2
  y = y + 1;
  plot(ax, (mean(Ex)+[-1 1]*std(Ex))/1e9, [y y], '-', 'Color', grey, 'LineWidth', 1.6);
  plot(ax, mean(Ex)/1e9, y, 'o', 'Color', grey, 'MarkerSize', 6);
  rows(end+1,:) = {'a(x) only', y};
end
if numel(Eb) >= 2
  y = y + 1;
  plot(ax, (mean(Eb)+[-1 1]*std(Eb))/1e9, [y y], '-', 'Color', grey, 'LineWidth', 1.6);
  plot(ax, mean(Eb)/1e9, y, 'd', 'Color', grey, 'MarkerSize', 6, ...
    'MarkerFaceColor', grey);
  rows(end+1,:) = {'BedMachine h', y};
end

% published estimates
y = y + 1; ysep = y - 0.5;
plot(ax, [1.1 6.1], [y y], 'k-', 'LineWidth', 1.4);
plot(ax, 3.6, y, 'kd', 'MarkerSize', 7);
rows(end+1,:) = {'Elgart+ 2025', y};

y = y + 1;
plot(ax, 1.8,  y, 'kv', 'MarkerSize', 7);
plot(ax, 9.33, y, 'k^', 'MarkerSize', 7);
rows(end+1,:) = {'Sayag & Worster 2013', y};

y = y + 1;
plot(ax, [0.53 1.23], [y y], 'k-', 'LineWidth', 1.4);
plot(ax, 0.88, y, 'kd', 'MarkerSize', 7);
rows(end+1,:) = {'Vaughan 1995', y};

y = y + 1;
plot(ax, 9, y, 'kp', 'MarkerSize', 9, 'MarkerFaceColor', 'k');
rows(end+1,:) = {'laboratory ice', y};

plot(ax, XL, [ysep ysep], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);

set(ax,'XScale','log','YDir','reverse');
xlim(ax, XL); ylim(ax, [0.4 y+0.6]);
set(ax,'YTick', cell2mat(rows(:,2)), 'YTickLabel', rows(:,1));
grid(ax,'on'); box(ax,'on');
set(ax,'YGrid','off');
xlabel(ax,'E* (GPa)');
title(ax,'(c)  Effective Young''s modulus');

if ~exist(opts.out_dir,'dir'), mkdir(opts.out_dir); end
out_fn = fullfile(opts.out_dir,'EAGER_2022_elasticity_results.png');
print(hf, out_fn, '-dpng', '-r140');
fprintf('Wrote %s\n', out_fn);

end

%% ========================================================================
function s = shortname(name)
s = strrep(name, 'EAGER_2022_', '');
end

%% ========================================================================
function xs = site_xsea(L, slat, slon)
%SITE_XSEA Seaward coordinate of a site, projected onto the line's blocks.
%   Nearest-segment projection in a local metric frame, interpolating
%   x_sea along the segment - not nearest-block snapping, whose 250 m
%   quantisation would be visible at this panel's scale.
ok = isfinite(L.lat) & isfinite(L.lon) & isfinite(L.x_sea);
la = L.lat(ok); lo = L.lon(ok); xv = L.x_sea(ok);
xs = NaN;
if numel(xv) < 2, return; end
lat0 = mean(la); lon0 = mean(lo);
bx = (lo - lon0)*cosd(lat0)*111320; by = (la - lat0)*110540;
px = (slon - lon0)*cosd(lat0)*111320; py = (slat - lat0)*110540;
best = inf;
for k = 1:numel(bx)-1
  vx = bx(k+1)-bx(k); vy = by(k+1)-by(k);
  L2 = vx^2 + vy^2;
  if L2 <= 0, continue; end
  t = ((px-bx(k))*vx + (py-by(k))*vy) / L2;
  t = min(max(t,0),1);
  d = hypot(px - (bx(k)+t*vx), py - (by(k)+t*vy));
  if d < best
    best = d;
    xs = xv(k) + t*(xv(k+1)-xv(k));
  end
end
end
