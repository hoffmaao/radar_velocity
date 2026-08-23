function OUT = flexure_inversion(opts)
%FLEXURE_INVERSION Figure: inverting tidal flexure for the effective modulus.
%
%   Four panels, in the order the argument has to be made:
%
%     (a) THE OBSERVATION AND THE FIT. Surface tidal admittance a(x) per
%         line with its fit errors, the fitted beam through it, and the
%         clamp each fit chose. This is the panel that says whether an
%         elastic beam describes the data at all.
%     (b) WHERE THE ANSWER IS NOT UNIQUE. The misfit surface over E* and
%         clamp position for one line. It is a narrow diagonal valley, not
%         a bowl: a stiffer beam clamped further landward reproduces much
%         the same curve over a short window. Drawing it is the honest way
%         to present an ill-posed inversion, and it is why a grid search
%         plus a pattern search is used rather than a local optimiser
%         started somewhere hopeful.
%     (c) WHAT THE DATA DO BOUND. The misfit profiled over clamp position,
%         in units of the fitted data variance, so the one-sigma interval
%         is where each curve crosses 1. A curve that never crosses is a
%         modulus this geometry does not constrain, and says so.
%     (d) THE THICKNESS LEVER. Flexure constrains E*h^3, so the modulus
%         moves as h^-3 and nothing in the flexure data can separate the
%         two. Plotted against the published range, this is usually the
%         largest term in the error budget and it belongs on the figure
%         rather than in a footnote.
%
%   Fits come from scripts/diagnostics/elastic_modulus.m so there is one
%   implementation of the inversion and the figure cannot drift from the
%   diagnostic. Panel (d) refits over a thickness sweep on a coarser
%   modulus grid. The two files are deliberately NOT both called
%   elastic_modulus - same-named files on the MATLAB path shadow each
%   other, which this project has already lost time to once.
%
%   opts, all optional:
%     .mp_dir, .products, .block   passed through to the diagnostic
%     .out_dir     where the png goes
%     .h_sweep     thicknesses for panel (d) [m]
%     .fit         a previously returned OUT, to redraw without refitting
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "addpath('<code>/scripts/figures'); flexure_inversion"

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root);                                   % +vdef
addpath(fullfile(root,'scripts','diagnostics')); % the inversion driver

if ~isfield(opts,'out_dir') || isempty(opts.out_dir)
  opts.out_dir = '/kucresis/scratch/hoffmana_sta/vvel/figures_flexure';
end
if ~isfield(opts,'h_sweep') || isempty(opts.h_sweep)
  opts.h_sweep = 200:25:400;
end
if ~isfield(opts,'fit'), opts.fit = []; end

%% Fits
if isempty(opts.fit)
  dopts = struct();
  for f = {'mp_dir','products','block'}
    if isfield(opts, f{1}), dopts.(f{1}) = opts.(f{1}); end
  end
  D = elastic_modulus(dopts);
else
  D = opts.fit;
end
L = D.lines; R = D.fits; CASE = 1;               % primary (joint) case
nL = numel(L);
OUT = D;
avg_w = D.block * 2.5;                           % the width the fits used

%% Thickness sweep for panel (d)
% A coarser modulus grid than the diagnostic uses: this panel is about the
% SLOPE of E* against h, which is a factor of 30 across the sweep, so it
% does not need the diagnostic's resolution and would cost four minutes if
% it had it. The sweep repeats the SAME fit as the primary case - strain
% included when the line has it - so the curve passes through the
% headline point rather than through a shape-only cousin of it.
E_SWEEP = logspace(log10(0.05e9), log10(60e9), 31);
Esw = nan(numel(opts.h_sweep), nL);
for i = 1:nL
  if isempty(R{CASE,i}), continue; end
  ok = isfinite(L(i).a) & isfinite(L(i).a_std) & L(i).a_std > 0;
  for q = 1:numel(opts.h_sweep)
    qopts = struct('h', opts.h_sweep(q), 'sigma', L(i).a_std(ok), ...
                   'avg_width', avg_w, 'E_grid', E_SWEEP, ...
                   'x0_grid', R{CASE,i}.x0_grid);
    qopts = with_line_strain(qopts, L(i));
    Rq = vdef.invertElasticModulus(L(i).x_sea(ok), L(i).a(ok), qopts);
    Esw(q,i) = Rq.E;
  end
  fprintf('thickness sweep done for %s\n', shortname(L(i).name));
end
OUT.h_sweep = opts.h_sweep;
OUT.E_sweep = Esw;

%% Draw
cols = [0.11 0.42 0.69; 0.89 0.47 0.10; 0.20 0.60 0.25; 0.75 0.20 0.30];
hf = figure('Visible','off','Position',[100 100 1150 780]);
set(0,'CurrentFigure',hf);

% ---- (a) observation and fit
ax = axes('parent',hf,'Position',[0.07 0.58 0.40 0.36]);
hold(ax,'on'); hleg = []; lleg = {}; xlo = inf; xhi = -inf;
for i = 1:nL
  cc = cols(mod(i-1,size(cols,1))+1,:);
  ok = isfinite(L(i).a);
  for b = find(ok(:).')
    if ~isfinite(L(i).a_std(b)), continue; end
    plot(ax, [1 1]*L(i).x_sea(b)/1e3, ...
      L(i).a(b)+[-1 1]*L(i).a_std(b), '-', 'Color', cc, 'LineWidth', 1);
  end
  h1 = plot(ax, L(i).x_sea(ok)/1e3, L(i).a(ok), 'o', 'Color', cc, ...
    'MarkerSize', 4, 'MarkerFaceColor', cc);
  hleg(end+1) = h1; lleg{end+1} = shortname(L(i).name); %#ok<AGROW>
  Ri = R{CASE,i};
  if isempty(Ri), continue; end
  plot(ax, Ri.x_grid/1e3, Ri.w_grid, '-', 'Color', cc, 'LineWidth', 1.4);
  plot(ax, Ri.x0/1e3, 0, '^', 'Color', cc, 'MarkerSize', 7, 'MarkerFaceColor', cc);
  xlo = min([xlo; Ri.x0/1e3 - 0.4]);
  xhi = max([xhi; max(L(i).x_sea(ok))/1e3 + 0.4]);
end
if isfinite(xlo) && xhi > xlo, xlim(ax, [xlo xhi]); end
grid(ax,'on'); box(ax,'on');
xlabel(ax,'Seaward distance (km)');
ylabel(ax,'a(x) (line-mean units)');
title(ax,'(a)  Tidal admittance and fitted beam');
legend(ax, hleg, lleg, 'Location','SouthEast');
yl = get(ax,'YLim');
text(xlo+0.15, yl(1)+0.92*diff(yl), 'triangles: fitted clamp', ...
  'parent', ax, 'FontSize', 9);

% ---- (b) misfit surface for one line
ib = find(~cellfun(@isempty, R(CASE,:)), 1);
ax = axes('parent',hf,'Position',[0.57 0.58 0.36 0.36]);
if ~isempty(ib)
  % Refit this one line on a grid fine enough to RESOLVE the valley. The
  % diagnostic's 250 m clamp steps are ample for finding the minimum, which
  % it refines off-grid afterwards, but they are wider than the valley
  % itself - so drawing that surface directly gives disconnected blocks of
  % colour wherever the valley slips between grid rows. That is aliasing,
  % and on this panel it would be read as rival minima.
  R0 = R{CASE,ib};
  okb = isfinite(L(ib).a) & isfinite(L(ib).a_std) & L(ib).a_std > 0;
  bopts = struct('h', R0.h_spec, 'sigma', L(ib).a_std(okb), 'avg_width', avg_w, ...
                 'E_grid', logspace(log10(R0.E/8), log10(R0.E*8), 61), ...
                 'x0_grid', R0.x0 + (-1500:50:1500));
  bopts = with_line_strain(bopts, L(ib));
  Rb = vdef.invertElasticModulus(L(ib).x_sea(okb), L(ib).a(okb), bopts);
  % Relative to the minimum and in units of the data variance, capped so the
  % valley is legible rather than one dark pixel beside a saturated field.
  ZCAP = 30;
  Z = (Rb.J - min(Rb.J(:))) / max(Rb.s2, eps);
  imagesc(log10(Rb.E_grid/1e9), Rb.x0_grid/1e3, min(Z, ZCAP), 'parent', ax);
  set(ax,'YDir','normal'); hold(ax,'on');
  % The trade-off curve itself: the best clamp at each modulus. It runs
  % straight down the floor of the valley, which is what makes the point.
  plot(ax, log10(Rb.E_grid/1e9), Rb.x0_profile/1e3, 'w-', 'LineWidth', 1.2);
  plot(ax, log10(Rb.E/1e9), Rb.x0/1e3, 'w+', 'MarkerSize', 12, 'LineWidth', 2);
  caxis(ax,[0 ZCAP]); colorbar('peer', ax);
  % Zoom to where the misfit is still on scale, padded. The search grid is
  % deliberately far wider than the answer - it has to be, to show the
  % minimum is interior - but drawn in full a precise fit is one dark pixel
  % in a saturated field, and the valley that is the point of the panel
  % disappears. The full grid was searched either way.
  [ry, rx] = find(Z < ZCAP);
  if numel(rx) >= 2
    ex = log10(Rb.E_grid([min(rx) max(rx)])/1e9);
    xy = Rb.x0_grid([min(ry) max(ry)])/1e3;
    px = max(0.15, 0.35*diff(ex)); py = max(0.4, 0.35*diff(xy));
    xlim(ax, [ex(1)-px, ex(2)+px]);
    ylim(ax, [xy(1)-py, xy(2)+py]);
  end
  xlabel(ax,'log_{10} E* (GPa)'); ylabel(ax,'Clamp position (km)');
  title(ax, sprintf('(b)  Misfit surface, %s', shortname(L(ib).name)));
end

% ---- (c) profiled misfit
ax = axes('parent',hf,'Position',[0.07 0.09 0.40 0.36]);
hold(ax,'on'); elo = inf; ehi = -inf;
for i = 1:nL
  Ri = R{CASE,i};
  if isempty(Ri) || ~isfinite(Ri.s2) || Ri.s2 <= 0, continue; end
  cc = cols(mod(i-1,size(cols,1))+1,:);
  plot(ax, Ri.E_grid/1e9, (Ri.J_profile - min(Ri.J_profile))/Ri.s2, '-', ...
    'Color', cc, 'LineWidth', 1.4);
  plot(ax, Ri.E/1e9, 0, 'v', 'Color', cc, 'MarkerSize', 7, 'MarkerFaceColor', cc);
  lo = Ri.E_lo; if ~isfinite(lo), lo = Ri.E/3; end
  hi = Ri.E_hi; if ~isfinite(hi), hi = Ri.E*3; end
  elo = min([elo; lo/2e9]); ehi = max([ehi; hi*2/1e9]);
end
if isfinite(elo) && ehi > elo
  plot(ax, [elo ehi], [1 1], 'k--', 'LineWidth', 1);
  text(elo*1.1, 1.5, 'one sigma', 'parent', ax, 'FontSize', 9);
  xlim(ax, [elo ehi]);
end
set(ax,'XScale','log'); grid(ax,'on'); box(ax,'on'); ylim(ax,[0 9]);
xlabel(ax,'E* (GPa)'); ylabel(ax,'\Delta misfit / data variance');
title(ax,'(c)  What the data bound');

% ---- (d) the thickness lever
ax = axes('parent',hf,'Position',[0.57 0.09 0.36 0.36]);
hold(ax,'on');
% Published range for context: Elgart, Minchew and Meyer (2025) report
% 0.6-9 GPa across three Ross Ice Shelf sites, mean 3.6 +/- 2.5.
xr = [min(opts.h_sweep) max(opts.h_sweep)];
fill(ax, [xr fliplr(xr)], [0.6 0.6 9 9], [0.85 0.85 0.85], ...
  'EdgeColor','none');
for i = 1:nL
  if all(~isfinite(Esw(:,i))), continue; end
  cc = cols(mod(i-1,size(cols,1))+1,:);
  plot(ax, opts.h_sweep, Esw(:,i)/1e9, '-', 'Color', cc, 'LineWidth', 1.4);
end
set(ax,'YScale','log'); grid(ax,'on'); box(ax,'on');
xlim(ax, xr);
ylim(ax, [max(0.05, min(Esw(:))/2e9) max(60, max(Esw(:))*2/1e9)]);
yl = get(ax,'YLim');
for hh = [265 295]
  plot(ax, [hh hh], yl, 'k:', 'LineWidth', 1);
end
text(266, yl(2)/1.5, 'ApRES', 'parent', ax, 'FontSize', 9);
text(296, yl(2)/3.5, 'radar', 'parent', ax, 'FontSize', 9);
xlabel(ax,'Assumed ice thickness (m)'); ylabel(ax,'E* (GPa)');
title(ax,'(d)  E* moves as h^{-3}');
text(xr(2)-5, 0.75, 'Elgart+ 2025, three Ross sites', 'parent', ax, ...
  'FontSize', 9, 'HorizontalAlignment', 'right');

if ~exist(opts.out_dir,'dir'), mkdir(opts.out_dir); end
out_fn = fullfile(opts.out_dir,'EAGER_2022_elastic_modulus.png');
print(hf, out_fn, '-dpng', '-r140');
fprintf('Wrote %s\n', out_fn);

end

%% ========================================================================
function s = shortname(name)
s = strrep(name, 'EAGER_2022_', '');
end

%% ========================================================================
function o = with_line_strain(o, Li)
%WITH_LINE_STRAIN Attach a line's englacial strain dataset, when it has one.
if isfield(Li,'s_x') && numel(Li.s_x) >= 3
  o.strain = struct('x', Li.s_x, 'y', Li.s_y, 'sigma', Li.s_sig, ...
                    'ref_depth', 100, 'avg_width', o.avg_width);
end
end
