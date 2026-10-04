function tidal_phase_figure(opts)
%TIDAL_PHASE_FIGURE Is the strain in phase with the tide? The quadrature fit, drawn.
%
%   (a) Every near-hinge cell of the reference-5 m stack as its in-phase
%       a and quadrature b response (vdef.tidalStack, opts.dq), with the
%       jackknife sigmas. An elastic column puts every cell on b = 0; a
%       response lagging the tide by tau puts them on the line at angle
%       -2*pi*tau/T, drawn for +/-4 h.
%   (b) The lag per line (chi2 interval, too narrow: cells share passes)
%       and across the four lines (mean and standard error from their
%       scatter, the value to quote), for both references
%       (scripts/diagnostics/tidal_phase.m). Relative to 60 m the in-phase
%       response is barely above noise and the lag is not constrained.
%   (c) What the reference-5 m lag means over one diurnal cycle: the tide, the
%       strain (in phase with it if elastic) and the strain rate, a quarter
%       period ahead of the strain.
%
%   Reads tidal_stack_nozc.mat and tidal_stack_nozc_ref60.mat from
%   vdef.figureDir; writes EAGER_2022_tidal_phase.png there.

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here))); addpath(here); addpath(fullfile(fileparts(here), 'diagnostics'));
if ~isfield(opts, 'out_dir') || isempty(opts.out_dir), opts.out_dir = vdef.figureDir(); end
T = 24.84;
P = tidal_phase(struct('out_dir', opts.out_dir, 'period_h', T));
ink = [0 0 0]; soft = [0.5 0.5 0.5];
lc = [0.165 0.471 0.839; 0.85 0.45 0.10; 0.20 0.60 0.30; 0.698 0.094 0.169];   % GL1-GL4
rc = [0.55 0.55 0.55; 0 0 0];                                                     % ref 5, ref 60

[h, GRL] = grl_figure(170, 70); set(0, 'CurrentFigure', h);

% (a) b against a, reference 5 m
ax = axes('parent', h, 'Position', [0.07 0.17 0.25 0.70]); hold(ax, 'on');
lim = 8; tl = [-4 4];
for t = tl
  ang = -2*pi*t/T; plot(ax, lim*[-1 1]*cos(ang), lim*[-1 1]*sin(ang), ':', 'Color', soft, 'LineWidth', 0.6);
  e = sign(sin(ang));                                                 % label the upper end of each line
  text(ax, 0.8*lim*cos(ang)*e, 0.8*lim*sin(ang)*e, sprintf('%+d h', t), 'FontSize', 6.5, 'Color', soft, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'BackgroundColor', 'w', 'Margin', 0.5);
end
plot(ax, lim*[-1 1], [0 0], '-', 'Color', ink, 'LineWidth', 0.8);
plot(ax, [0 0], lim*[-1 1], '-', 'Color', soft, 'LineWidth', 0.5);
S = load(fullfile(opts.out_dir, 'tidal_stack_nozc.mat')); OUT = S.OUT;
hl = gobjects(1, numel(OUT));
for k = 1:numel(OUT)
  O = OUT(k);
  sel = O.zsel(:) >= 30 & (O.x_sea(:).'/1e3 >= -0.6 & O.x_sea(:).'/1e3 <= 1.5);
  a = 1e3*O.a_q(sel); b = 1e3*O.b_q(sel); sa = 1e3*O.a_q_sd(sel); sb = 1e3*O.b_q_sd(sel);
  g = isfinite(a) & isfinite(b) & sa < 3 & sb < 3;                 % the cells that constrain the angle
  for i = find(g).'
    plot(ax, a(i) + sa(i)*[-1 1], b(i)*[1 1], '-', 'Color', [lc(k,:) 0.35], 'LineWidth', 0.5);
    plot(ax, a(i)*[1 1], b(i) + sb(i)*[-1 1], '-', 'Color', [lc(k,:) 0.35], 'LineWidth', 0.5);
  end
  hl(k) = plot(ax, a(g), b(g), 'o', 'MarkerSize', 3, 'MarkerFaceColor', lc(k,:), 'MarkerEdgeColor', 'none');
end
set(ax, 'XLim', lim*[-1 1], 'YLim', lim*[-1 1], 'DataAspectRatio', [1 1 1], 'FontSize', 7.5, 'Box', 'off', ...
  'TickDir', 'out', 'Layer', 'top');
xlabel(ax, 'In-phase response (mm/m)', 'FontSize', 7.5); ylabel(ax, 'Quadrature response (mm/m)', 'FontSize', 7.5);
title(ax, '(a)  Near hinge, reference 5 m', 'FontWeight', 'normal', 'FontSize', 8);
legend(ax, hl, strrep({OUT.name}, 'EAGER_2022_', ''), 'Location', 'south', 'NumColumns', 4, 'FontSize', 6.5, 'Box', 'off');

% (b) lag per line and across the lines, both references
ax = axes('parent', h, 'Position', [0.42 0.17 0.22 0.70]); hold(ax, 'on');
yl = [-6.5 6.5];
plot(ax, [0.4 6.6], [0 0], '-', 'Color', ink, 'LineWidth', 0.8);
names = [{P(1).line.name} {'mean'}];
hr = gobjects(1, numel(P));
for f = 1:numel(P)
  E = P(f).line; dx = (f - 1.5)*0.25;
  for k = 1:numel(E)
    lo = max(E(k).tau_lo, yl(1)); hi = min(E(k).tau_hi, yl(2));
    plot(ax, k + dx*[1 1], [lo hi], '-', 'Color', rc(f,:), 'LineWidth', 1.1);
    hr(f) = plot(ax, k + dx, E(k).tau_h, 'o', 'MarkerSize', 4, 'MarkerFaceColor', rc(f,:), 'MarkerEdgeColor', 'w');
  end
  x = numel(E) + 2 + dx; c = P(f).across;                        % gap before the across-line mean
  plot(ax, [x x], c.tau_h + c.se_h*[-1 1], '-', 'Color', rc(f,:), 'LineWidth', 1.6);
  plot(ax, x, c.tau_h, 's', 'MarkerSize', 5, 'MarkerFaceColor', rc(f,:), 'MarkerEdgeColor', 'w');
end
set(ax, 'XLim', [0.4 6.6], 'YLim', yl, 'XTick', [1:numel(names)-1 numel(names)+1], 'XTickLabel', names, ...
  'FontSize', 7.5, 'Box', 'off', 'TickDir', 'out', 'Layer', 'top');
ylabel(ax, 'Response lag (h)', 'FontSize', 7.5);
title(ax, '(b)  Lag behind the tide', 'FontWeight', 'normal', 'FontSize', 8);
legend(ax, hr, {P.label}, 'Location', 'southwest', 'FontSize', 6.5, 'Box', 'off');

% (c) one diurnal cycle with the reference-5 m across-line lag
tau = P(1).across.tau_h; tlo = tau - P(1).across.se_h; thi = tau + P(1).across.se_h;
ax = axes('parent', h, 'Position', [0.72 0.17 0.26 0.70]); hold(ax, 'on');
t = linspace(0, T, 400); w = 2*pi/T;
fill(ax, [t fliplr(t)], [sin(w*(t - tlo)) fliplr(sin(w*(t - thi)))], [0.165 0.471 0.839], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
h2 = plot(ax, t, sin(w*(t - tau)), '-', 'Color', [0.165 0.471 0.839], 'LineWidth', 2.2);
h1 = plot(ax, t, sin(w*t), ':', 'Color', ink, 'LineWidth', 1.0);         % on top: the strain hides it if elastic
h3 = plot(ax, t, cos(w*(t - tau)), '--', 'Color', [0.698 0.094 0.169], 'LineWidth', 1.2);
plot(ax, [0 T], [0 0], '-', 'Color', soft, 'LineWidth', 0.5);
set(ax, 'XLim', [0 T], 'YLim', [-1.25 1.25], 'XTick', 0:6:24, 'YTick', -1:1:1, 'FontSize', 7.5, 'Box', 'off', ...
  'TickDir', 'out', 'Layer', 'top');
xlabel(ax, 'Elapsed time (h)', 'FontSize', 7.5); ylabel(ax, 'Normalized amplitude', 'FontSize', 7.5);
title(ax, sprintf('(c)  Strain rate leads the tide by %.1f h', T/4 - tau), 'FontWeight', 'normal', 'FontSize', 8);
legend(ax, [h1 h2 h3], {'tide', 'strain', 'strain rate'}, 'Location', 'southwest', 'FontSize', 6.5, 'Box', 'off');

out_fn = fullfile(opts.out_dir, 'EAGER_2022_tidal_phase.png');
print(h, out_fn, '-dpng', sprintf('-r%d', GRL.dpi)); close(h);
fprintf('Wrote %s\n', out_fn);
end
