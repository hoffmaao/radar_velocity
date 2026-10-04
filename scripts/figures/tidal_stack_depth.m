function tidal_stack_depth(opts)
%TIDAL_STACK_DEPTH The tidal response against depth, and what the reference does to it.
%
%   WHAT THIS SHOWS. Per line, the line-mean coherent-stack response a(z)
%   (mm of column change per metre of tide, trend co-estimated, weighted by
%   the per-pass jackknife sigma) for the far field (1.5 km or more seaward
%   of the hinge, where bending is a few percent of its peak) and the near
%   field (the rest), twice: with the chain's default phase reference 50 ns
%   below the surface pick (~5 m of firn) and with the reference moved to
%   60 m.
%
%   WHY IT MATTERS. With the 5 m reference every line shows a response that
%   is positive at 20 m, zero near 30 m and about -2 mm/m from 40 m to the
%   bottom of the column, in the far field as well as the near field. A
%   displacement that is the same at 40 m and at 250 m is not deformation
%   of the column: it is a change between the reference and everything
%   below it, and it accumulates between 5 and 60 m. Referencing at 60 m
%   removes it: the far field is then zero at every depth on GL1, GL3 and
%   GL4, and what remains near the hinge is the small, depth-shaped
%   relative motion of the ice column that bending predicts.
%
%   WHAT THE SHALLOW TERM IS NOT. Not a trend leak (it survives the trend
%   term), not a brine layer (no reflector in 10-60 m), not pass direction
%   (every pass of a line runs the same way), and not noise (the same
%   numbers on three independent lines). What it is cannot be settled from
%   these passes: each line's passes span under three days and the McMurdo
%   tide is diurnal, so the CATS tide at the pass times is a function of
%   the hour of day (R^2 = 0.99), and a diurnal process in the surface snow,
%   the sled, or the instrument is indistinguishable from a tidal one.
%
%   Reads tidal_stack_nozc.mat and tidal_stack_nozc_ref60.mat from
%   vdef.figureDir; writes EAGER_2022_tidal_stack_depth.png there.

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));             % +vdef
addpath(here);                                   % grl_figure
if ~isfield(opts,'out_dir') || isempty(opts.out_dir), opts.out_dir = vdef.figureDir(); end
if ~isfield(opts,'files') || isempty(opts.files), opts.files = {'tidal_stack_nozc.mat', 'tidal_stack_nozc_ref60.mat'}; end
if ~isfield(opts,'split_km') || isempty(opts.split_km), opts.split_km = 1.5; end
A = load(fullfile(opts.out_dir, opts.files{1})); S5 = A.OUT;
B = load(fullfile(opts.out_dir, opts.files{2})); S60 = B.OUT;
if ~isfield(opts,'ref_m') || isempty(opts.ref_m), opts.ref_m = [recorded_ref(S5, 5) recorded_ref(S60, 60)]; end
order = {'EAGER_2022_GL1','EAGER_2022_GL2','EAGER_2022_GL3','EAGER_2022_GL4'};
lab = {'GL1','GL2','GL3','GL4'};

col_far = [0.165 0.471 0.839]; col_near = [0.698 0.094 0.169];
ink = [0 0 0]; soft = [0.5 0.5 0.5]; firn = [0.93 0.93 0.90];
XL = [-4.5 4.5]; ZL = [0 260];

[h, GRL] = grl_figure(170, 95); set(0, 'CurrentFigure', h);
w = 0.195; x0 = 0.075; gap = 0.235; y0 = 0.19; hh = 0.66;
hl = gobjects(1, 4);
for k = 1:4
  O5 = S5(strcmp({S5.name}, order{k})); O60 = S60(strcmp({S60.name}, order{k}));
  ax = axes('parent', h, 'Position', [x0 + (k-1)*gap, y0, w, hh]); hold(ax, 'on');
  % the shallow band the 5 m reference cannot see past, and the two reference depths
  fill(ax, [XL(1) XL(2) XL(2) XL(1)], [opts.ref_m(1) opts.ref_m(1) opts.ref_m(2) opts.ref_m(2)], firn, 'EdgeColor', 'none');
  plot(ax, XL, [0 0], '-', 'Color', soft, 'LineWidth', 0.5);
  plot(ax, [0 0], ZL, '-', 'Color', soft, 'LineWidth', 0.5);
  for r = opts.ref_m
    plot(ax, XL, [r r], ':', 'Color', ink, 'LineWidth', 0.6);
  end
  for pass_ = 1:2
    if pass_ == 1, O = O5; ls_ = '-'; lw = 1.4; alpha = 0.18; else, O = O60; ls_ = '--'; lw = 1.4; alpha = 0.10; end
    x = O.x_sea(:)/1e3; z = O.zsel(:);
    for sel = 1:2
      if sel == 1, s = x >= opts.split_km; c = col_far; else, s = x < opts.split_km; c = col_near; end
      a = 1e3*O.a_t(:, s); sd = 1e3*O.a_t_sd(:, s); wgt = 1./sd.^2; wgt(~isfinite(wgt)) = 0;
      m = sum(wgt.*a, 2, 'omitnan') ./ max(sum(wgt, 2), eps); e = 1./sqrt(max(sum(wgt, 2), eps));
      ok = isfinite(m) & sum(wgt, 2) > 0;
      if pass_ == 2, ok = ok & z > opts.ref_m(2); end          % relative to 60 m: show the column below it
      fill(ax, [m(ok) - e(ok); flipud(m(ok) + e(ok))], [z(ok); flipud(z(ok))], c, 'FaceAlpha', alpha, 'EdgeColor', 'none');
      hl(sel + 2*(pass_-1)) = plot(ax, m(ok), z(ok), ls_, 'Color', c, 'LineWidth', lw);
    end
  end
  set(ax, 'YDir', 'reverse', 'XLim', XL, 'YLim', ZL, 'FontSize', 7.5, 'Box', 'off', 'TickDir', 'out', ...
    'XColor', ink, 'YColor', ink, 'Layer', 'top', 'XTick', -4:2:4);
  grid(ax, 'on'); set(ax, 'GridAlpha', 0.10);
  title(ax, sprintf('(%s)  %s,  %d pairs', char('a'+k-1), lab{k}, O5.npair), 'FontWeight', 'normal', 'FontSize', 8);
  xlabel(ax, 'Response (mm per m of tide)', 'FontSize', 7.5);
  if k == 1, ylabel(ax, 'Depth below the surface (m)', 'FontSize', 7.5); end
  text(ax, XL(1) + 0.15, mean(opts.ref_m), 'shallow firn', 'FontSize', 6.5, 'Color', soft, 'VerticalAlignment', 'middle');
  text(ax, XL(2) - 0.15, opts.ref_m(1), sprintf('ref %d m', opts.ref_m(1)), 'FontSize', 6.5, 'Color', ink, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
  text(ax, XL(2) - 0.15, opts.ref_m(2), sprintf('ref %d m', opts.ref_m(2)), 'FontSize', 6.5, 'Color', ink, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
end
lg = legend(hl, {sprintf('far field (\\geq %.1f km), reference %d m', opts.split_km, opts.ref_m(1)), ...
                 sprintf('near field (< %.1f km), reference %d m', opts.split_km, opts.ref_m(1)), ...
                 sprintf('far field, reference %d m', opts.ref_m(2)), ...
                 sprintf('near field, reference %d m', opts.ref_m(2))}, ...
  'Position', [0.30 0.905 0.42 0.085], 'FontSize', 7, 'NumColumns', 2);
set(lg, 'Box', 'off');
annotation(h, 'textbox', [x0 0.0 0.90 0.085], 'String', ...
  ['Line means of the trend-controlled coherent stack, weighted by the per-pass jackknife sigma; bands are 1 sigma. ' ...
   'With the 5 m reference every line carries a response that is the same from 40 m to the bottom of the column, ' ...
   'in the far field as in the near field: a change between the reference and the column, confined to the shallow firn. ' ...
   'Referenced at 60 m, the far field is zero and only the near-hinge column moves.'], ...
  'EdgeColor', 'none', 'FontSize', 6.8, 'Color', soft, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
out_fn = fullfile(opts.out_dir, 'EAGER_2022_tidal_stack_depth.png');
print(h, out_fn, '-dpng', sprintf('-r%d', GRL.dpi)); close(h);
fprintf('Wrote %s\n', out_fn);
end

function r = recorded_ref(OUT, legacy)
% the phase reference tidal_stack saved with the file [m]; NaN is its 50 ns
% default (~5 m), and files saved before it was recorded take the legacy value
if ~isfield(OUT, 'ref_depth'), r = legacy; return; end
r = OUT(1).ref_depth;
if isnan(r), r = 5; end
end
