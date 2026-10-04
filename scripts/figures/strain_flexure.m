function strain_flexure(opts)
%STRAIN_FLEXURE Figure: the englacial strain admittance against the fitted beam.
%
%   The companion of flexure_inversion's panel (a) for the OTHER observable
%   of the beam. Per line: the radar's own column change at REF_DEPTH per
%   metre of tide, block by block with its 1-sigma, and what the joint fit's
%   beam predicts for it - dh = amp * K2(x) * w''(x), with the SAME
%   amplitude the surface admittance fixed. The surface panel says whether
%   an elastic beam describes the deflection; this one says whether the
%   same beam describes the bending strain inside the ice, and by how much
%   it does not.
%
%   READ THE CHI-SQUARED, NOT THE EYE. The strain rows misfit the beam at
%   reduced chi-squared 2-4 on three of the four lines (against 0.1-1.6
%   for the shape), so the joint fit scaled their sigmas up by 1.6-2.0
%   before weighting them (opts.rescale in vdef.invertElasticModulus).
%   The bars drawn here are the INPUT sigmas; the annotation carries the
%   raw chi-squared and the scale, so the excess over thin-plate bending
%   is visible rather than absorbed. GL2's strain fits at ~1 - it is the
%   line whose surface data needed the pass gate, not its radar.
%
%   Reads the driver's saved output rather than refitting, so the curves
%   are exactly the fits the headline numbers came from.
%
%   opts, all optional:
%     .fit       an OUT struct from scripts/diagnostics/elastic_modulus.m
%     .fit_file  a .mat holding one (default figs/flexure_fit_cats.mat)
%     .out_dir   where the png goes
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "addpath('<code>/scripts/figures'); strain_flexure"

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));             % +vdef
addpath(here);                                   % grl_figure

if ~isfield(opts,'out_dir') || isempty(opts.out_dir)
  opts.out_dir = vdef.figureDir();
end
if ~isfield(opts,'fit_file') || isempty(opts.fit_file)
  opts.fit_file = fullfile(opts.out_dir, 'flexure_fit_cats.mat');
end
if isfield(opts,'fit') && ~isempty(opts.fit)
  D = opts.fit;
else
  S = load(opts.fit_file);
  fn = fieldnames(S);
  D = S.(fn{1});
end
L = D.lines; R = D.fits; CASE = 1;               % primary (joint) case
nL = numel(L);
tide_kind = 'line-mean height';
if isfield(D,'tide') && strcmp(D.tide,'cats'), tide_kind = 'CATS2008 tide'; end

cols = [0.11 0.42 0.69; 0.89 0.47 0.10; 0.20 0.60 0.25; 0.75 0.20 0.30];
[hf, GRL] = grl_figure(170, 115.3);
set(0,'CurrentFigure',hf);
pos = {[0.08 0.60 0.40 0.31], [0.57 0.60 0.40 0.31], ...
       [0.08 0.14 0.40 0.31], [0.57 0.14 0.40 0.31]};
ylim_all = [0 0];
axs = [];
for i = 1:min(nL, 4)
  Ri = R{CASE,i};
  ax = axes('parent', hf, 'Position', pos{i}); axs(end+1) = ax; %#ok<AGROW>
  hold(ax,'on');
  cc = cols(mod(i-1,size(cols,1))+1,:);
  if isempty(Ri) || ~Ri.has_strain
    text(0.5, 0.5, 'no strain dataset', 'parent', ax, 'Units', 'normalized', ...
      'HorizontalAlignment', 'center');
    title(ax, sprintf('(%s)  %s', char('a'+i-1), shortname(L(i).name)));
    continue;
  end
  x  = Ri.strain_x/1e3;
  y  = 1e3*Ri.strain_obs;
  sg = 1e3*Ri.strain_sigma;
  plot(ax, [x(1)-1 x(end)+1], [0 0], '-', 'Color', [0.7 0.7 0.7]);
  for b = 1:numel(x)
    plot(ax, [1 1]*x(b), y(b)+[-1 1]*sg(b), '-', 'Color', cc, 'LineWidth', 1);
  end
  plot(ax, x, y, 'o', 'Color', cc, 'MarkerSize', 4, 'MarkerFaceColor', cc);
  % The beam's prediction: the continuous curve on the beam section of
  % the plotting grid, and its block averages at the observation points -
  % the block average is what the fit compared against, and on a 500 m
  % block over a ~1.2 km flexural length the two differ visibly.
  ng = numel(Ri.strain_grid);
  xg = Ri.x_grid(end-ng+1:end)/1e3;
  plot(ax, xg, 1e3*Ri.strain_grid, '-', 'Color', cc, 'LineWidth', 1.4);
  plot(ax, x, 1e3*Ri.strain_model, 's', 'Color', cc, 'MarkerSize', 5);
  plot(ax, Ri.x0/1e3, 0, '^', 'Color', cc, 'MarkerSize', 7, 'MarkerFaceColor', cc);
  xlim(ax, [min(Ri.x0/1e3, x(1)) - 0.3, x(end) + 0.3]);
  yl = [min([y - sg; 1e3*Ri.strain_grid(:)]), max([y + sg; 1e3*Ri.strain_grid(:)])];
  ylim_all = [min(ylim_all(1), yl(1)), max(ylim_all(2), yl(2))];
  grid(ax,'on'); box(ax,'on');
  sc = 1; if isfield(Ri,'sigma_scale_strain'), sc = Ri.sigma_scale_strain; end
  % Two short lines rather than one long one: at GRL width the panels
  % are 65 mm across and a single-line title runs into its neighbour.
  title(ax, {sprintf('(%s)  %s,  E* = %.2f GPa', char('a'+i-1), shortname(L(i).name), Ri.E/1e9), ...
             sprintf('\\chi^2_{strain} = %.2f (input \\sigma), scale %.2f', Ri.chi2_strain, sc)}, ...
    'FontSize', 8, 'FontWeight', 'normal');
  if i >= 3, xlabel(ax, 'Seaward distance (km)'); end
  if mod(i,2) == 1
    ylabel(ax, sprintf('dh(%.0f m) (mm per m of tide)', Ri.strain_spec.ref_depth));
  end
end
pad = 0.08*diff(ylim_all);
for ax = axs
  if diff(ylim_all) > 0, ylim(ax, ylim_all + [-pad pad]); end
end
if ~isempty(axs)
  % Figure-level key along the bottom margin, clear of every panel.
  annotation(hf, 'textbox', [0.08 0.0 0.9 0.035], 'String', ...
    ['circles: radar, bars input 1\sigma;  line: fitted beam, squares: its block means;  ' ...
     'triangle: clamp.  Tide: ' tide_kind '.'], ...
    'EdgeColor', 'none', 'FontSize', 7, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');
end

if ~exist(opts.out_dir,'dir'), mkdir(opts.out_dir); end
out_fn = fullfile(opts.out_dir, 'EAGER_2022_strain_flexure.png');
print(hf, out_fn, '-dpng', sprintf('-r%d', GRL.dpi));
fprintf('Wrote %s\n', out_fn);
end

%% ========================================================================
function s = shortname(name)
s = strrep(name, 'EAGER_2022_', '');
end
