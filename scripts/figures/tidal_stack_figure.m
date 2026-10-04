function tidal_stack_figure(opts)
%TIDAL_STACK_FIGURE The coherent all-pairs tidal response, every profile.
%   Top row: a(z, x) from the coherent stack per leg, the full depth
%   section, on the project's diverging scale. Bottom row: at 100 m, the
%   coherent estimate with its bootstrap band, the unwrapped chain on the
%   same pairs, and the artefact-controlled estimate drawn ONLY where the
%   measured residual is separable from the tide (|r| < 0.5); the panel
%   title carries the collinearity so a ridge is never mistaken for a
%   result. Reads scripts/diagnostics/tidal_stack.m's saved output.
if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));             % +vdef
addpath(here);                                   % grl_figure
if ~isfield(opts,'out_dir') || isempty(opts.out_dir)
  opts.out_dir = vdef.figureDir();
end
if ~isfield(opts,'fit') || isempty(opts.fit)
  S = load(fullfile(opts.out_dir, 'tidal_stack.mat')); OUT = S.OUT;
else, OUT = opts.fit; end
order = {'EAGER_2022_GL4','EAGER_2022_GL3','EAGER_2022_GL2','EAGER_2022_GL1'};
col = containers.Map(order, {[0.835 0.369 0],[0 0.620 0.451],[0.902 0.624 0],[0 0.447 0.698]});
hkey = []; ink = [0 0 0]; soft = [0.45 0.45 0.45]; red = [0.69 0 0.13];
dneg = [0.698 0.094 0.169]; dmid = [0.941 0.937 0.925]; dpos = [0.165 0.471 0.839];
ndiv = 256; half = ndiv/2;
dmap = [interp1([0 1],[dneg; dmid],linspace(0,1,half)); interp1([0 1],[dmid; dpos],linspace(0,1,ndiv-half))];
CL = 6;
[h, GRL] = grl_figure(170, 128); set(0,'CurrentFigure',h);
nL = numel(order); w = 0.205; x0 = 0.075; gap = 0.235;
for k = 1:nL
  ii = find(strcmp({OUT.name}, order{k}), 1);
  if isempty(ii), continue; end
  O = OUT(ii); [xs, o] = sort(O.x_sea/1e3); A = 1e3*O.a(:, o);
  ax = axes('parent',h,'Position',[x0+(k-1)*gap 0.56 w 0.29]); hold(ax,'on');
  imagesc(ax, xs, O.zsel, A); set(ax,'YDir','reverse'); colormap(ax, dmap); caxis(ax,[-CL CL]);
  plot(ax, xs([1 end]), [100 100], '--', 'Color', ink, 'LineWidth', 0.8);
  xlim(ax, [xs(1)-0.25 xs(end)+0.25]); ylim(ax, [O.zsel(1)-5 O.zsel(end)+5]);
  title(ax, sprintf('(%s)  %s: a(z, x)', char('a'+k-1), strrep(O.name,'EAGER_2022_','')), 'FontWeight','normal','FontSize',8);
  if k == 1, ylabel(ax,'Depth (m)'); else, set(ax,'YTickLabel',[]); end
  set(ax,'Box','off');
  ax2 = axes('parent',h,'Position',[x0+(k-1)*gap 0.08 w 0.32]); hold(ax2,'on');
  [~, i100] = min(abs(O.zsel-100));
  a = 1e3*O.a(i100,o); sd = 1e3*O.a_sd(i100,o); a2 = 1e3*O.a2(i100,o); rc = O.rcol(o).'; cu = 1e3*O.c_unw(o).';
  plot(ax2, xs([1 end]), [0 0], '-', 'Color', soft, 'LineWidth', 0.8);
  fill(ax2, [xs; flipud(xs)], [(a-sd).'; flipud((a+sd).')], col(O.name), 'FaceAlpha', 0.18, 'EdgeColor','none');
  hs = plot(ax2, xs, a, '-', 'Color', col(O.name), 'LineWidth', 1.8, 'Marker','o','MarkerSize',3,'MarkerFaceColor',col(O.name),'MarkerEdgeColor','w');
  hu = plot(ax2, xs, cu, '-', 'Color', ink, 'LineWidth', 1.0, 'Marker','s','MarkerSize',3,'MarkerFaceColor','w');
  sep = abs(rc) < 0.5;
  hc = plot(ax2, xs(sep), a2(sep), 'd', 'Color', col(O.name), 'MarkerSize', 5, 'MarkerFaceColor','w','LineWidth',1.2);
  if ~isempty(hc) && any(sep), hkey = [hs hu hc]; end
  xlim(ax2, [xs(1)-0.25 xs(end)+0.25]); ylim(ax2, [-9 6]); grid(ax2,'on'); set(ax2,'GridAlpha',0.15,'Box','off');
  xlabel(ax2, 'Seaward distance (km)');
  if k == 1, ylabel(ax2, 'Response at 100 m (mm per m)'); else, set(ax2,'YTickLabel',[]); end
  rmed = median(abs(rc), 'omitnan');
  if any(sep), tc = ink; ts = 'artefact separable';
  else, tc = red; ts = 'artefact NOT separable'; end
  title(ax2, {sprintf('|r(d,T)| median %.2f', rmed), ts}, 'FontWeight','normal','FontSize',7.5,'Color',tc);
end
cb = colorbar(ax, 'Position', [x0 0.95 0.25 0.014], 'Orientation','horizontal');
set(cb.Label, 'String', 'a (mm per m of tide)');
if isempty(hkey), hkey = [hs hu]; end
keys = {'coherent stack \pm bootstrap sd', 'unwrapped chain, same pairs', ...
  'artefact-controlled stack, where |r(d,T)| < 0.5'};
lg = legend(hkey, keys(1:numel(hkey)), 'Position', [0.42 0.905 0.55 0.085], 'FontSize', 7);
set(lg, 'Box', 'off');
out_fn = fullfile(opts.out_dir, 'EAGER_2022_tidal_stack.png');
print(h, out_fn, '-dpng', sprintf('-r%d', GRL.dpi)); close(h);
fprintf('Wrote %s\n', out_fn);
end
