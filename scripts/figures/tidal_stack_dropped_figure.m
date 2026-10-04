function tidal_stack_dropped_figure(opts)
%TIDAL_STACK_DROPPED_FIGURE What excluding the bad acquisitions does to the answer.
%
%   THE COMPARISON THIS DRAWS. The same coherent stack, run twice: with
%   every pass the existing gates allow, and again with the acquisitions
%   that scripts/diagnostics/pass_quality.m flagged removed entirely.
%   Dropping a pass removes pairs that ALREADY PASSED the coalignment
%   floor, which is the point - a pass that cannot be registered against
%   most of its partners has no reason to be correctly registered against
%   the rest, and a marginal alignment is applied rather than refused,
%   which is worse than a gap.
%
%   HOW TO READ IT. Per line, the response at 100 m along track, baseline
%   and dropped, each with its bootstrap band. Where the bands overlap
%   everywhere, the bad passes were not driving the answer and the
%   baseline number stands. Where they separate, the published number
%   depended on acquisitions that cannot be registered, and the dropped
%   curve is the better estimate - at the cost of the pairs it removed,
%   which the panel title reports.
%
%   A curve that moves by less than its own band is not evidence of a
%   change. The title carries the median shift in units of the combined
%   sigma so that judgement does not have to be made by eye.
%
%   Reads tidal_stack.mat (baseline) and tidal_stack_dropped.mat.

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));             % +vdef
addpath(here);                                   % grl_figure
if ~isfield(opts,'out_dir') || isempty(opts.out_dir), opts.out_dir = vdef.figureDir(); end
A = load(fullfile(opts.out_dir,'tidal_stack.mat'));          BASE = A.OUT;
B = load(fullfile(opts.out_dir,'tidal_stack_dropped.mat'));  DROP = B.OUT;

order = {'EAGER_2022_GL1','EAGER_2022_GL2','EAGER_2022_GL3','EAGER_2022_GL4'};
lab   = {'GL1','GL2','GL3','GL4'};
have  = arrayfun(@(k) any(strcmp({BASE.name}, order{k})) && any(strcmp({DROP.name}, order{k})), 1:numel(order));
order = order(have); lab = lab(have); nL = numel(order);
cbase = [0.45 0.45 0.45];                         % kept: neutral, it is the reference
cdrop = [0 0.447 0.698];                          % the new estimate carries the hue
ink = [0 0 0]; soft = [0.62 0.62 0.62];

[h, GRL] = grl_figure(170, 82); set(0,'CurrentFigure',h);
w = 0.196; x0 = 0.068; gap = 0.234;
hb = []; hd = [];
for k = 1:nL
  O = BASE(strcmp({BASE.name}, order{k}));
  P = DROP(strcmp({DROP.name}, order{k}));
  [~, i100b] = min(abs(O.zsel-100)); [~, i100d] = min(abs(P.zsel-100));
  [xb, ob] = sort(O.x_sea/1e3); [xd, od] = sort(P.x_sea/1e3);
  ab = 1e3*O.a(i100b,ob); sb = 1e3*O.a_sd(i100b,ob);
  ad = 1e3*P.a(i100d,od); sd = 1e3*P.a_sd(i100d,od);

  ax = axes('parent',h,'Position',[x0+(k-1)*gap 0.235 w 0.50]); hold(ax,'on');
  plot(ax, [min(xb) max(xb)], [0 0], '-', 'Color', soft, 'LineWidth', 0.8);
  fill(ax, [xb; flipud(xb)], [(ab-sb).'; flipud((ab+sb).')], cbase, ...
    'FaceAlpha', 0.16, 'EdgeColor','none');
  fill(ax, [xd; flipud(xd)], [(ad-sd).'; flipud((ad+sd).')], cdrop, ...
    'FaceAlpha', 0.16, 'EdgeColor','none');
  hb = plot(ax, xb, ab, '-', 'Color', cbase, 'LineWidth', 1.3);
  hd = plot(ax, xd, ad, '-', 'Color', cdrop, 'LineWidth', 1.8, ...
    'Marker','o','MarkerSize',3,'MarkerFaceColor',cdrop,'MarkerEdgeColor','w');
  xlim(ax, [min(xb)-0.25 max(xb)+0.25]); ylim(ax, [-9 6]);
  grid(ax,'on'); set(ax,'GridAlpha',0.15,'Box','off','FontSize',7.5);
  xlabel(ax,'Seaward distance (km)','FontSize',7.5);
  if k == 1, ylabel(ax,'a at 100 m (mm per m)','FontSize',7.5); end

  % median shift in units of the combined sigma, on the shared grid
  n = min(numel(xb), numel(xd));
  if numel(xb) == numel(xd) && max(abs(xb-xd)) < 1e-6
    z = abs(ab-ad) ./ max(sqrt(sb.^2+sd.^2), eps);
    zs = sprintf('shift %.2f sigma', median(z,'omitnan'));
  else
    zs = 'grids differ';
  end
  dn = O.npair - P.npair;
  title(ax, {sprintf('(%s)  %s', char('a'+k-1), lab{k}), ...
             sprintf('\\rm\\fontsize{6.5}%d of %d pairs cut', dn, O.npair), ...
             sprintf('\\rm\\fontsize{6.5}%s', zs)}, ...
    'FontWeight','normal','FontSize',8);
end
lg = legend([hb hd], {'all passes the gates allow', 'flagged acquisitions dropped'}, ...
  'Position', [0.30 0.875 0.42 0.085], 'FontSize', 7.5);
set(lg,'Box','off');
annotation(h, 'textbox', [0.068 0.0 0.90 0.09], 'String', ...
  ['Bands are bootstrap 1 sigma over pairs.  Dropping a pass removes pairs that already PASSED the ' ...
   'coalignment floor, so this tests whether a badly registered acquisition contaminated the pairs it did pass.'], ...
  'EdgeColor','none','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','bottom');
out_fn = fullfile(opts.out_dir, 'EAGER_2022_tidal_stack_dropped.png');
print(h, out_fn, '-dpng', sprintf('-r%d', GRL.dpi)); close(h);
fprintf('Wrote %s\n', out_fn);
end
