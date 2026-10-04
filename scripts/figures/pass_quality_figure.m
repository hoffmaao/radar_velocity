function pass_quality_figure(opts)
%PASS_QUALITY_FIGURE Which passes are bad, and how much of the network they cost.
%
%   Top row, one panel per line: the pass-by-pass coalignment correlation,
%   the quantity the pair gate acts on. A bad PAIR is one dark cell; a bad
%   PASS is a dark row crossed with a dark column, and the difference is
%   the whole point - it says whether to fix a pair or drop an acquisition.
%   Cells below the 0.85 floor are outlined, so the gate's actual cut is
%   visible rather than inferred.
%
%   Bottom row: each pass's median correlation against all its partners,
%   with the floor drawn. Bars are ordered as the passes are, so the
%   x axis doubles as the acquisition sequence. A flagged pass is drawn in
%   the warning colour and named.
%
%   Reads scripts/diagnostics/pass_quality.m's saved output.

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));             % +vdef
addpath(here);                                   % grl_figure
if ~isfield(opts,'out_dir') || isempty(opts.out_dir), opts.out_dir = vdef.figureDir(); end
if ~isfield(opts,'fit') || isempty(opts.fit)
  S = load(fullfile(opts.out_dir, 'pass_quality.mat')); OUT = S.OUT;
else, OUT = opts.fit; end
if ~isfield(opts,'floor') || isempty(opts.floor), opts.floor = 0.85; end

order = {'EAGER_2022_GL1','EAGER_2022_GL2','EAGER_2022_GL3','EAGER_2022_GL4'};
lab   = {'GL1','GL2','GL3','GL4'};
have  = arrayfun(@(k) any(strcmp({OUT.name}, order{k})), 1:numel(order));
order = order(have); lab = lab(have); nL = numel(order);
ink = [0 0 0]; soft = [0.45 0.45 0.45]; warn = [0.835 0.369 0]; ok = [0 0.447 0.698];
% sequential ramp: magnitude, one hue, light to dark
ramp = interp1([0 1], [0.96 0.97 0.99; 0.03 0.19 0.38], linspace(0,1,256));

[h, GRL] = grl_figure(170, 118); set(0,'CurrentFigure',h);
w = 0.196; x0 = 0.062; gap = 0.235;
for k = 1:nL
  O = OUT(strcmp({OUT.name}, order{k}));
  Np = numel(O.q_med); Q = O.Q;
  ax = axes('parent',h,'Position',[x0+(k-1)*gap 0.585 w 0.30]); hold(ax,'on');
  imagesc(ax, 1:Np, 1:Np, Q, 'AlphaData', isfinite(Q));
  colormap(ax, ramp); caxis(ax, [0.5 1]);
  set(ax,'YDir','reverse'); axis(ax,'square');
  % outline the cells the gate actually cuts
  [bi, bj] = find(isfinite(Q) & Q < opts.floor);
  for q = 1:numel(bi)
    plot(ax, bj(q)+[-.5 .5 .5 -.5 -.5], bi(q)+[-.5 -.5 .5 .5 -.5], '-', ...
      'Color', warn, 'LineWidth', 0.6);
  end
  xlim(ax,[0.5 Np+0.5]); ylim(ax,[0.5 Np+0.5]);
  set(ax,'XTick',[1 5:5:Np],'YTick',[1 5:5:Np],'FontSize',7,'Box','off','TickLength',[0 0]);
  title(ax, sprintf('(%s)  %s', char('a'+k-1), lab{k}), 'FontWeight','normal','FontSize',8);
  xlabel(ax,'pass','FontSize',7.5);
  if k == 1, ylabel(ax,'pass','FontSize',7.5); end

  ax2 = axes('parent',h,'Position',[x0+(k-1)*gap 0.185 w 0.235]); hold(ax2,'on');
  for p = 1:Np
    c = ok; if O.flag(p), c = warn; end
    bar(ax2, p, O.q_med(p), 0.75, 'FaceColor', c, 'EdgeColor','none');
  end
  plot(ax2, [0.3 Np+0.7], [1 1]*opts.floor, '--', 'Color', ink, 'LineWidth', 0.9);
  % name the flagged acquisitions above the panel rather than on the bars,
  % where rotated date strings collide with everything
  fl = find(O.flag);
  if isempty(fl)
    ttl2 = {'none flagged'};  tc = soft;
  else
    % drop the year: the panels are 33 mm wide and three full date strings
    % will not fit on one line without running into the neighbour
    nmz = cellfun(@(c) c(5:end), O.seg(fl), 'uni', 0);
    if numel(nmz) > 2
      ttl2 = {strjoin(nmz(1:2), ', '), strjoin(nmz(3:end), ', ')};
    else
      ttl2 = {strjoin(nmz, ', ')};
    end
    tc = warn;
  end
  title(ax2, ttl2, 'FontWeight','normal','FontSize',6.5,'Color',tc,'Interpreter','none');
  xlim(ax2,[0.3 Np+0.7]); ylim(ax2,[0 1.05]);
  set(ax2,'XTick',[1 5:5:Np],'FontSize',7.5,'Box','off');
  grid(ax2,'on'); set(ax2,'GridAlpha',0.15);
  xlabel(ax2,'pass','FontSize',7.5);
  if k == 1, ylabel(ax2,'median correlation','FontSize',7.5); end
end
cb = colorbar(ax, 'Position', [x0 0.945 0.20 0.014], 'Orientation','horizontal');
set(cb, 'FontSize', 6.5);
annotation(h, 'textbox', [0.272 0.938 0.30 0.028], 'String', 'coalignment correlation', ...
  'EdgeColor','none', 'FontSize', 7, 'VerticalAlignment','middle');
annotation(h, 'textbox', [0.062 0.0 0.90 0.055], 'String', ...
  ['Orange outline: pair below the 0.85 floor (dashed) that the gate cuts.  Orange bar: flagged pass, named above.  ' ...
   'A pale row crossed with a pale column is a bad acquisition; a scattered pale cell is a bad pair.'], ...
  'EdgeColor','none', 'FontSize', 7, 'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
  'Interpreter','none');
out_fn = fullfile(opts.out_dir, 'EAGER_2022_pass_quality.png');
print(h, out_fn, '-dpng', sprintf('-r%d', GRL.dpi)); close(h);
fprintf('Wrote %s\n', out_fn);
end
