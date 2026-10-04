function tidal_stack_blocksweep_figure(opts)
%TIDAL_STACK_BLOCKSWEEP_FIGURE How much of the tidal response is the block length?
%
%   Top row, one panel per profile: the response at 100 m along the line,
%   drawn once per along-track block length. The curves are one hue,
%   light for the shortest block and dark for the longest, because block
%   length is a magnitude and not an identity. Where the curves cross zero
%   in different places, the reported sign is a property of the averaging
%   window rather than of the ice.
%
%   Bottom row: (f) the fraction of the along-track grid whose sign agrees
%   with the 500 m blocks the products use, (g) the line-median response
%   with the amplitude weighting removed for comparison, and (h) the
%   median collinearity between the measured residual misalignment and the
%   tide. Panel (h) is the mechanism: averaging over a longer block
%   removes the along-track variation in the residual, which is the only
%   thing that separates the artefact from the signal, so the artefact
%   becomes less separable as the blocks grow.
%
%   Reads scripts/diagnostics/tidal_stack_blocksweep.m's saved output.

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));             % +vdef
addpath(here);                                   % grl_figure
if ~isfield(opts,'out_dir') || isempty(opts.out_dir), opts.out_dir = vdef.figureDir(); end
if ~isfield(opts,'fit') || isempty(opts.fit)
  S = load(fullfile(opts.out_dir, 'tidal_stack_blocksweep.mat')); OUT = S.OUT;
else, OUT = opts.fit; end

% Four lines in the section panels. EAGER_2022 is a second build of GL1,
% not a fifth profile (vdef.surveyLines), so it is drawn only in the
% summary panels and labelled as the build it is - its disagreement with
% GL1 is a measurement of the method's own error, and putting it in the
% top row would present one line twice.
order = {'EAGER_2022_GL4','EAGER_2022_GL3','EAGER_2022_GL2','EAGER_2022_GL1','EAGER_2022'};
lab   = {'GL4','GL3','GL2','GL1','GL1 build 2'};
nsect = 4;                       % how many of those get a section panel
legcol = [0.835 0.369 0; 0 0.620 0.451; 0.902 0.624 0; 0 0.447 0.698; 0 0.447 0.698];
ink = [0 0 0]; soft = [0.45 0.45 0.45]; red = [0.69 0 0.13];
have = arrayfun(@(k) any(strcmp({OUT.name}, order{k})), 1:numel(order));
order = order(have); lab = lab(have); legcol = legcol(have,:);
nL = numel(order);
nS = min(nsect, nL);
blocks = OUT(1).blocks; nk = numel(blocks);
% sequential ramp, light to dark, for a magnitude
ramp = interp1([0 1], [0.78 0.87 0.94; 0.03 0.19 0.38], linspace(0,1,nk));

[h, GRL] = grl_figure(170, 128); set(0,'CurrentFigure',h);
w = 0.200; x0 = 0.075; gap = 0.234;
for k = 1:nS
  O = OUT(strcmp({OUT.name}, order{k}));
  ax = axes('parent',h,'Position',[x0+(k-1)*gap 0.600 w 0.285]); hold(ax,'on');
  plot(ax, [O.xgrid(1) O.xgrid(end)], [0 0], '-', 'Color', soft, 'LineWidth', 0.8);
  for q = 1:nk
    plot(ax, O.xgrid, O.a_grid(:,q), '-', 'Color', ramp(q,:), 'LineWidth', 1.3);
  end
  xlim(ax, [O.xgrid(1)-0.1 O.xgrid(end)+0.1]); ylim(ax, [-8 8]);
  grid(ax,'on'); set(ax,'GridAlpha',0.15,'Box','off');
  title(ax, sprintf('(%s)  %s', char('a'+k-1), lab{k}), 'FontWeight','normal','FontSize',8);
  xlabel(ax, 'Seaward (km)', 'FontSize', 7.5);
  set(ax, 'FontSize', 7.5);
  if k == 1, ylabel(ax, 'a at 100 m (mm per m)'); else, set(ax,'YTickLabel',[]); end
end
% colour key for the ramp, clear of every panel
axk = axes('parent',h,'Position',[x0 0.955 0.30 0.014]); hold(axk,'on');
for q = 1:nk
  fill(axk, [q-1 q q q-1], [0 0 1 1], ramp(q,:), 'EdgeColor','none');
end
xlim(axk,[0 nk]); ylim(axk,[0 1]);
set(axk,'YTick',[],'XTick',(1:nk)-0.5, 'XTickLabel', ...
  arrayfun(@(b) sprintf('%d', round(b*2.5)), blocks, 'uni', 0), ...
  'FontSize', 6.5, 'TickLength',[0 0], 'Box','off', 'XTickLabelRotation', 0, ...
  'XAxisLocation','top');
text(nk*1.03, 0.5, 'block length (m)', 'parent', axk, 'FontSize', 7, ...
  'VerticalAlignment','middle', 'HorizontalAlignment','left');

pos2 = {[0.068 0.150 0.235 0.295], [0.398 0.150 0.235 0.295], [0.728 0.150 0.235 0.295]};
ttl = {'(f)  Sign agreement, 500 m ref', '(g)  Line-median response', ...
       '(h)  Artefact collinearity'};
for pnl = 1:3
  ax = axes('parent',h,'Position',pos2{pnl}); hold(ax,'on');
  for k = 1:nL
    O = OUT(strcmp({OUT.name}, order{k}));
    switch pnl
      case 1, yv = O.sign_agree;
      case 2, yv = O.a_med;
      case 3, yv = O.rcol_med;
    end
    sty = '-o'; lw = 1.5;
    if k > nS, sty = ':o'; lw = 1.2; end          % the duplicate build
    plot(ax, blocks*2.5, yv, sty, 'Color', legcol(k,:), 'LineWidth', lw, ...
      'MarkerSize', 4, 'MarkerFaceColor', legcol(k,:), 'MarkerEdgeColor','w');
    if pnl == 2
      plot(ax, blocks*2.5, O.au_med, '--', 'Color', legcol(k,:), 'LineWidth', 1.0);
    end
  end
  set(ax,'XScale','log'); xlim(ax, [blocks(1)*2.5*0.8 blocks(end)*2.5*1.25]);
  set(ax,'XTick', [63 250 1000], 'XTickLabel', {'63','250','1000'}, ...
    'XTickLabelRotation', 0, 'FontSize', 7.5);
  grid(ax,'on'); set(ax,'GridAlpha',0.15,'Box','off');
  xlabel(ax, 'Block length (m)', 'FontSize', 7.5);
  title(ax, ttl{pnl}, 'FontWeight','normal','FontSize',8);
  switch pnl
    case 1
      ylim(ax,[0 1.05]); ylabel(ax,'fraction of grid');
      plot(ax, get(ax,'XLim'), [1 1], ':', 'Color', ink, 'LineWidth', 0.8);
    case 2
      ylabel(ax,'median a (mm per m)');
      plot(ax, get(ax,'XLim'), [0 0], '-', 'Color', soft, 'LineWidth', 0.8);
    case 3
      ylim(ax,[0 1.05]); ylabel(ax,'median |r(d, tide)|');
      plot(ax, get(ax,'XLim'), [0.9 0.9], ':', 'Color', red, 'LineWidth', 1);
      text(blocks(1)*2.5*0.9, 0.80, 'not separable above', 'parent', ax, 'FontSize', 6.5, 'Color', red);
  end
  if pnl == 1
    for k = 1:nL
      text(0.03 + 0.175*(k-1), 0.09, lab{k}, 'parent', ax, 'Units','normalized', ...
        'Color', legcol(k,:), 'FontSize', 6.5, 'FontWeight','bold');
    end
  end
end
annotation(h, 'textbox', [0.068 0.0 0.90 0.042], 'String', ...
  ['Solid: the estimator, an amplitude-weighted complex block mean.  Dashed in (g): the same with unit ' ...
   'phasors.  Dotted in (f)-(h): the second build of GL1, a control on the method, not a fifth line.'], ...
   'EdgeColor','none', 'FontSize', 7, ...
  'HorizontalAlignment','left', 'VerticalAlignment','bottom');
out_fn = fullfile(opts.out_dir, 'EAGER_2022_tidal_stack_blocksweep.png');
print(h, out_fn, '-dpng', sprintf('-r%d', GRL.dpi)); close(h);
fprintf('Wrote %s\n', out_fn);
end
