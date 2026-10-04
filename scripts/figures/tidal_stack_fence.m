function tidal_stack_fence(opts)
%TIDAL_STACK_FENCE The tidal response a(x, y, z) as a 3-D fence diagram.
%
%   Each survey line is a vertical curtain standing on its own track:
%   along-survey distance and across-survey offset are the horizontal axes
%   (EPSG:3031 rotated so the survey runs along x, origin at the mean
%   grounding-end hinge), depth hangs down, and colour is the coherent-stack
%   response of every block at every depth. Colour saturates at 2 sigma
%   (per-pass jackknife) and washes toward neutral below, as on the maps.
%
%   Two curtains sets by default: the chain's 5 m reference, which carries
%   the shallow-firn term, and the 60 m reference, which shows the ice
%   column alone. Depth and across-track are exaggerated so the 250 m deep,
%   120-210 m spaced curtains read; the exaggerations are printed on the
%   axes.
%
%   opts.files {..} .mat names in out_dir (tidal_stack_nozc.mat,
%   tidal_stack_nozc_ref60.mat), opts.labels, opts.field ('a_t'),
%   opts.view [az el]. Writes EAGER_2022_tidal_stack_fence.png.

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here))); addpath(here);
def = struct('out_dir', vdef.figureDir(), ...
  'files', {{'tidal_stack_nozc.mat', 'tidal_stack_nozc_ref60.mat'}}, ...
  'labels', {{'Reference 5 m: firn and ice column', 'Reference 60 m: ice column alone'}}, ...
  'field', 'a_t', 'sd_field', 'a_t_sd', 'view', [-38 22], 'clim', 4, 'nsig', 2, ...
  'zmax', 320, 'xex', 1, 'yex', 4, 'zex', 8, 'rg', 'radargrams_nozc.mat');
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end, end

PAL.neg = [0.698 0.094 0.169]; PAL.mid = [0.941 0.937 0.925]; PAL.pos = [0.165 0.471 0.839];
ndiv = 256; half = ndiv/2;
dmap = [interp1([0 1], [PAL.neg; PAL.mid], linspace(0, 1, half)); interp1([0 1], [PAL.mid; PAL.pos], linspace(0, 1, ndiv - half))];
ps = projcrs(3031);

% common rotated frame from the first file's block positions
A = load(fullfile(opts.out_dir, opts.files{1})); OUT0 = A.OUT;
[u, org] = survey_frame(OUT0, ps);

RG = []; if exist(fullfile(opts.out_dir, opts.rg), 'file'), RG = load(fullfile(opts.out_dir, opts.rg)); end
np = numel(opts.files);
[h, GRL] = grl_figure(170, 75*np + 12); set(0, 'CurrentFigure', h);
ph = (1 - 0.10) / np;
for p = 1:np
  B = load(fullfile(opts.out_dir, opts.files{p})); OUT = B.OUT;
  ax = axes('parent', h, 'Position', [0.02 1 - p*ph + 0.01 0.84 ph - 0.03]); hold(ax, 'on');
  for k = 1:numel(OUT)
    O = OUT(k); [s, n] = rotate(O, ps, u, org);
    [s, o] = sort(s); n = n(o);
    se = [s(1) - (s(2)-s(1))/2; (s(1:end-1) + s(2:end))/2; s(end) + (s(end)-s(end-1))/2];
    ne = interp1(s, n, se, 'linear', 'extrap');
    zall = O.zsel(:); zi = zall <= opts.zmax; z = zall(zi);
    ze = [max(0, z(1) - (z(2)-z(1))/2); (z(1:end-1) + z(2:end))/2; z(end) + (z(end)-z(end-1))/2];
    a = 1e3*O.(opts.field)(zi, o); sd = 1e3*O.(opts.sd_field)(zi, o);
    rgb = wash(a, sd, dmap, opts.clim, opts.nsig, PAL.mid);
    rgb(repmat(isnan(a), 1, 1, 3)) = NaN;                                 % below the ice base: not drawn
    rgb = cat(1, rgb, rgb(end,:,:)); rgb = cat(2, rgb, rgb(:,end,:));        % vertex-padded for flat faces
    X = repmat(se.', numel(ze), 1); Y = repmat(ne.', numel(ze), 1); Z = repmat(-ze, 1, numel(se));
    surf(ax, X*opts.xex, Y*opts.yex, Z/1e3*opts.zex, rgb, 'FaceColor', 'flat', 'EdgeColor', 'none');
    plot3(ax, se*opts.xex, ne*opts.yex, zeros(size(se)), '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.8);   % the track on the surface
    nm = strrep(O.name, 'EAGER_2022_', '');
    if ~isempty(RG) && isfield(RG, nm)                                    % ice base, picked on the radargram
      G = RG.(nm); [xg, yg] = projfwd(ps, G.lat(:), G.lon(:));
      sg = ((xg - org(1))*u(1) + (yg - org(2))*u(2))/1e3; ng = ((xg - org(1))*(-u(2)) + (yg - org(2))*u(1))/1e3;
      hb = plot3(ax, sg*opts.xex, ng*opts.yex, -G.base(:)/1e3*opts.zex, '-', 'Color', [0.05 0.25 0.55], 'LineWidth', 1.6);
    end
    text(ax, se(1)*opts.xex - 0.15, ne(1)*opts.yex, 0.05, strrep(O.name, 'EAGER_2022_', ''), 'FontSize', 7, ...
      'HorizontalAlignment', 'right');
  end
  % hinge line and the 60 m reference level
  yl = ylim(ax);
  plot3(ax, [0 0], yl, [0 0], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.8);
  text(ax, 0, yl(2), 0.08, 'hinge', 'FontSize', 6.5, 'HorizontalAlignment', 'center');
  set(ax, 'DataAspectRatio', [1 1 1], 'FontSize', 7, 'Box', 'off', 'TickDir', 'out');
  view(ax, opts.view); grid(ax, 'on'); set(ax, 'GridAlpha', 0.12);
  zt = 0:50:300; set(ax, 'ZTick', -fliplr(zt)/1e3*opts.zex, 'ZTickLabel', arrayfun(@num2str, fliplr(zt), 'uni', 0));
  xt = -0.5:0.5:4.5; set(ax, 'XTick', xt*opts.xex, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), xt, 'uni', 0));
  yt = -0.2:0.2:0.2; set(ax, 'YTick', yt*opts.yex, 'YTickLabel', arrayfun(@(v) sprintf('%g', v), yt, 'uni', 0));
  xlabel(ax, 'Along survey (km)', 'FontSize', 7); ylabel(ax, sprintf('Across (km, \\times%d)', opts.yex), 'FontSize', 7);
  zlabel(ax, sprintf('Depth (m, \\times%d)', opts.zex), 'FontSize', 7);
  title(ax, sprintf('(%s)  %s', char('a' + p - 1), opts.labels{p}), 'FontWeight', 'normal', 'FontSize', 8);
  colormap(ax, dmap); clim(ax, [-opts.clim opts.clim]);
  if exist('hb', 'var'), legend(ax, hb, 'ice base', 'Location', 'southwest', 'FontSize', 7, 'Box', 'off'); end
end
cb = colorbar(ax, 'Position', [0.90 0.25 0.015 0.5]);
set(cb, 'FontSize', 7, 'TickDirection', 'out'); set(cb.Label, 'String', 'mm per m of tide', 'FontSize', 7);
annotation(h, 'textbox', [0.03 0 0.86 0.06], 'String', ...
  sprintf(['Each curtain is one survey line at its true position: response of every 500 m block at every depth, trend-controlled coherent stack. ' ...
   'Full colour at %d sigma (per-pass jackknife), washed toward neutral below. Across-track and depth exaggerated.'], opts.nsig), ...
  'EdgeColor', 'none', 'FontSize', 6.5, 'Color', [0.5 0.5 0.5], 'VerticalAlignment', 'bottom');
out_fn = fullfile(opts.out_dir, 'EAGER_2022_tidal_stack_fence.png');
print(h, out_fn, '-dpng', sprintf('-r%d', GRL.dpi)); close(h);
fprintf('Wrote %s\n', out_fn);
end

function [u, org] = survey_frame(OUT, ps)
X = []; Y = []; X0 = []; Y0 = [];
for k = 1:numel(OUT)
  [x, y] = projfwd(ps, OUT(k).lat(:), OUT(k).lon(:)); X = [X; x]; Y = [Y; y]; %#ok<AGROW>
  [~, i0] = min(abs(OUT(k).x_sea)); X0(end+1) = x(i0); Y0(end+1) = y(i0); %#ok<AGROW>
end
C = cov([X Y]); [V, D] = eig(C); [~, m] = max(diag(D)); u = V(:, m);
O = OUT(1); [~, i0] = min(O.x_sea); [~, i1] = max(O.x_sea);
[xa, ya] = projfwd(ps, O.lat([i0 i1]), O.lon([i0 i1]));
if u.' * [diff(xa); diff(ya)] < 0, u = -u; end      % +x points seaward
org = [mean(X0); mean(Y0)];
end

function [s, n] = rotate(O, ps, u, org)
[x, y] = projfwd(ps, O.lat(:), O.lon(:));
s = ((x - org(1))*u(1) + (y - org(2))*u(2))/1e3;
n = ((x - org(1))*(-u(2)) + (y - org(2))*u(1))/1e3;
end

function rgb = wash(a, sd, dmap, clim, nsig, mid)
ndiv = size(dmap, 1);
ci = max(1, min(ndiv, round((a + clim)/(2*clim)*(ndiv - 1)) + 1)); ci(~isfinite(ci)) = round(ndiv/2);
w = min(1, abs(a)./(nsig*sd)); w(~isfinite(w)) = 0;
rgb = zeros([size(a) 3]);
for c = 1:3
  col = dmap(:, c); rgb(:,:,c) = mid(c) + w.*(col(ci) - mid(c));
end
end
