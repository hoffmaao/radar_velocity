function tidal_fence_movie(opts)
%TIDAL_FENCE_MOVIE The shelf through a tidal cycle: radargrams on 3-D curtains,
%   moving with the tide, coloured by the tidal vertical strain rate we measure.
%
%   Each survey line is a vertical curtain at its true position (EPSG:3031
%   rotated so the survey runs along x, origin at the mean hinge) carrying
%   its own echo image (main pass of the surface-coupled build, log power,
%   radargrams_nozc.mat). Through one cycle of an ILLUSTRATIVE diurnal tide
%   with the survey's 1.1 m range:
%
%     - the curtains rise and fall by w(x)*eta(t), w the thin-plate hinge
%       profile (0 grounded, 1 floating, lambda = opts.lambda), so the hinge
%       stays put and the floating shelf heaves; the motion is exaggerated
%       (opts.vex) because 1 m is invisible against 250 m of ice,
%     - the colour over the echo image is the vertical strain RATE
%       d eps_zz/dt (x, z, t) = d a/d z (x, z) * d eta/dt, with a the in-phase
%       response of the phase-lag fit (a_q: tide, quadrature tide and trend
%       fitted together; tidal_stack_nozc.mat, 5 m reference, so the
%       shallow-firn term is included), smoothed over ~50 m in depth; its
%       opacity grows with the rate so the echoes show through weak strain.
%       The quadrature response is left out: the measured lag is 0.15 +/-
%       0.42 h (scripts/diagnostics/tidal_phase.m), so the strain is in phase
%       with the tide and its rate peaks at mid-tide, a quarter period
%       before high water, and vanishes at high and low water.
%       opts.colour = 'strain' draws eps_zz = d a/d z * eta instead,
%     - side panels track eta(t) and the tide rate that drives the colour.
%
%   The shallow-firn part of the strain cannot be told from a daily process
%   in this pass set (tide and hour of day correlate at R^2 = 0.99); the
%   caption says so.
%
%   Writes EAGER_2022_tidal_fence_movie.mp4 and a still where the colour
%   peaks (mid rising tide for the rate, high tide for the strain).

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here))); addpath(here);
def = struct('out_dir', vdef.figureDir(), 'stack', 'tidal_stack_nozc.mat', 'rg', 'radargrams_nozc.mat', ...
  'nframes', 72, 'fps', 12, 'period_h', 24.84, 'amp', 0.55, 'lambda', 1.3, ...
  'zex', 8, 'yex', 4, 'vex', 0.35, 'eclim', 10, 'rclim', 25, 'zmax', 320, 'view', [-38 22], ...
  'colour', 'rate', 'field', 'a_q');
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end, end
A = load(fullfile(opts.out_dir, opts.stack)); OUT = A.OUT;
RG = load(fullfile(opts.out_dir, opts.rg));
ps = projcrs(3031);

% ---- rotated frame (as the maps and the fence)
X = []; Y = []; X0 = []; Y0 = [];
for k = 1:numel(OUT)
  [x, y] = projfwd(ps, OUT(k).lat(:), OUT(k).lon(:)); X = [X; x]; Y = [Y; y]; %#ok<AGROW>
  [~, i0] = min(abs(OUT(k).x_sea)); X0(end+1) = x(i0); Y0(end+1) = y(i0); %#ok<AGROW>
end
[V, D] = eig(cov([X Y])); [~, m] = max(diag(D)); u = V(:, m);
[~, i0] = min(OUT(1).x_sea); [~, i1] = max(OUT(1).x_sea);
[xa, ya] = projfwd(ps, OUT(1).lat([i0 i1]), OUT(1).lon([i0 i1]));
if u.' * [diff(xa); diff(ya)] < 0, u = -u; end
org = [mean(X0); mean(Y0)];
rot = @(x, y) deal(((x - org(1))*u(1) + (y - org(2))*u(2))/1e3, ((x - org(1))*(-u(2)) + (y - org(2))*u(1))/1e3);

% ---- per line: radargram geometry, grey image, strain per metre of tide
PAL.neg = [0.698 0.094 0.169]; PAL.mid = [0.941 0.937 0.925]; PAL.pos = [0.165 0.471 0.839];
ndiv = 256; half = ndiv/2;
dmap = [interp1([0 1], [PAL.neg; PAL.mid], linspace(0, 1, half)); interp1([0 1], [PAL.mid; PAL.pos], linspace(0, 1, ndiv - half))];
Lc = struct('s', {}, 'n', {}, 'z', {}, 'grey', {}, 'eps', {}, 'w', {}, 'base', {});
for k = 1:numel(OUT)
  O = OUT(k); nm = strrep(O.name, 'EAGER_2022_', ''); R = RG.(nm);
  zk = R.z(:); keepz = zk <= opts.zmax; zk = zk(keepz);
  [xr, yr] = projfwd(ps, R.lat(:), R.lon(:)); [s, n] = rot(xr, yr);
  img = double(R.img(keepz, :));
  lo = prctile(img(:), 5); hi = prctile(img(:), 99.5);
  grey = 0.97 - 0.82*min(max((img - lo)/(hi - lo), 0), 1);              % strong echo = dark
  % strain per metre of tide: d a/dz on the stack's depth levels, smoothed
  [xb, yb] = projfwd(ps, O.lat(:), O.lon(:)); [sb, ~] = rot(xb, yb);
  zs = O.zsel(:); a = O.(opts.field);                                   % m per m of tide
  as = movmean(a, 5, 1, 'omitnan');
  e = [diff(as, 1, 1) ./ diff(zs); nan(1, size(a, 2))];                 % per m of tide
  zc = [ (zs(1:end-1) + zs(2:end))/2; zs(end)];
  [sbs, ob] = sort(sb); e = e(:, ob);
  % down to the base: the stack stops 15 m above each block's base (clear of
  % the basal echo) and the derivative loses a level more, so in each block
  % the deepest valid strain is held down the column, and the field is then
  % cut at the base picked on this radargram, column by column
  for b = 1:size(e, 2)
    v = find(isfinite(e(:, b)));
    if isempty(v), continue; end
    e(1:v(1)-1, b) = e(v(1), b); e(v(end)+1:end, b) = e(v(end), b);
    e(:, b) = fillmissing(e(:, b), 'linear');
  end
  zq = [0; zc; max(zk) + 1]; sq = [min(s) - 1; sbs(:); max(s) + 1];     % cover the whole curtain
  epad = e([1 1:end end], [1 1:end end]);                                % padded in depth and along track
  E = interp2(sq.', zq, epad, repmat(s.', numel(zk), 1), repmat(zk, 1, numel(s)), 'linear');
  E(zk > R.base(:).') = NaN;                                             % nothing below the ice base
  xs = s; w = zeros(size(xs)); fl = xs > 0; bx = xs(fl)/opts.lambda;
  w(fl) = 1 - exp(-bx).*(cos(bx) + sin(bx));                            % thin-plate hinge profile
  Lc(k) = struct('s', s(:).', 'n', n(:).', 'z', zk, 'grey', grey, 'eps', E, 'w', w(:).', 'base', R.base(:).');
end

% ---- tide
t = linspace(0, opts.period_h, opts.nframes + 1); t = t(1:end-1);
eta = opts.amp * sin(2*pi*t/opts.period_h);
deta = opts.amp * 2*pi/opts.period_h * cos(2*pi*t/opts.period_h);   % m per hour
rate = strcmp(opts.colour, 'rate');
if rate   % strain rate in 1e-6 per hour
  drive = deta * 1e6; clim_ = opts.rclim; cblab = 'Strain rate (10^{-6}/h)'; still_k = find(deta == max(deta), 1);
else      % strain in 1e-5
  drive = eta * 1e5;  clim_ = opts.eclim; cblab = 'Vertical strain (10^{-5})'; still_k = find(eta == max(eta), 1);
end

% ---- frames
vw = VideoWriter(fullfile(opts.out_dir, 'EAGER_2022_tidal_fence_movie.mp4'), 'MPEG-4');
vw.FrameRate = opts.fps; vw.Quality = 95; open(vw);
for f = 1:opts.nframes
  h = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', 'Position', [0 0 1600 900]);
  if exist('theme', 'file'), theme(h, 'light'); end
  ax = axes('parent', h, 'Position', [0.02 0.10 0.70 0.84]); hold(ax, 'on');
  for k = 1:numel(Lc)
    C = Lc(k); ep = C.eps * drive(f);                                   % colour units
    ep(~isfinite(ep)) = 0;     % no estimate: plain radargram (min/max ignore NaN and would saturate it)
    ci = max(1, min(ndiv, round((ep + clim_)/(2*clim_)*(ndiv - 1)) + 1)); ci(~isfinite(ci)) = half;
    al = 0.85*min(1, abs(ep)/clim_); al(~isfinite(al)) = 0;
    rgb = zeros([size(C.grey) 3]);
    for c = 1:3, col = dmap(:, c); rgb(:,:,c) = (1 - al).*C.grey + al.*col(ci); end
    rgb(repmat(~isfinite(C.grey), 1, 1, 3)) = 1;
    Zs = -C.z/1e3*opts.zex + C.w*eta(f)*opts.vex;
    surf(ax, repmat(C.s, numel(C.z), 1), repmat(C.n*opts.yex, numel(C.z), 1), Zs, rgb, ...
      'FaceColor', 'texturemap', 'EdgeColor', 'none');
    plot3(ax, C.s, C.n*opts.yex, C.w*eta(f)*opts.vex, '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 1);
    hbase = plot3(ax, C.s, C.n*opts.yex, -C.base/1e3*opts.zex + C.w*eta(f)*opts.vex, '-', ...
      'Color', [0.05 0.25 0.55], 'LineWidth', 2.2);                       % ice base, heaving with the shelf
  end
  set(ax, 'DataAspectRatio', [1 1 1], 'XLim', [-0.8 4.8], 'YLim', [-0.35 0.35]*opts.yex, ...
    'ZLim', [-opts.zmax/1e3*opts.zex - 0.1, opts.amp*opts.vex + 0.3], 'FontSize', 10, 'TickDir', 'out');
  view(ax, opts.view); grid(ax, 'on'); set(ax, 'GridAlpha', 0.12);
  zt = 0:50:300; set(ax, 'ZTick', -fliplr(zt)/1e3*opts.zex, 'ZTickLabel', arrayfun(@num2str, fliplr(zt), 'uni', 0));
  set(ax, 'YTick', (-0.2:0.2:0.2)*opts.yex, 'YTickLabel', {'-0.2', '0', '0.2'});
  xlabel(ax, 'Along track (km)'); ylabel(ax, sprintf('Across track (km, \\times%d)', opts.yex)); zlabel(ax, sprintf('Ice depth (m, \\times%d)', opts.zex));
  colormap(ax, dmap); clim(ax, [-clim_ clim_]);
  legend(ax, hbase, 'ice base (picked)', 'Location', 'southwest', 'FontSize', 10, 'Box', 'off');
  cb = colorbar(ax, 'Position', [0.715 0.55 0.012 0.33]); cb.Label.String = cblab; cb.FontSize = 10;
  % tide panel
  a2 = axes('parent', h, 'Position', [0.815 0.58 0.165 0.30]); hold(a2, 'on');
  plot(a2, t, eta, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.4);
  plot(a2, t(f), eta(f), 'o', 'MarkerFaceColor', PAL.neg, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
  set(a2, 'XLim', [0 opts.period_h], 'YLim', [-0.7 0.7], 'FontSize', 9, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:6:24);
  grid(a2, 'on'); xlabel(a2, 'Elapsed time (h)'); ylabel(a2, 'Tide height (m)'); title(a2, 'Illustrative diurnal tide', 'FontWeight', 'normal', 'FontSize', 10);
  a3 = axes('parent', h, 'Position', [0.815 0.16 0.165 0.30]); hold(a3, 'on');
  plot(a3, t, deta*100, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.4);
  plot(a3, t(f), deta(f)*100, 'o', 'MarkerFaceColor', PAL.pos, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
  set(a3, 'XLim', [0 opts.period_h], 'FontSize', 9, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:6:24);
  grid(a3, 'on'); xlabel(a3, 'Elapsed time (h)'); ylabel(a3, 'Tide rate (cm/h)');
  title(a3, 'Strain rate follows the tide rate', 'FontWeight', 'normal', 'FontSize', 10);
  img = print(h, '-RGBImage', '-r96');            % off-screen: getframe needs a live display and can block when it sleeps
  img = img(1:2*floor(end/2), 1:2*floor(end/2), :);
  writeVideo(vw, img);
  if f == still_k, imwrite(img, fullfile(opts.out_dir, 'EAGER_2022_tidal_fence_still.png')); end
  close(h);
end
close(vw);
fprintf('Wrote %s (%d frames)\n', fullfile(opts.out_dir, 'EAGER_2022_tidal_fence_movie.mp4'), opts.nframes);
end
