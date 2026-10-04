function OUT = tidal_depth_profile(opts)
%TIDAL_DEPTH_PROFILE The tidal admittance of the column as a FUNCTION OF DEPTH.
%
%   The single-depth tidal analysis reads the column change at one depth
%   (100 m) per block. The products carry the whole profile dh(z) at 0.3 m
%   spacing, and a Legendre expansion in depth is already this project's
%   representation of a column profile (vdef.invertStrainRate, order 2 =
%   strain linear in depth). This driver takes that expansion to order
%   opts.order (default 3) on every pair's dh(z), so the tidal response
%   may be nonlinear in depth, and pushes it through the same two steps
%   the single-depth admittance takes:
%
%     per pair, per block   dh(z) = sum_k c_k P_k(x),  x = 2 z/H - 1
%     network inversion     every coefficient row over all pairs
%     joint tide fit        c_k = a_k + b_k t + adm_k tide per row
%                           (CATS2008 tide and the pass gate, pass_tide)
%
%   The tidal admittance of the column is then c(z) = sum_k adm_k P_k(x).
%
%   THE PLATE SHAPE IS THREE LEGENDRE TERMS, EXACTLY. Thin-plate bending
%   predicts dh(z) = A (z_n z - z^2/2), and that polynomial is finite in
%   this basis: with x = 2z/H - 1,
%
%       m_0 = A (z_n H/2 - H^2/6)
%       m_1 = A (z_n H/2 - H^2/4)
%       m_2 = -A H^2/12
%       m_k = 0   for k >= 3
%
%   Three consequences, and they are why the expansion is the right
%   representation here rather than a convenience:
%     * the AMPLITUDE is carried by P2 alone, A = -12 adm_2 / H^2;
%     * the NEUTRAL PLANE is a ratio of two fitted numbers,
%       z_n = H/2 - (H/6)(adm_1/adm_2), so it is measurable BLOCK BY
%       BLOCK rather than only as a line average;
%     * P3 and above are ZERO for any plate, whatever its thickness or
%       neutral plane, so the order-3 admittance is a direct test of
%       departure from thin-plate bending - it cannot be absorbed by
%       re-fitting A or z_n.
%
%   FIT QUALITY IS MEASURED AGAINST A BASIS-FREE PROFILE. The same chain
%   is run a second time with no expansion at all: dh interpolated onto a
%   depth grid, each depth network-inverted and tide-fitted on its own.
%   That pointwise admittance a_pt(z) is what the Legendre curve has to
%   pass through, and the normalised residual (c - a_pt)/sigma_pt says
%   whether the polynomial is representing the column or smoothing it
%   away. The same comparison repeated over orders 1..opts.order_max is
%   the order sweep: where it stops improving is the order the data
%   support.
%
%   opts, all optional:
%     .products, .net_dir, .mp_dir   as elastic_modulus
%     .order       Legendre order for the primary fit (default 3)
%     .order_max   highest order in the sweep (default 4)
%     .z_range     fitted depth range [m] (default [20 250])
%     .norm_depth  H for the basis (default 250, the products' own)
%     .z_step      pointwise depth grid step [m] (default 5)
%     .zn_grid     neutral-plane trial depths [m] (default 20:2:400)
%     .zmax_bend   blocks landward of this x_sea carry the bending
%                  analysis [m] (default 2500; beyond it the curvature,
%                  and with it the signal, is gone)
%     .fit_file    flexure_fit_cats.mat, for the block frame and h(x)
%     .out_dir     where the figure goes
%
%   READ THIS BEFORE BELIEVING A BENDING AMPLITUDE. Every quality test in
%   this driver - along-track coherence, the plate shape, a tight z_n
%   interval, a high formal significance, a clean order sweep - is passed
%   MORE convincingly by EAGER_2022, the uncalibrated build, than by
%   EAGER_2022_GL1, which is the SAME LEG processed correctly. Run on the
%   two of them (4 Sep 2026) the uncalibrated build returns a median
%   column response of 6.41 mm/m against 1.65, a bending amplitude up to
%   15 sigma against 6.5, and a neutral plane of 148-176 m marching
%   monotonically along the line. The difference between the two builds
%   is the project's systematic floor, and it EXCEEDS the calibrated
%   build's whole signal. The formal sigmas here come from the tide fit
%   and do not see that; nothing in this file does. So a bending
%   amplitude from this driver is a formal statement about one build,
%   not evidence of flexure, until the same quantity has been shown to
%   reproduce between independent builds of the same ice.
%
%   Prints, and returns in OUT, the coefficient admittances, the
%   per-block neutral plane and bending amplitude with their intervals,
%   and the fit-quality residuals.

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));                       % +vdef
addpath(here);                                             % pass_tide
addpath(fullfile(fileparts(here), 'figures'));             % grl_figure

def = struct( ...
  'products', {{'EAGER_2022_GL1','EAGER_2022_GL2','EAGER_2022_GL3','EAGER_2022_GL4'}}, ...
  'net_dir', '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel_net', ...
  'mp_dir',  '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass', ...
  'order', 3, 'order_max', 4, 'z_range', [20 250], 'norm_depth', 250, ...
  'z_step', 5, 'zn_grid', 20:2:400, 'zmax_bend', 2500, ...
  'max_baseline', 10, 'min_pairs', 10, ...
  'fit_file', vdef.figureDir('flexure_fit_cats.mat'), ...
  'out_dir', vdef.figureDir());
fn = fieldnames(def);
for i = 1:numel(fn)
  if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end
end
K = opts.order; Kmax = max(opts.order_max, K); H = opts.norm_depth;
C = vdef.constants();

FIT = [];
if exist(opts.fit_file, 'file')
  S = load(opts.fit_file); f = fieldnames(S); FIT = S.(f{1});
end

zpt = (opts.z_range(1):opts.z_step:opts.z_range(2)).';    % pointwise grid
nz  = numel(zpt);
Pz  = vdef.legendreBasis(2*zpt/H - 1, Kmax);

OUT = [];
for n = 1:numel(opts.products)
  pn = opts.products{n}; nm = strrep(pn, 'EAGER_2022_', '');
  f = dir(fullfile(opts.net_dir, [pn '_vvel_*.mat']));
  keep = ~cellfun('isempty', regexp({f.name}, ...
    ['^' regexptranslate('escape',pn) '_vvel_\d+_\d+\.mat$'], 'once'));
  f = f(keep);
  if isempty(f), fprintf('%s: no products\n', nm); continue; end

  [tide, pass_ok, tinfo] = pass_tide(pn, opts.mp_dir);
  Np = numel(tide); tday = (tinfo.tmid - min(tinfo.tmid))/86400;

  %% One pass over the pair files: every order, and the pointwise grid
  P = []; W = []; along = []; Nblk = 0;
  Dk = cell(1, Kmax); Dpt = [];
  sopts = struct('order', K, 'fit_top_depth', opts.z_range(1), ...
    'fit_bot_depth', opts.z_range(2), 'norm_depth', H, 'reg', 0, ...
    'bins_per_look', 5, 'min_samples', 50);
  for q = 1:numel(f)
    tok = regexp(f(q).name, ...
      ['^' regexptranslate('escape',pn) '_vvel_(\d+)_(\d+)\.mat$'], 'tokens','once');
    o = load(fullfile(opts.net_dir, f(q).name));
    if ~vdef.pairAligned(o), continue; end
    if max(abs(o.baseline_y)) > opts.max_baseline, continue; end
    if Nblk == 0
      Nblk = size(o.dh_blk, 2); along = o.Along_track(:);
      if isfield(o,'param_vvel') && isfield(o.param_vvel.vvel,'bins_per_look')
        sopts.bins_per_look = o.param_vvel.vvel.bins_per_look;
      end
    end
    % dh(z) stands in for v(d): the fit is identical, the unit is metres
    % of column change rather than metres per year.
    V = struct('v', o.dh_blk, 'depth', o.depth_blk, ...
               'v_scatter', o.dtau_scatter_blk * C.c/2 ./ o.n_local);
    blk = struct('coh', o.coh_blk);
    cf = cell(1, Kmax); okall = true;
    for kk = 1:Kmax
      so = sopts; so.order = kk;
      Sk = vdef.invertStrainRate(V, blk, so);
      cf{kk} = Sk.coef;
      if all(~isfinite(Sk.coef(:))), okall = false; end
    end
    if ~okall, continue; end
    % Basis-free: the same profile read straight off the product
    pt = nan(nz, Nblk);
    for b = 1:Nblk
      d = o.depth_blk(:,b); y = o.dh_blk(:,b);
      g = isfinite(d) & isfinite(y);
      if nnz(g) < 50, continue; end
      pt(:,b) = interp1(d(g), y(g), zpt, 'linear', NaN);
    end
    for kk = 1:Kmax, Dk{kk}(:, end+1) = cf{kk}(:); end
    Dpt(:, end+1) = pt(:); %#ok<AGROW>
    P(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
    W(end+1) = max(mean(o.coh_blk(:), 'omitnan'), 1e-3);   %#ok<AGROW>
  end
  if size(P,1) < opts.min_pairs
    fprintf('%s: only %d usable pairs\n', nm, size(P,1)); continue;
  end
  nopts = struct('n_sigma', 3, 'weights', W, 'n_epoch', Np);

  %% Admittance of every coefficient, at every order
  adm_k = cell(1, Kmax); ads_k = cell(1, Kmax);
  for kk = 1:Kmax
    N = vdef.invertNetwork(P, Dk{kk}, nopts);
    A = vdef.fitTideAdmittance(N.x, tday, tide);
    adm_k{kk} = reshape(A.admittance(:), kk+1, Nblk);
    ads_k{kk} = reshape(A.admittance_std(:), kk+1, Nblk);
  end
  adm = adm_k{K}; ads = ads_k{K};

  %% The basis-free admittance, through the identical two steps
  Npt = vdef.invertNetwork(P, Dpt, nopts);
  Apt = vdef.fitTideAdmittance(Npt.x, tday, tide);
  a_pt = reshape(Apt.admittance(:), nz, Nblk);
  s_pt = reshape(Apt.admittance_std(:), nz, Nblk);

  %% Block frame and thickness, from the flexure fit
  x_sea = along - along(1); h_half = nan(Nblk,1);
  if ~isempty(FIT)
    i = find(strcmp({FIT.lines.name}, pn), 1);
    if ~isempty(i)
      Li = FIT.lines(i); Ri = FIT.fits{1, i};
      x_sea = Li.x_flip_sign * (along - Li.x_flip_ref);
      if isnumeric(Ri.h_spec) && ~isscalar(Ri.h_spec)
        hb = interp1(Ri.h_spec(:,1), Ri.h_spec(:,2), ...
               min(max(x_sea, Ri.h_spec(1,1)), Ri.h_spec(end,1)));
      else
        hb = repmat(Ri.h_spec, Nblk, 1);
      end
      h_half = hb(:)/2;
    end
  end

  %% Fit quality: the Legendre curve against the basis-free profile
  % WHICH BLOCKS COUNT. A block whose network solution is broken comes
  % back with coefficient admittances in the hundreds or thousands of
  % mm/m and sigmas to match, and it is neither a measurement nor a test
  % of the basis. Those are gated out on the coefficients themselves, not
  % on the pointwise error, which is what an earlier version did and it
  % let them through. The statistic over the survivors is then a MEDIAN
  % over blocks: even after gating, one bad block in ten would dominate a
  % mean and turn an order sweep into noise.
  MAXADM = 50e-3;                       % mm/m ceiling on a sane coefficient
  meas = isfinite(x_sea) & median(s_pt, 1, 'omitnan').' < 5e-3;
  for kk = 1:Kmax
    meas = meas & all(isfinite(adm_k{kk}), 1).' & all(ads_k{kk} > 0, 1).' & ...
           all(abs(adm_k{kk}) < MAXADM, 1).';
  end
  good = find(meas);
  chi_ord = nan(1, Kmax); res_n = nan(nz, Nblk); frac1 = nan(1, Kmax);
  for kk = 1:Kmax
    ck = Pz(:,1:kk+1) * adm_k{kk};
    r  = (ck - a_pt) ./ s_pt;
    per_blk = mean(r(:, good).^2, 1, 'omitnan');
    chi_ord(kk) = median(per_blk, 'omitnan');
    frac1(kk) = mean(abs(r(:, good)) < 1, 'all', 'omitnan');
    if kk == K, res_n = r; end
  end
  c = Pz(:,1:K+1) * adm;
  c_std = sqrt(Pz(:,1:K+1).^2 * ads.^2);

  %% Plate shape, block by block
  % The plate has exactly three Legendre terms (see the header), so the
  % amplitude comes from P2 alone and the neutral plane from the ratio
  % adm_1/adm_2. The interval is profiled rather than propagated: a ratio
  % whose denominator is not significant has no upper bound, and the
  % profile says so instead of returning a small number.
  m_of = @(zn) [zn*H/2 - H^2/6; zn*H/2 - H^2/4; -H^2/12; zeros(K-2,1)];
  kk = 2:K+1;                       % P1..PK; P0 absorbs surface reference
  zn_grid = opts.zn_grid;
  zn = nan(Nblk,1); zn_lo = nan(Nblk,1); zn_hi = nan(Nblk,1);
  Ab = nan(Nblk,1); Ab_sig = nan(Nblk,1); chi_b = nan(Nblk,1);
  for b = 1:Nblk
    if ~all(isfinite(adm(:,b))) || ~all(ads(:,b) > 0), continue; end
    w = 1 ./ ads(kk,b).^2;
    J = nan(size(zn_grid)); Ag = nan(size(zn_grid));
    for iz = 1:numel(zn_grid)
      mf = m_of(zn_grid(iz)); m = mf(kk);
      Ag(iz) = sum(w .* m .* adm(kk,b)) / sum(w .* m.^2);
      J(iz)  = sum(w .* (adm(kk,b) - Ag(iz)*m).^2);
    end
    [Jm, iz0] = min(J);
    zn(b) = zn_grid(iz0); Ab(b) = Ag(iz0);
    % delta-misfit of one in the INPUT sigmas: the tide fit's own errors,
    % not rescaled, because a per-block fit has one degree of freedom and
    % a variance estimated from it would be noise.
    ins = J <= Jm + 1;
    if iz0 > 1 && iz0 < numel(zn_grid)
      zn_lo(b) = zn_grid(find(ins,1,'first')); zn_hi(b) = zn_grid(find(ins,1,'last'));
      if zn_lo(b) == zn_grid(1),   zn_lo(b) = NaN; end
      if zn_hi(b) == zn_grid(end), zn_hi(b) = NaN; end
    end
    mf0 = m_of(zn(b)); m0 = mf0(kk);
    Ab_sig(b) = Ab(b) * sqrt(sum(w .* m0.^2));
    chi_b(b)  = Jm / max(numel(kk) - 1, 1);
  end

  %% Print
  % Every block the chain actually measured, along the WHOLE line. An
  % earlier version cut this at opts.zmax_bend on the assumption that
  % curvature dies seaward of it; the bending amplitude is significant
  % out to 4.5 km, so the cut was hiding measured blocks. Significance
  % now gates the reading, not distance.
  use = find(meas & isfinite(zn));
  fprintf('\n===== %s: %d pairs, %d passes, %d blocks, order %d over %d-%d m =====\n', ...
    nm, size(P,1), nnz(pass_ok), Nblk, K, opts.z_range);
  fprintf('%6s %6s %9s %9s %9s %8s %8s %7s\n', ...
    'x km','h/2','adm P1','adm P2','adm P3','z_n m','A/sig','chi2');
  for b = 1:Nblk
    fprintf('%6.2f %6.0f %+5.2f/%4.2f %+5.2f/%4.2f %+5.2f/%4.2f %8s %8.1f %7.2f\n', ...
      x_sea(b)/1e3, h_half(b), 1e3*adm(2,b), 1e3*ads(2,b), 1e3*adm(3,b), 1e3*ads(3,b), ...
      1e3*adm(min(4,K+1),b), 1e3*ads(min(4,K+1),b), ...
      zn_str(zn(b), zn_lo(b), zn_hi(b)), Ab_sig(b), chi_b(b));
  end
  fprintf('(admittances mm per m of tide on P_k(2z/H-1), H = %.0f m; z_n from adm1/adm2)\n', H);
  fprintf('fit quality over %d measured blocks, median ((fit - pointwise)/sigma)^2: %s\n', ...
    numel(good), sprintf('K%d %.3f  ', [1:Kmax; chi_ord]));
  fprintf('                                fraction of depths within 1 sigma: %s\n', ...
    sprintf('K%d %.2f  ', [1:Kmax; frac1]));
  sig = find(abs(Ab_sig) >= 2 & isfinite(zn));
  fprintf('blocks with bending over 2 sigma: %d (%s km); z_n there %s m\n', ...
    numel(sig), strtrim(sprintf('%.1f ', x_sea(sig)/1e3)), strtrim(sprintf('%.0f ', zn(sig))));

  O = struct('name', pn, 'x_sea', x_sea, 'h_half', h_half, 'zpt', zpt, ...
    'a_pt', a_pt, 's_pt', s_pt, 'c', c, 'c_std', c_std, 'res_n', res_n, ...
    'adm', adm, 'adm_std', ads, 'chi_ord', chi_ord, 'frac1', frac1, 'zn', zn, 'zn_lo', zn_lo, ...
    'zn_hi', zn_hi, 'A', Ab, 'A_sig', Ab_sig, 'chi_b', chi_b, 'use', use, ...
    'good', good, 'n_pair', size(P,1));
  if isempty(OUT), OUT = O; else, OUT(end+1) = O; end %#ok<AGROW>
end
if isempty(OUT), return; end

%% Figure
nL = numel(OUT);
cols = [0.11 0.42 0.69; 0.89 0.47 0.10; 0.20 0.60 0.25; 0.75 0.20 0.30];
mks  = {'o','s','d','v'};
[hf, GRL] = grl_figure(170, 125);
set(0, 'CurrentFigure', hf);

% (a) the fit, on the block with the strongest bending anywhere
best_i = 1; best_b = []; best_s = 0;
for i = 1:nL
  [v, b] = max(abs(OUT(i).A_sig));
  if isfinite(v) && v > best_s, best_s = v; best_i = i; best_b = b; end
end
ax = axes('parent', hf, 'Position', [0.07 0.58 0.24 0.33]); hold(ax,'on');
O = OUT(best_i); b = best_b;
plot(ax, [0 0], [0 O.zpt(end)], '-', 'Color', [0.8 0.8 0.8]);
for j = 1:numel(O.zpt)
  plot(ax, 1e3*(O.a_pt(j,b) + [-1 1]*O.s_pt(j,b)), [1 1]*O.zpt(j), '-', ...
    'Color', [0.65 0.65 0.65], 'LineWidth', 0.6);
end
plot(ax, 1e3*O.a_pt(:,b), O.zpt, '.', 'Color', [0.35 0.35 0.35], 'MarkerSize', 7);
plot(ax, 1e3*O.c(:,b), O.zpt, '-', 'Color', cols(best_i,:), 'LineWidth', 1.8);
gp = O.zn(b)*O.zpt - O.zpt.^2/2;
plot(ax, 1e3*O.A(b)*gp, O.zpt, '--', 'Color', cols(best_i,:), 'LineWidth', 1.2);
plot(ax, get(ax,'XLim'), [1 1]*O.zn(b), 'k-', 'LineWidth', 0.8);
set(ax,'YDir','reverse'); grid(ax,'on'); box(ax,'on'); ylim(ax,[0 O.zpt(end)]);
ylabel(ax,'Depth (m)'); xlabel(ax,'Admittance (mm per m)');
title(ax, {sprintf('(a) fit at %s, %.1f km', strrep(O.name,'EAGER_2022_',''), O.x_sea(b)/1e3), ...
  'grey: pointwise; line: order 3; dashed: plate'}, 'FontWeight','normal','FontSize',7);

% (b) normalised residual against the basis-free profile
ax = axes('parent', hf, 'Position', [0.40 0.58 0.24 0.33]); hold(ax,'on');
for i = 1:nL
  O = OUT(i);
  R = O.res_n(:, O.use);
  if isempty(R), continue; end
  q = quantile(R, [0.25 0.5 0.75], 2);
  fill(ax, [q(:,1); flipud(q(:,3))], [O.zpt; flipud(O.zpt)], cols(i,:), ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none');
  plot(ax, q(:,2), O.zpt, '-', 'Color', cols(i,:), 'LineWidth', 1.4);
end
plot(ax, [0 0], [0 zpt(end)], 'k-', 'LineWidth', 0.8);
plot(ax, [-1 -1], [0 zpt(end)], 'k--'); plot(ax, [1 1], [0 zpt(end)], 'k--');
set(ax,'YDir','reverse'); grid(ax,'on'); box(ax,'on');
xlim(ax,[-3 3]); ylim(ax,[0 zpt(end)]);
xlabel(ax,'(order 3 - pointwise) / \sigma'); set(ax,'YTickLabel',[]);
title(ax, {'(b) what the basis discards', 'median and quartiles; dashed: one sigma'}, ...
  'FontWeight','normal','FontSize',7);

% (c) order sweep
ax = axes('parent', hf, 'Position', [0.73 0.58 0.24 0.33]); hold(ax,'on');
for i = 1:nL
  plot(ax, 1:numel(OUT(i).chi_ord), OUT(i).chi_ord, '-', 'Color', cols(i,:), ...
    'LineWidth', 1.4, 'Marker', mks{i}, 'MarkerFaceColor', cols(i,:), 'MarkerSize', 4);
end
plot(ax, [1 numel(OUT(1).chi_ord)], [1 1], 'k--');
set(ax,'YScale','log'); grid(ax,'on'); box(ax,'on');
xlim(ax,[0.8 numel(OUT(1).chi_ord)+0.2]); set(ax,'XTick',1:numel(OUT(1).chi_ord));
xlabel(ax,'Legendre order'); ylabel(ax,'median ((fit - pointwise)/\sigma)^2');
title(ax, {'(c) order the data support', 'over measured blocks; flat = nothing left'}, ...
  'FontWeight','normal','FontSize',7);
lg = legend(ax, arrayfun(@(o) strrep(o.name,'EAGER_2022_',''), OUT, 'uni', 0), ...
  'Location','northeast','FontSize',6); set(lg,'Box','off');

% (d) neutral plane ALONG THE LINE
ax = axes('parent', hf, 'Position', [0.07 0.12 0.40 0.33]); hold(ax,'on');
for i = 1:nL
  O = OUT(i);
  plot(ax, O.x_sea/1e3, O.h_half, ':', 'Color', cols(i,:), 'LineWidth', 1);
  for bb = O.use.'
    sigb = abs(O.A_sig(bb)) >= 2;
    % A bar is drawn only where BOTH bounds close. A one-sided or
    % unbounded z_n is a block whose P2 admittance is not significant,
    % and drawing it as a full-height bar says "measured, badly" when
    % the truth is "not measured"; it gets an open marker instead.
    if sigb && isfinite(O.zn_lo(bb)) && isfinite(O.zn_hi(bb))
      plot(ax, [1 1]*O.x_sea(bb)/1e3, [O.zn_lo(bb) O.zn_hi(bb)], '-', ...
        'Color', cols(i,:), 'LineWidth', 1.1);
      plot(ax, O.x_sea(bb)/1e3, O.zn(bb), mks{i}, 'Color', cols(i,:), ...
        'MarkerFaceColor', cols(i,:), 'MarkerSize', 5);
    else
      plot(ax, O.x_sea(bb)/1e3, O.zn(bb), mks{i}, 'Color', 0.5*cols(i,:)+0.5, ...
        'MarkerSize', 4);
    end
  end
end
set(ax,'YDir','reverse'); grid(ax,'on'); box(ax,'on');
ylim(ax,[0 300]); xlabel(ax,'Seaward distance (km)'); ylabel(ax,'z_n (m)');
title(ax, {'(d) neutral plane along the line', ...
  'filled with bar: bending over 2\sigma and z_n bounded; dotted: h(x)/2'}, ...
  'FontWeight','normal','FontSize',7);

% (e) bending amplitude and the order-3 term along the line
ax = axes('parent', hf, 'Position', [0.57 0.12 0.40 0.33]); hold(ax,'on');
plot(ax, [-1 6], [0 0], '-', 'Color', [0.8 0.8 0.8]);
plot(ax, [-1 6], [2 2], 'k--'); plot(ax, [-1 6], [-2 -2], 'k--');
for i = 1:nL
  O = OUT(i);
  gd = isfinite(O.A_sig);
  plot(ax, O.x_sea(gd)/1e3, O.A_sig(gd), '-', 'Color', cols(i,:), 'LineWidth', 1.3, ...
    'Marker', mks{i}, 'MarkerFaceColor', cols(i,:), 'MarkerSize', 4);
  t3 = O.adm(min(4,K+1),:).' ./ O.adm_std(min(4,K+1),:).';
  gd3 = isfinite(t3);
  plot(ax, O.x_sea(gd3)/1e3, t3(gd3), ':', 'Color', cols(i,:), 'LineWidth', 1.1);
end
grid(ax,'on'); box(ax,'on'); ylim(ax,[-8 8]);
xl = [min(arrayfun(@(o) min(o.x_sea), OUT)) max(arrayfun(@(o) max(o.x_sea), OUT))]/1e3;
xlim(ax, xl + [-0.2 0.2]);
xlabel(ax,'Seaward distance (km)'); ylabel(ax,'Significance (\sigma)');
title(ax, {'(e) bending amplitude (solid) and the P3 term (dotted)', ...
  'P3 is zero for any plate, so it tests the shape'}, 'FontWeight','normal','FontSize',7);

if ~exist(opts.out_dir, 'dir'), mkdir(opts.out_dir); end
out_fn = fullfile(opts.out_dir, 'EAGER_2022_tidal_depth_profile.png');
print(hf, out_fn, '-dpng', sprintf('-r%d', GRL.dpi));
fprintf('\nWrote %s\n', out_fn);
save(fullfile(opts.out_dir, 'tidal_depth_profile.mat'), 'OUT');
end

%% ========================================================================
function s = zn_str(zn, lo, hi)
if ~isfinite(zn), s = '-'; return; end
if isfinite(lo) && isfinite(hi)
  s = sprintf('%.0f[%.0f,%.0f]', zn, lo, hi);
elseif isfinite(lo)
  s = sprintf('%.0f[%.0f,-]', zn, lo);
elseif isfinite(hi)
  s = sprintf('%.0f[-,%.0f]', zn, hi);
else
  s = sprintf('%.0f[unb]', zn);
end
end
