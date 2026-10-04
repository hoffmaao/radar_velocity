function R = tidalStack(I, W, dtide, kz, opts)
%TIDALSTACK Coherent all-pairs tidal matched filter on complex block fields.
%   R = TIDALSTACK(I, W, dtide, kz, opts) estimates the tidal column
%   response a(z, x) by counter-rotating every pair's surface-referenced
%   complex block field by the phase the tidal model predicts,
%
%       Phi_p(z) = kz(z) * a(z, x) * dtide_p,   kz = phase_sign*4*pi*fc*n(z)/c
%
%   and summing over pairs. The trial a at which the pairs add coherently
%   is the estimate; the normalised peak F says how tide-locked the phase
%   is (1 = every pair phase-locked). Every pair is interfered against
%   every other through the shared model, and nothing is unwrapped, so
%   this estimator is independent of the unwrap -> dh -> regression chain.
%
%   WHAT ELSE IS IN THE PHASE. The tide is not the only thing that differs
%   between two passes, and whatever else does and happens to correlate
%   with the tide difference across the pairs lands in a. Two such terms
%   are handled, each as a second regressor in a joint 2-D scan:
%
%   opts.dt - the SECULAR TREND. Firn compaction and dynamic strain move
%   every layer relative to the surface at a steady rate, so a pair's
%   phase carries kz*v*dt_p as well as the tidal term. The tide and the
%   pair interval are not orthogonal over a few days of passes, and
%   without the trend the stack attributes v*cov(dt,dtide)/var(dtide) to
%   the tide. With opts.dt the scan fits Phi = kz*(a*dtide + v*dt)
%   jointly: R.a_t is the trend-controlled response, R.v_t the rate, and
%   R.r_tt = corr(dt, dtide) says how separable the two were.
%
%   opts.dres - the MIS-REGISTRATION ARTEFACT of a product built with the
%   z-motion compensation (vdef.zmotionApplied): each pair's measured
%   residual misalignment per block (ns), fitted as Phi = kz*(a*dtide +
%   g*dres). Where dres is collinear with dtide (|r| near 1) the surface is
%   a ridge and a2 is not identified; R.rcol reports that per block.
%
%   opts.dq - the QUADRATURE TIDE. An elastic column strains in phase with
%   the tide, so its strain RATE leads the tide by a quarter period. A
%   response that lags the tide (viscous, or delayed by anything) carries
%   a component in phase with the tide's rate instead, which the model
%   above cannot represent. With opts.dq - each pair's difference of the
%   quadrature tide, the tide rate scaled to metres (T/(2*pi) * deta/dt for
%   the diurnal period T) - the response is fitted as Phi = kz*(a*dtide +
%   b*dq [+ v*dt]): R.a_q, R.b_q (and R.v_q with opts.dt). b = 0 is the
%   elastic quarter-period; for a sinusoidal tide a response lagging by
%   tau gives b/a = -tan(2*pi*tau/T). Three parameters make the grid scan
%   too large, so this fit starts at the grid peak of the in-phase model
%   (b = 0) and climbs the same coherent sum continuously (Gauss-Newton on
%   the pair phases, with the common phase free; every step is accepted
%   only if F rises). R.r_qt = corr(dq, dtide) says how separable the two
%   tides were. The phases are a fraction of a radian, so the coherent sum
%   has one peak near the start.
%
%   THE ERROR BAR. Pairs are not independent samples: a pass's own error -
%   its registration, its surface, its day's conditions - enters every pair
%   that pass is in, and a bootstrap over pairs treats those as fresh draws
%   and understates the scatter. With opts.pairs the standard error is the
%   delete-one-PASS jackknife (each pass left out with all of its pairs),
%   returned as R.a_sd_jk and, because it is the honest one, as R.a_sd;
%   the pair bootstrap is kept as R.a_sd_boot. Without opts.pairs, R.a_sd
%   is the bootstrap.
%
%   I, W    nz x nb x np complex block means and their weights
%           (vdef.complexBlocks), one slab per pair
%   dtide   np-vector, tide difference of each pair [m]
%   kz      nz-vector, rad per metre of column change at each depth
%   opts    .agrid     trial responses [m per m] (default -25..25 mm/m by 0.05)
%           .pairs     np x 2 pass indices of each pair: enables the jackknife
%           .dt        np-vector, pair time separation [days]: enables the
%                      trend-controlled scan
%           .vgrid     trial rates [m per day] (default -6..6 mm/day by 0.05)
%           .dres      nb x np residual misalignment [ns] (optional)
%           .dq        np-vector, pair quadrature-tide difference [m]:
%                      enables the quadrature fit
%           .ggrid     trial artefact gains [m/m per ns] (default -20..20 by 0.25)
%           .nboot     bootstrap resamples over pairs (default 150)
%           .min_pairs cells with fewer usable pairs return NaN (default 6)
%           .self_test inject this response [m/m] as a phase rotation and
%                      check it comes back as an exact shift (default
%                      1.2e-3; 0 disables). R.inj_err is the worst error
%                      over interior cells; anything above a grid step
%                      means the phase sign or the arithmetic is wrong.
%           .seed      seed of the bootstrap's own generator (default 3)
%
%   R.a, R.F, R.n_used           nz x nb, the 1-D (tide-only) scan
%   R.a_sd, R.a_sd_boot, R.a_sd_jk  nz x nb standard errors (see above)
%   R.a_t, R.v_t, R.F_t, R.a_t_sd   nz x nb (only with opts.dt; a_t_sd
%                                needs opts.pairs too), and R.r_tt
%   R.a_q, R.b_q, R.v_q, R.F_q   nz x nb (only with opts.dq; v_q needs opts.dt),
%   R.a_q_sd, R.b_q_sd           their jackknife errors (needs opts.pairs), R.r_qt
%   R.a2, R.g2, R.F2             nz x nb (only with opts.dres)
%   R.rcol                       nb x 1, corr(dres, dtide) per block
%   R.inj_err, R.n_edge          self-test error and grid-edge cell count
%
%   See also vdef.complexBlocks, scripts/diagnostics/tidal_stack.m.

if nargin < 5 || isempty(opts), opts = struct(); end
def = struct('agrid', (-25:0.05:25)*1e-3, 'ggrid', (-20:0.25:20)*1e-3, 'dres', [], ...
             'pairs', [], 'dt', [], 'dq', [], 'vgrid', (-6:0.05:6)*1e-3, ...
             'nboot', 150, 'min_pairs', 6, 'self_test', 1.2e-3, 'seed', 3);
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end, end
[nz, nb, np] = size(I);
dtide = dtide(:).'; kz = kz(:);
assert(numel(dtide) == np, 'dtide has %d entries for %d pairs', numel(dtide), np);
assert(numel(kz) == nz, 'kz has %d entries for %d depths', numel(kz), nz);
ag = opts.agrid(:).';
pairs = opts.pairs;
if ~isempty(pairs)
  assert(isequal(size(pairs), [np 2]), 'opts.pairs must be np x 2');
end

[a, F, nu] = scan1(I, W, dtide, kz, ag, opts.min_pairs);
R = struct('a', a, 'F', F, 'n_used', nu, 'agrid', ag);

% self-test: a known response injected as a phase rotation must return as
% an exact shift of the estimate, cell by cell, away from the grid edge
R.inj_err = NaN; R.n_edge = 0;
if opts.self_test > 0
  Ii = I;
  for p = 1:np, Ii(:,:,p) = I(:,:,p) .* exp(1i * kz * opts.self_test * dtide(p)); end
  ai = scan1(Ii, W, dtide, kz, ag, opts.min_pairs);
  sh = ai - a;
  edge = abs(a) > max(ag) - 1.5*opts.self_test;
  ok = isfinite(sh) & ~edge;
  R.inj_err = max(abs(sh(ok) - opts.self_test));
  R.n_edge = nnz(edge);
end

% bootstrap over pairs: the within-pair scatter, kept for comparison with
% the jackknife. It draws from a generator of its own: reseeding the
% global stream (rng) would reset the caller's random state as a side
% effect, and MATLAB refuses rng() once legacy seeding (rand('seed', ...))
% is active, which the Octave-compatible tests use
bidx = boot_indices(opts.seed, np, opts.nboot);
ab = nan(nz, nb, opts.nboot);
for r = 1:opts.nboot
  s = bidx(:, r);
  ab(:,:,r) = scan1(I(:,:,s), W(:,:,s), dtide(s), kz, ag, opts.min_pairs);
end
R.a_sd_boot = std_omitnan(ab, 3);     % std(ab, 0, 3, 'omitnan'), which Octave 8 does not accept
R.a_sd = R.a_sd_boot;

% delete-one-pass jackknife: the error bar that respects shared passes
R.a_sd_jk = [];
if ~isempty(pairs)
  R.a_sd_jk = jackknife_sd(@(keep) scan1(I(:,:,keep), W(:,:,keep), dtide(keep), kz, ag, opts.min_pairs), pairs);
  R.a_sd = R.a_sd_jk;
end

% trend-controlled scan: Phi = kz*(a*dtide + v*dt)
if ~isempty(opts.dt)
  dt = opts.dt(:).';
  assert(numel(dt) == np, 'opts.dt has %d entries for %d pairs', numel(dt), np);
  vg = opts.vgrid(:).';
  reg = repmat(dt, nb, 1);
  [R.a_t, R.v_t, R.F_t, rtt] = scan2(I, W, dtide, kz, ag, reg, vg, opts.min_pairs);
  R.r_tt = rtt(1);                    % the same for every block
  R.vgrid = vg;
  R.a_t_sd = [];
  if ~isempty(pairs)
    R.a_t_sd = jackknife_sd(@(keep) scan2(I(:,:,keep), W(:,:,keep), dtide(keep), kz, ag, ...
      reg(:,keep), vg, opts.min_pairs), pairs);
  end
end

% quadrature-controlled fit: Phi = kz*(a*dtide + b*dq [+ v*dt]), started at
% the in-phase grid peak with b = 0
if ~isempty(opts.dq)
  dq = opts.dq(:).';
  assert(numel(dq) == np, 'opts.dq has %d entries for %d pairs', numel(dq), np);
  cq = corrcoef(dq, dtide); R.r_qt = cq(1,2);
  if ~isempty(opts.dt)
    X = [dtide; dq; opts.dt(:).']; TH0 = cat(3, R.a_t, zeros(nz, nb), R.v_t);
  else
    X = [dtide; dq]; TH0 = cat(3, R.a, zeros(nz, nb));
  end
  [TH, R.F_q] = fitlin(I, W, X, kz, TH0, opts.min_pairs);
  R.a_q = TH(:,:,1); R.b_q = TH(:,:,2);
  if size(TH, 3) > 2, R.v_q = TH(:,:,3); end
  R.a_q_sd = []; R.b_q_sd = [];
  if ~isempty(pairs)
    % each delete-one fit starts from the full-data peak
    sd = jackknife_sd(@(keep) fitlin(I(:,:,keep), W(:,:,keep), X(:,keep), kz, TH, opts.min_pairs), pairs);
    R.a_q_sd = sd(:,:,1); R.b_q_sd = sd(:,:,2);
  end
end

% artefact-controlled scan: Phi = kz*(a*dtide + g*dres)
if ~isempty(opts.dres)
  dres = opts.dres;
  assert(isequal(size(dres), [nb np]), 'opts.dres must be nb x np');
  gg = opts.ggrid(:).';
  [R.a2, R.g2, R.F2, R.rcol] = scan2(I, W, dtide, kz, ag, dres, gg, opts.min_pairs);
  R.ggrid = gg;
end
end

%% ========================================================================
function [a, F, nu] = scan1(I, W, dtide, kz, ag, min_pairs)
[nz, nb, ~] = size(I);
a = nan(nz, nb); F = nan(nz, nb); nu = zeros(nz, nb);
for z = 1:nz
  ph = exp(-1i * kz(z) * (dtide(:) * ag));            % np x na
  for b = 1:nb
    Ip = reshape(I(z,b,:), 1, []); Wp = reshape(W(z,b,:), 1, []);
    g = Wp > 0 & isfinite(Ip);
    nu(z,b) = nnz(g);
    if nu(z,b) < min_pairs, continue; end
    Fv = abs((Ip(g) .* Wp(g)) * ph(g,:)) / sum(abs(Ip(g)) .* Wp(g));
    [F(z,b), im] = max(Fv);
    a(z,b) = ag(im);
  end
end
end

%% ========================================================================
function [a2, g2, F2, rcol] = scan2(I, W, dtide, kz, ag, reg, gg, min_pairs)
%SCAN2 Joint scan over the response a and the gain g of a second regressor
%   reg (nb x np): the trial pair at which the pairs add most coherently.
%   rcol is corr(reg, dtide) per block - near +/-1 the surface is a ridge
%   and a2 is not identified.
[nz, nb, ~] = size(I);
a2 = nan(nz, nb); g2 = nan(nz, nb); F2 = nan(nz, nb); rcol = nan(nb, 1);
for b = 1:nb
  d = reg(b,:); gd = isfinite(d);
  if nnz(gd) >= min_pairs && std(d(gd)) > 0 && std(dtide(gd)) > 0
    cc = corrcoef(d(gd), dtide(gd)); rcol(b) = cc(1,2);
  end
  for z = 1:nz
    Ip = reshape(I(z,b,:), 1, []); Wp = reshape(W(z,b,:), 1, []);
    g = Wp > 0 & isfinite(Ip) & gd;
    if nnz(g) < max(min_pairs, 8), continue; end
    PA = exp(-1i * kz(z) * (dtide(g).' * ag));      % np x na
    PG = exp(-1i * kz(z) * (d(g).' * gg));          % np x ng
    Fm = abs(PG.' * ((Ip(g).*Wp(g)).' .* PA)) / sum(abs(Ip(g)).*Wp(g));
    [F2(z,b), im] = max(Fm(:)); [ig, ia] = ind2sub(size(Fm), im);
    a2(z,b) = ag(ia); g2(z,b) = gg(ig);
  end
end
end

%% ========================================================================
function [TH, F] = fitlin(I, W, X, kz, TH0, min_pairs)
%FITLIN Continuous maximum of the coherent sum over k linear regressors.
%   Phi_p = kz(z) * X(:,p).' * th; maximises F = |sum_p I_p W_p exp(-i Phi_p)|
%   / sum_p |I_p| W_p from TH0 (nz x nb x k) by Gauss-Newton on the pair
%   phases about the current common phase, accepting a step only if F
%   rises (halving it otherwise). Cells with no finite start or too few
%   pairs are NaN.
[nz, nb, ~] = size(I); k = size(X, 1);
TH = nan(nz, nb, k); F = nan(nz, nb);
for b = 1:nb
  for z = 1:nz
    th = reshape(TH0(z,b,:), [], 1);
    Ip = reshape(I(z,b,:), [], 1); Wp = reshape(W(z,b,:), [], 1);
    g = Wp > 0 & isfinite(Ip) & all(isfinite(X), 1).';
    if nnz(g) < max(min_pairs, 8) || any(~isfinite(th)), continue; end
    Zp = Ip(g) .* Wp(g); m = abs(Zp); A = kz(z) * X(:, g).'; Aa = [ones(nnz(g), 1) A];
    f = abs(sum(Zp .* exp(-1i * A * th)));
    for it = 1:50
      S = sum(Zp .* exp(-1i * A * th));
      r = angle(Zp .* exp(-1i * (A * th + angle(S))));
      w = m .* max(cos(r), 0.05);
      d = (Aa.' * (w .* Aa)) \ (Aa.' * (m .* sin(r)));
      st = d(2:end); s = 1; moved = false;
      while s > 1e-4
        f1 = abs(sum(Zp .* exp(-1i * A * (th + s*st))));
        if f1 >= f, th = th + s*st; f = f1; moved = true; break; end
        s = s/2;
      end
      if ~moved || max(abs(s*st)) < 1e-9, break; end
    end
    TH(z,b,:) = th; F(z,b) = f / sum(m);
  end
end
end

%% ========================================================================
function sd = jackknife_sd(fun, pairs)
%JACKKNIFE_SD Delete-one-pass jackknife standard error of fun's estimate.
%   fun(keep) returns the estimate (nz x nb, or nz x nb x k for k
%   parameters at once) from the pairs flagged in the
%   logical row keep; pairs is np x 2 pass indices. Each pass is left out
%   in turn together with every pair it belongs to, so an error that a
%   pass carries into all of its pairs is seen as one draw, not as many.
%   Efron's (n-1)/n scaling; a cell estimated in fewer than three
%   delete-one subsets is NaN.
passes = unique(pairs(:)).'; n = numel(passes);
leave = @(q) ~any(pairs == passes(q), 2).';
e = cell(1, n);
for q = 1:n, e{q} = fun(leave(q)); end
est = cat(4, e{:});                      % the subsets along dim 4
ok = ~isnan(est); nn = sum(ok, 4);
e0 = est; e0(~ok) = 0; m = sum(e0, 4) ./ max(nn, 1);
d = est - m; d(~ok) = 0;
sd = sqrt((nn - 1) ./ max(nn, 1) .* sum(d.^2, 4));
sd(nn < 3) = NaN;
end

%% ========================================================================
function idx = boot_indices(seed, n, nboot)
%BOOT_INDICES Bootstrap resampling indices (n x nboot, each in 1..n) from a
%   Park-Miller minimal-standard generator: identical in MATLAB and Octave
%   (every product stays below 2^53, so the double arithmetic is exact) and
%   independent of the global random stream.
m = 2147483647; a = 16807;
x = mod(max(1, round(abs(seed))), m - 1) + 1;      % state in 1..m-1
u = zeros(n*nboot, 1);
for k = 1:n*nboot
  x = mod(a*x, m); u(k) = x/m;
end
idx = reshape(min(n, floor(u*n) + 1), n, nboot);
end

%% ========================================================================
function s = std_omitnan(x, dim)
%STD_OMITNAN Standard deviation along dim ignoring NaN, normalised by n-1:
%   std(x, 0, dim, 'omitnan') exactly (0 for a single value, NaN for none),
%   written out because Octave 8's std and var take no NaN flag.
ok = ~isnan(x); n = sum(ok, dim);
x0 = x; x0(~ok) = 0;
d = x - sum(x0, dim) ./ max(n, 1); d(~ok) = 0;
s = sqrt(sum(d.^2, dim) ./ max(n - 1, 1));
s(n == 0) = NaN;
end
