function N = invertNetwork(pairs, d, opts)
%INVERTNETWORK Per-epoch displacement from a redundant network of pairs.
%   N = INVERTNETWORK(pairs, d, opts) solves the overdetermined system
%
%       d(p) = x(j_p) - x(i_p)
%
%   for the per-epoch values x, given measurements d over pairs (i,j). Each
%   ROW of d is solved independently, so a whole depth profile or a set of
%   along-track blocks can be inverted in one call.
%
%   WHY THIS EXISTS. Pairing every pass against ONE reference uses N-1 of
%   the N(N-1)/2 available pairs - 12 of 78 for 13 passes. That throws away
%   the redundancy, and with it both the averaging and the only internal
%   check the data admits. Worse, it makes the answer depend on which pass
%   was chosen as reference: for an end-referenced product the interval
%   sampling is one-sided, and any error that grows with separation lands
%   entirely in the fitted rate.
%
%   A network inversion fixes all three. It uses every pair, it has no
%   privileged epoch, and the fit residuals ARE the closure errors, so
%   inconsistent pairs identify themselves instead of having to be guessed
%   at from quality metrics.
%
%   THE DATUM. x is determined only up to an additive constant, since the
%   data constrain differences alone. The constraint applied here is
%   sum(x) = 0 rather than x(ref) = 0, so no epoch is privileged. It makes
%   no difference downstream: a joint strain = a + b*t + c*tide fit absorbs
%   any constant into the intercept, which is the whole reason that fit is
%   reference-invariant.
%
%   ROBUSTNESS. Solve, measure residuals, reject pairs beyond
%   opts.n_sigma robust sigmas, re-solve, repeat. Rejection can disconnect
%   the network - an epoch reachable through no surviving pair cannot be
%   recovered and is returned NaN rather than silently pinned to the datum.
%
%   pairs  Npair x 2 of epoch indices; the measurement is x(j) - x(i)
%   d      Nrow x Npair measurements (NaN where a pair has no value)
%   opts   .n_sigma      rejection threshold in robust sigmas (default 3)
%          .max_iter     reweighting passes (default 4)
%          .min_pairs    fewest surviving pairs to attempt a row (default 6)
%          .weights      1 x Npair relative weights (default all ones)
%
%   Returns
%     N.x          Nrow x Nepoch solved values, NaN where unreachable
%     N.x_std      Nrow x Nepoch 1-sigma from the residual scatter
%     N.resid      Nrow x Npair fit residuals - the CLOSURE errors
%     N.used       Nrow x Npair logical, false where rejected
%     N.rms        Nrow x 1 residual rms per row
%     N.n_used     Nrow x 1 pairs surviving
%     N.n_epoch    Nrow x 1 epochs recovered
%
%   See also vdef.fitTideAdmittance, scripts/diagnostics/closure.m.

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts,'n_sigma')  || isempty(opts.n_sigma),  opts.n_sigma  = 3; end
if ~isfield(opts,'max_iter') || isempty(opts.max_iter), opts.max_iter = 4; end
if ~isfield(opts,'min_pairs')|| isempty(opts.min_pairs),opts.min_pairs= 6; end

pairs = round(pairs);
Npair = size(pairs,1);
assert(size(pairs,2) == 2, 'pairs must be Npair x 2');
assert(size(d,2) == Npair, 'd has %d columns for %d pairs', size(d,2), Npair);
Nep  = max(pairs(:));
Nrow = size(d,1);

if isfield(opts,'weights') && ~isempty(opts.weights)
  w0 = opts.weights(:).';
  assert(numel(w0) == Npair, 'weights must be 1 x Npair');
else
  w0 = ones(1,Npair);
end

% Incidence matrix: one row per pair, -1 at i, +1 at j
A0 = zeros(Npair, Nep);
for p = 1:Npair
  A0(p, pairs(p,1)) = -1;
  A0(p, pairs(p,2)) = +1;
end

N = struct();
N.x       = nan(Nrow, Nep);
N.x_std   = nan(Nrow, Nep);
N.resid   = nan(Nrow, Npair);
N.used    = false(Nrow, Npair);
N.rms     = nan(Nrow, 1);
N.n_used  = zeros(Nrow, 1);
N.n_epoch = zeros(Nrow, 1);

for r = 1:Nrow
  y  = d(r,:).';
  ok = isfinite(y) & isfinite(w0(:)) & w0(:) > 0;
  if nnz(ok) < opts.min_pairs, continue; end

  for it = 1:opts.max_iter
    [x, xs, res, reach] = solve_once(A0, y, ok, w0(:), pairs, Nep);
    if isempty(x), break; end
    rr = res(ok);
    s  = 1.4826 * median(abs(rr - median(rr)));
    if ~(s > 0), break; end
    bad = ok & (abs(res) > opts.n_sigma * s);
    if ~any(bad), break; end
    ok(bad) = false;
    if nnz(ok) < opts.min_pairs, break; end
  end

  if nnz(ok) < opts.min_pairs, continue; end
  [x, xs, res, reach] = solve_once(A0, y, ok, w0(:), pairs, Nep);
  if isempty(x), continue; end

  x(~reach)  = NaN;
  xs(~reach) = NaN;
  N.x(r,:)      = x(:).';
  N.x_std(r,:)  = xs(:).';
  N.resid(r,:)  = res(:).';
  N.used(r,:)   = ok(:).';
  N.rms(r)      = sqrt(mean(res(ok).^2));
  N.n_used(r)   = nnz(ok);
  N.n_epoch(r)  = nnz(reach);
end

end

%% ========================================================================
function [x, xs, res, reach] = solve_once(A0, y, ok, w, pairs, Nep)
% Weighted least squares with the sum(x) = 0 datum appended as one more
% equation. Epochs that no surviving pair touches are reported as
% unreachable rather than being pinned to the datum by the constraint.
x = []; xs = []; res = nan(size(y)); reach = false(1,Nep);

reach = connected(pairs(ok,:), Nep);
if ~any(reach), return; end

A = A0(ok,:);
W = sqrt(w(ok));
% datum row, weighted like a single average observation so it fixes the
% constant without competing with the data
Ad = [bsxfun(@times, W, A); ones(1,Nep)];
yd = [W .* y(ok); 0];

% columns for unreachable epochs carry no information; drop them so the
% system is not rank deficient for a reason the caller cannot see
keep = reach;
Ak = Ad(:, keep);
if rank(Ak) < nnz(keep), return; end

xk = Ak \ yd;
x = zeros(Nep,1); x(keep) = xk;

res_ok = A*x - y(ok);
res(ok) = res_ok;

dof = max(1, nnz(ok) - (nnz(keep) - 1));
s2  = (res_ok' * res_ok) / dof;
C   = s2 * pinv(Ak' * Ak);
xs  = zeros(Nep,1);
xs(keep) = sqrt(abs(diag(C)));
end

%% ========================================================================
function reach = connected(pr, Nep)
% Epochs in the largest connected component of the surviving pair graph.
% A network split into pieces has no common datum between them, so only
% the largest piece is solvable against one constraint.
reach = false(1,Nep);
if isempty(pr), return; end
adj = false(Nep);
for p = 1:size(pr,1)
  adj(pr(p,1),pr(p,2)) = true;
  adj(pr(p,2),pr(p,1)) = true;
end
best = false(1,Nep);
seen = false(1,Nep);
for s = 1:Nep
  if seen(s) || ~any(adj(s,:)), continue; end
  comp = false(1,Nep); stack = s;
  while ~isempty(stack)
    v = stack(end); stack(end) = [];
    if comp(v), continue; end
    comp(v) = true;
    nb = find(adj(v,:) & ~comp);
    stack = [stack nb]; %#ok<AGROW>
  end
  seen = seen | comp;
  if nnz(comp) > nnz(best), best = comp; end
end
reach = best;
end
