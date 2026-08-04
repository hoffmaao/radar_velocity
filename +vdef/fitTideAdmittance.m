function A = fitTideAdmittance(strain, t, tide, opts)
%FITTIDEADMITTANCE Reference-invariant tidal admittance of block strain.
%   A = FITTIDEADMITTANCE(strain, t, tide, opts) fits, independently in
%   each along-track block, the joint model
%
%     strain(b,k) = a(b) + b(b)*t(k) + c(b)*tide(k)
%
%   by ordinary least squares over the pairs k with finite data, and
%   reports the tide admittance c - the strain response per metre of tide
%   with the secular trend modelled out.
%
%   WHY A JOINT FIT. Every strain series here is measured relative to one
%   reference pass, and so are t and tide, so changing the reference
%   subtracts a constant from all three series. A constant is absorbed
%   entirely by the intercept, which makes the trend b and the admittance
%   c invariant to the reference choice (exactly so when the pair set is
%   unchanged). A PLAIN correlation of strain with tide has no term for
%   the trend, so the trend aliases into it through the sample covariance
%   of t with tide - which depends on which pairs the build happens to
%   contain, and gave the same leg different absolute r in builds
%   referencing different epochs. See scripts/figures/leg1_merge_check.m,
%   whose acceptance metric this fit provides.
%
%   INPUTS. strain is Nblk x Npair (blocks down, pairs across); t [days]
%   and tide [m] are Npair-vectors, any origin - both are centred
%   internally, so only differences matter. opts fields, all optional:
%     .min_pairs         finite pairs required per block (default 6; the
%                        hard floor is 4, leaving one degree of freedom)
%     .max_collinearity  skip a block when |corr(t,tide)| over its usable
%                        pairs exceeds this (default 0.99): b and c are
%                        then jointly undetermined and both come back
%                        meaningless rather than visibly broken
%
%   Returns per-block row vectors (1 x Nblk):
%     A.intercept   fitted strain at the mean t and mean tide of the pairs
%                   used - a nuisance parameter that absorbs the reference
%                   constants, reported only for reconstruction
%     A.trend       b [strain/day]
%     A.admittance  c [strain per metre of tide]
%     A.trend_std, A.admittance_std   1-sigma from the fit residuals
%     A.r_partial   partial correlation of strain with tide given the
%                   trend, r = tc/sqrt(tc^2 + dof) with tc the t-statistic
%                   of c. Bounded [-1,1] like the plain correlation it
%                   replaces, so change-point thresholds carry over; +/-1
%                   when the fit is exact.
%     A.p_adm       two-sided Student-t p-value for c
%     A.n           pairs used
%     A.dof         n - 3
%     A.r_tt        corr(t, tide) over the pairs used - the collinearity
%                   that decides how separable trend and tide are
%
%   See also vdef.invertStrainRate, scripts/figures/tidal_deformation.m.

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts,'min_pairs') || isempty(opts.min_pairs)
  opts.min_pairs = 6;
end
if ~isfield(opts,'max_collinearity') || isempty(opts.max_collinearity)
  opts.max_collinearity = 0.99;
end

t    = t(:).';
tide = tide(:).';
Nblk  = size(strain, 1);
Npair = size(strain, 2);
assert(numel(t) == Npair && numel(tide) == Npair, ...
  'strain has %d pairs but t has %d and tide has %d', Npair, numel(t), numel(tide));

A = [];
A.intercept      = nan(1, Nblk);
A.trend          = nan(1, Nblk);
A.admittance     = nan(1, Nblk);
A.trend_std      = nan(1, Nblk);
A.admittance_std = nan(1, Nblk);
A.r_partial      = nan(1, Nblk);
A.p_adm          = nan(1, Nblk);
A.n              = zeros(1, Nblk);
A.dof            = nan(1, Nblk);
A.r_tt           = nan(1, Nblk);

for b = 1:Nblk
  y  = strain(b,:);
  ok = isfinite(y) & isfinite(t) & isfinite(tide);
  n  = nnz(ok);
  A.n(b) = n;
  if n < max(opts.min_pairs, 4)
    continue;
  end

  % Centring the predictors makes the normal equations well scaled and the
  % intercept the value at the sample mean; the slopes are unaffected.
  tc = t(ok)  - mean(t(ok));
  hc = tide(ok) - mean(tide(ok));
  yv = y(ok).';

  st = sqrt(sum(tc.^2));
  sh = sqrt(sum(hc.^2));
  if st == 0 || sh == 0
    % All usable pairs share one time or one tide value: with no spread
    % there is nothing to regress against. r_tt stays NaN.
    continue;
  end
  r_tt = (tc*hc.') / (st*sh);
  A.r_tt(b) = r_tt;
  if abs(r_tt) > opts.max_collinearity
    continue;
  end

  X    = [ones(n,1), tc(:), hc(:)];
  beta = X \ yv;

  resid  = yv - X*beta;
  dof    = n - 3;
  sigma2 = sum(resid.^2) / dof;
  Cov    = sigma2 * ((X'*X) \ eye(3));
  bstd   = sqrt(abs(diag(Cov)));

  A.intercept(b)      = beta(1);
  A.trend(b)          = beta(2);
  A.admittance(b)     = beta(3);
  A.trend_std(b)      = bstd(2);
  A.admittance_std(b) = bstd(3);
  A.dof(b)            = dof;

  if bstd(3) > 0
    tstat = beta(3) / bstd(3);
    A.r_partial(b) = tstat / sqrt(tstat^2 + dof);
    % Two-sided Student-t p-value via the exact incomplete-beta form
    % (betainc is base MATLAB and Octave; tcdf would need a toolbox)
    A.p_adm(b) = betainc(dof/(dof + tstat^2), dof/2, 0.5);
  else
    % An exact fit: the partial correlation saturates at the sign of c
    A.r_partial(b) = sign(beta(3));
    A.p_adm(b) = 0;
  end
end

end
