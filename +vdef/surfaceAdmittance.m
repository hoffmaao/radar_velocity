function S = surfaceAdmittance(Z, tide, opts)
%SURFACEADMITTANCE Tidal admittance of the ice surface from repeat-pass heights.
%   S = SURFACEADMITTANCE(Z, tide, opts) regresses, in each along-track
%   block, the block-mean platform height of every pass on the tide across
%   passes, and returns the slope a(x): the local vertical admittance of
%   the surface to the tide, ~1 where the shelf floats freely and falling
%   toward 0 into the grounding zone. It is the observable that the beam
%   inversion (vdef.invertElasticModulus) fits.
%
%   THE REGRESSOR IS THE POINT OF THIS FUNCTION. The first version of this
%   chain used each pass's own line mean of Z as its tide. That is not the
%   tide: on the EAGER 2022 lines the per-pass line means depart from the
%   CATS2008 prediction by 13-15 cm rms, and about half of that is a
%   height error UNIFORM along the line (a per-pass platform/GPS offset).
%   A uniform per-pass error enters every block AND the line-mean
%   regressor with coefficient one, so the regression slope is pulled
%   toward one - errors in the regressor - and the a(x) profile is
%   compressed toward a flat line. The beam fit cannot tell a compressed
%   profile from a stiffer beam: on synthetic data with the survey's
%   geometry, 14 cm rms of uniform offsets biased E* by +58% with the
%   line-mean regressor while the same fit against the true tide was
%   unbiased. So the tide is supplied EXTERNALLY here, and the line mean
%   is kept only for two nuisance roles:
%     * a PASS GATE: a pass whose line mean departs from the tide by more
%       than opts.max_pass_sigma robust sigmas is rejected. This is what
%       catches a pass with a wrong antenna height or a bad GPS solution -
%       GL2's first pass sat 1.13 m off the tide and, left in, produced a
%       block with negative admittance and a modulus four times the other
%       legs.
%     * a COMMON-MODE NUISANCE COLUMN: the residual of the line mean about
%       its fit to the tide is added as a third regressor, so the part of
%       each block's scatter that is shared by every block on the pass is
%       modelled rather than counted as block noise. It is orthogonal to
%       the tide by construction, so it does not change a(x); it changes
%       a_std, which otherwise carries the shared error and overstates the
%       block's own noise several-fold.
%   What this does NOT fix: a uniform offset still projects onto the
%   tide through the sample covariance over ~13 passes, which shifts the
%   whole profile by a random constant per line. No two-way model can
%   separate that from a uniform shift of a(x) (a uniform offset u_k is
%   exactly confounded with m*tide_k against a(x) - m). Its effect on E*
%   is quantified by a jackknife over passes in the driver, and removed
%   only by an independent height measurement.
%
%   INPUTS
%     Z      Nx x Np platform heights on the main-pass along-track axis
%            [m], one column per pass; a pass with no heights is a NaN
%            column. On the multipass products this is pass(k).ref_z.
%     tide   1 x Np external tide at each pass [m] (any datum), or []
%            to use the line mean of Z as the regressor - the legacy
%            estimator, kept for comparison and refused nothing, but see
%            above for why it is biased.
%     opts   .block           samples per block (default 200)
%            .min_pass        passes needed to regress a block (default 5)
%            .max_pass_sigma  pass-gate threshold in robust sigmas of the
%                             line-mean residual (default 4; Inf disables;
%                             only possible with an external tide)
%            .common_mode     include the line-mean residual as a nuisance
%                             regressor (default true; only meaningful with
%                             an external tide, where it is orthogonal to it)
%
%   RETURNS
%     S.a, S.a_std   nb x 1 admittance per block and its 1-sigma, in
%                    metres of surface per metre of tide (per metre of
%                    line mean when tide is [])
%     S.n            nb x 1 passes used per block
%     S.cols         1 x nb cell of sample indices per block
%     S.pass_ok      1 x Np passes accepted by the gate (and finite)
%     S.pass_resid   1 x Np line-mean residual about its fit to the tide
%     S.abar         slope of the line mean on the tide - the line-mean
%                    admittance, which is what the legacy normalisation
%                    divided by
%     S.tide         the regressor actually used, 1 x Np
%     S.common_mode  whether the nuisance column was used
%
%   See also vdef.invertElasticModulus, scripts/diagnostics/elastic_modulus.m.

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts,'block')          || isempty(opts.block),          opts.block = 200;        end
if ~isfield(opts,'min_pass')       || isempty(opts.min_pass),       opts.min_pass = 5;       end
if ~isfield(opts,'max_pass_sigma') || isempty(opts.max_pass_sigma), opts.max_pass_sigma = 4; end
if ~isfield(opts,'common_mode')    || isempty(opts.common_mode),    opts.common_mode = true; end

[Nx, Np] = size(Z);
linemean = mean(Z, 1, 'omitnan');
have_pass = isfinite(linemean);

external = ~isempty(tide);
if external
  tide = tide(:).';
  assert(numel(tide) == Np, 'tide has %d values for %d passes', numel(tide), Np);
else
  tide = linemean;
end
pass_ok = have_pass & isfinite(tide);
assert(nnz(pass_ok) >= opts.min_pass, ...
  'only %d passes have both heights and a tide; need %d', nnz(pass_ok), opts.min_pass);

%% Pass gate and common-mode residual, from the line mean against the tide
pass_resid = nan(1, Np);
abar = NaN;
if external
  [abar, pass_resid] = linefit(linemean, tide, pass_ok);
  if isfinite(opts.max_pass_sigma)
    r  = pass_resid(pass_ok);
    s  = 1.4826 * median(abs(r - median(r)));
    if s > 0
      bad = pass_ok & abs(pass_resid - median(r)) > opts.max_pass_sigma * s;
      if any(bad)
        pass_ok(bad) = false;
        assert(nnz(pass_ok) >= opts.min_pass, ...
          'the pass gate left %d passes; need %d', nnz(pass_ok), opts.min_pass);
        % Refit without the rejected passes so the nuisance column is not
        % shaped by the outlier it just removed.
        [abar, pass_resid] = linefit(linemean, tide, pass_ok);
      end
    end
  end
end
common_mode = external && opts.common_mode;

%% Per-block regression
nb = floor(Nx / opts.block);
a = nan(nb,1); a_std = nan(nb,1); n = zeros(nb,1); cols = cell(1,nb);
for b = 1:nb
  idx = (b-1)*opts.block+1 : b*opts.block;
  cols{b} = idx;
  zb = mean(Z(idx,:), 1, 'omitnan');
  ok = pass_ok & isfinite(zb);
  m  = nnz(ok);
  if common_mode
    X = [ones(m,1), tide(ok).', pass_resid(ok).'];
  else
    X = [ones(m,1), tide(ok).'];
  end
  p = size(X,2);
  if m < max(opts.min_pass, p+1), continue; end
  if rcond(X.'*X) < 1e-12, continue; end
  beta = X \ zb(ok).';
  res  = zb(ok).' - X*beta;
  G    = (X.'*X) \ eye(p);
  a(b)     = beta(2);
  a_std(b) = sqrt(sum(res.^2)/(m - p) * G(2,2));
  n(b)     = m;
end

S = struct('a', a, 'a_std', a_std, 'n', n, 'cols', {cols}, ...
  'pass_ok', pass_ok, 'pass_resid', pass_resid, 'abar', abar, ...
  'tide', tide, 'common_mode', common_mode);
end

%% ========================================================================
function [slope, resid] = linefit(y, x, ok)
% Fitted over the accepted passes; the residual is reported for every pass
% with data, so a rejected pass still shows how far off it was.
X = [ones(nnz(ok),1), x(ok).'];
b = X \ y(ok).';
slope = b(2);
resid = nan(size(y));
have  = isfinite(y) & isfinite(x);
resid(have) = y(have) - (b(1) + b(2)*x(have));
end
