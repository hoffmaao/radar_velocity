function R = invertElasticModulus(x_obs, w_obs, opts)
%INVERTELASTICMODULUS Effective Young's modulus from an observed flexure profile.
%   R = INVERTELASTICMODULUS(x_obs, w_obs, opts) fits the elastic-beam
%   flexure model of vdef.beamFlexure to an observed profile of tidal
%   deflection across a grounding zone, by least squares over the effective
%   Young's modulus E* and the position of the landward boundary. It follows
%   the method of Elgart, Minchew and Meyer (2025) for ICESat-2 flexure on
%   the Ross Ice Shelf, with the ice thickness supplied rather than fitted.
%
%   WHAT IS AND IS NOT BEING MEASURED. The flexure of a beam depends on its
%   rigidity D = E h^3 / (12(1-nu^2)) and on nothing else, so a flexure
%   profile constrains only the PRODUCT. E* here is what D implies once the
%   thickness is taken as known, and it inherits the thickness error
%   amplified threefold: dE/E = -3 dh/h. R.dlogE_dlogh states that
%   explicitly for the fitted solution rather than leaving it as a caveat.
%   Nor does the fit separate flexure of the ice from flexure of the bed
%   beneath it, or from the anelastic part of the response at tidal
%   frequency, which is why the result is an EFFECTIVE modulus and typically
%   lands well below the ~9 GPa of laboratory ice.
%
%   THREE PARAMETERS, TWO SEARCHED. For a candidate (E*, x0):
%     * x0 is the landward boundary, where the beam is clamped. It is
%       searched rather than fixed because the inward limit of tidal flexure
%       migrates on tidal timescales and is not the mapped grounding line -
%       ice is thick, so some flexure is carried upstream of it. Predicted
%       deflection is identically zero landward of x0, so observations there
%       stay in the fit and constrain x0 instead of being discarded.
%     * the amplitude - the far-field tidal deflection A0 - enters the model
%       linearly, so it is not searched. At each (E*, x0) the weighted least
%       squares optimum is written down in closed form and substituted back.
%       That is what makes this usable on a tidal ADMITTANCE profile a(x),
%       whose overall scale is set by an arbitrary normalisation (a line mean
%       for GPS-derived a(x)) rather than by the tide: the scale cancels, and
%       only the SHAPE of the profile is being fitted.
%
%   ILL-POSEDNESS. The problem is ill-posed in the sense of Lucchinetti and
%   Stussi (2002), and E* trades off against x0 - a stiffer beam whose clamp
%   sits further landward can reproduce much the same curve over a short
%   window. The full misfit surface is returned for exactly this reason, and
%   R.interior and R.n_local_min say whether the minimum found is a single
%   interior one. A result with R.interior false has run to the edge of a
%   search grid and means nothing on its own. Elgart and others handle the
%   same problem by requiring the minimum to be well behaved and consistent
%   between neighbouring ground tracks; the analogue here is consistency
%   between the survey legs.
%
%   INPUTS
%     x_obs  M-vector of along-track position [m], increasing SEAWARD. The
%            orientation is checked: w_obs must trend upward with x, since
%            deflection grows from zero at the clamp to the far-field tide.
%     w_obs  M-vector of observed deflection, in any units - metres of
%            deflection, or metres of deflection per metre of tide (an
%            admittance). NaNs are dropped.
%     opts   .h            thickness [m]: a scalar, an Nh x 2 table
%                          [x_h h] interpolated onto the model grid
%                          (linear inside, held constant outside - a linear
%                          extrapolation of thickness over kilometres is
%                          not defensible), or a handle @(x) h
%            .sigma        M-vector 1-sigma on w_obs, for weighting and for
%                          the reduced chi-squared (default: unweighted)
%            .E_grid       trial moduli [Pa] (default 61 points log-spaced
%                          over 0.05-30 GPa, spanning the published range)
%            .x0_grid      trial landward boundaries [m], absolute (default
%                          x0_init + (-4000:250:4000))
%            .x0_init      centre of the default x0 grid (default x_obs(1))
%            .seaward_pad  how far past the last observation the seaward
%                          boundary is placed, in flexural lengths
%                          (default 4). w = A0 and dw/dx = 0 there is a
%                          far-field statement, so too small a pad stiffens
%                          the beam and biases E* upward.
%            .avg_width    width [m] over which each observation is an
%                          AVERAGE, scalar or per-observation (default 0,
%                          point samples). Set it whenever the profile came
%                          from along-track blocks. The deflection is curved
%                          on the scale of a flexural length, so a block mean
%                          is not the value at the block centre, and the
%                          difference is a systematic that the fit cannot
%                          absorb: on 500 m blocks over a ~1 km flexural
%                          length it reaches ten to thirty times the formal
%                          error on a(x), leaves a misfit valley with several
%                          shallow bottoms, and makes the recovered E* jump
%                          around under a 1% change in assumed thickness.
%            .dx_max       coarsest model grid spacing [m] (default 50; the
%                          grid is the finer of this and lambda/30)
%            .max_nodes    node ceiling per solve (default 1500)
%            .fit_offset   also fit an additive constant (default false).
%                          A robustness check: a GPS baseline ramp or a
%                          mis-set reference shows up as a nonzero offset.
%                          Shape rows only - the strain admittance is
%                          reference-invariant and gets no offset.
%            .strain       OPTIONAL SECOND DATASET: the englacial strain
%                          admittance, fitted JOINTLY with the deflection.
%                          Fields:
%                            .x          positions [m], same seaward frame
%                            .y          dh at ref_depth per metre of tide
%                                        [m/m], surface-referenced, dh < 0
%                                        = column shortened
%                            .sigma      1-sigma on y
%                            .ref_depth  depth the column change is read
%                                        at [m] (default 100)
%                            .avg_width  block width [m] (default
%                                        opts.avg_width)
%                          WHY IT HELPS, and why it is not just more
%                          points: the shape fit eliminates the amplitude,
%                          so a line-mean-normalised a(x) constrains only
%                          the SHAPE of the beam. The strain admittance is
%                          absolute, and the beam predicts it with the
%                          SAME shared amplitude - dh = amp*K2(x)*w''(x) -
%                          whose curvature scale goes as D^(-1/2). It is
%                          therefore an amplitude equation the shape data
%                          cannot supply, and it is what closes the E*-x0
%                          trade-off valley on a window that never sees
%                          the far field. The prediction is thin-plate
%                          bending only (eps_zz = nu/(1-nu)*(z_n - z)*w'',
%                          z_n = h(x)/2); any non-flexural tidal strain in
%                          the data lands in chi2_strain, which is the
%                          number to check before believing the joint E*.
%            .E_patch      LOCAL mode: struct with fields lo, hi, E_ref
%                          [m, m, Pa]. The searched modulus then applies
%                          ONLY between lo and hi; everywhere else the
%                          beam keeps E_ref. Implemented through the exact
%                          equivalence D = E h^3/(12(1-nu^2)): inside the
%                          patch the solver's thickness is scaled by
%                          (E/E_ref)^(1/3) and the beam solved with E_ref,
%                          so no second solver is needed and the bending
%                          lever arm keeps the TRUE thickness. Fix the
%                          clamp (a one-point x0_grid at the global fit's
%                          x0) when scanning patches, or the clamp will
%                          chase each patch. R.E is then the PATCH
%                          modulus; R.D and R.h_mean still describe the
%                          whole beam and should be read from the global
%                          fit instead. Windows narrower than a flexural
%                          length report smoothed averages - flexure is
%                          nonlocal and that is a resolution limit, not an
%                          implementation one.
%            .beam         extra opts passed through to vdef.beamFlexure
%                          (.nu, .rho_w, .g)
%
%   RETURNS
%     R.E, R.E_lo, R.E_hi   best-fit effective Young's modulus and its
%                           profile-likelihood interval [Pa]. The interval
%                           is where the misfit profiled over x0 rises by
%                           one estimated data variance above its minimum -
%                           the delta-chi-squared = 1 rule with the variance
%                           taken from the fit residuals, so it reflects
%                           scatter about the model and not the input sigma.
%                           NaN when the bound runs off the E grid.
%     R.x0, R.amplitude, R.offset   the other fitted parameters
%     R.D, R.lambda         rigidity and flexural length at the mean fitted
%                           thickness [Pa m^3], [m]
%     R.dlogE_dlogh         -3, restated as the sensitivity that dominates
%     R.w_model             model deflection at x_obs, in data units
%     R.x_grid, R.w_grid    the full fitted solution, for plotting
%     R.rms, R.chi2red      residual rms in data units; reduced chi-squared
%                           (NaN when no sigma was supplied)
%     R.s2                  data variance implied by the fit residuals, and
%                           the unit the interval is measured in: the bounds
%                           are where the profiled misfit rises by one of it
%     R.J, R.E_grid, R.x0_grid   the misfit surface, n_x0 by n_E, on the
%                           search grid before any refinement
%     R.J_profile           misfit profiled over x0, one value per E, with
%                           x0 refined continuously at each E rather than
%                           taken off the grid
%     R.x0_profile          the clamp position each of those used - the
%                           trade-off curve itself, and worth plotting
%     R.interior            minimum is off every grid edge
%     R.n_local_min         distinct RIVAL minima in R.J_profile: separated
%                           from the best by a barrier of more than one data
%                           variance, and within 9 of them of it. More than
%                           one is the ill-posedness showing, and the caller
%                           has to choose between them on external grounds.
%     R.curvature           d2J/d(log10 E)^2 at the minimum, NaN when the
%                           bounds came from walking the profile instead
%     R.n_lambda            seaward pad achieved, in flexural lengths
%     R.n_obs, R.dof
%     R.has_strain, R.n_strain   whether a strain dataset was fitted
%     R.chi2_shape, R.chi2_strain   reduced chi-squared PER DATASET, each
%                           against its own input sigmas - the joint fit
%                           cannot be trusted unless both are order 1
%     R.strain_x/_obs/_sigma/_model   the strain data and the fit at them
%     R.strain_grid         predicted strain admittance on R.x_grid's
%                           beam section, for plotting
%
%   See also vdef.beamFlexure, scripts/diagnostics/elastic_modulus.m.

if nargin < 3 || isempty(opts), opts = struct(); end

x_obs = x_obs(:);
w_obs = w_obs(:);
assert(numel(x_obs) == numel(w_obs), ...
  'x_obs has %d points and w_obs %d', numel(x_obs), numel(w_obs));

if isfield(opts,'sigma') && ~isempty(opts.sigma)
  sigma = opts.sigma(:);
  assert(numel(sigma) == numel(x_obs), 'sigma has %d values for %d points', ...
    numel(sigma), numel(x_obs));
  have_sigma = true;
else
  sigma = ones(size(x_obs));
  have_sigma = false;
end

ok = isfinite(x_obs) & isfinite(w_obs) & isfinite(sigma) & sigma > 0;
x_obs = x_obs(ok); w_obs = w_obs(ok); sigma = sigma(ok);
M = numel(x_obs);
assert(M >= 6, 'need at least 6 usable observations, have %d', M);
assert(all(diff(x_obs) > 0), 'x_obs must be strictly increasing (seaward)');

% Orientation. Getting this backwards is the easy mistake - the survey
% coordinate often runs the other way - and it fails silently, returning
% whatever modulus best fits a mirrored curve. Catch it here.
sx = x_obs - mean(x_obs); sw = w_obs - mean(w_obs);
trend = (sx.'*sw) / max(sqrt((sx.'*sx)*(sw.'*sw)), eps);
assert(trend > 0, ['w_obs falls as x_obs increases, so x is running ' ...
  'landward. x_obs must increase SEAWARD, away from the clamped end ' ...
  '(corr = %.2f).'], trend);

if ~isfield(opts,'h') || isempty(opts.h)
  error('vdef:invertElasticModulus:noThickness', ...
    'opts.h is required: flexure constrains E*h^3, so E* is meaningless without h.');
end
if ~isfield(opts,'x0_init')     || isempty(opts.x0_init),     opts.x0_init = x_obs(1); end
if ~isfield(opts,'E_grid')      || isempty(opts.E_grid)
  opts.E_grid = logspace(log10(0.05e9), log10(30e9), 61);
end
if ~isfield(opts,'x0_grid')     || isempty(opts.x0_grid)
  opts.x0_grid = opts.x0_init + (-4000:250:4000);
end
if ~isfield(opts,'seaward_pad') || isempty(opts.seaward_pad), opts.seaward_pad = 4;  end
if ~isfield(opts,'avg_width')   || isempty(opts.avg_width),   opts.avg_width = 0;    end
if ~isscalar(opts.avg_width)
  opts.avg_width = opts.avg_width(:);
  assert(numel(opts.avg_width) == numel(ok), ...
    'avg_width must be scalar or one value per input observation');
  opts.avg_width = opts.avg_width(ok);
end
if ~isfield(opts,'dx_max')      || isempty(opts.dx_max),      opts.dx_max = 50;      end
if ~isfield(opts,'max_nodes')   || isempty(opts.max_nodes),   opts.max_nodes = 1500; end
if ~isfield(opts,'fit_offset')  || isempty(opts.fit_offset),  opts.fit_offset = false; end
if ~isfield(opts,'beam')        || isempty(opts.beam),        opts.beam = struct();  end

% Optional second dataset: englacial strain admittance. Validated here and
% then carried inside opts so the model evaluators see it everywhere.
has_strain = isfield(opts,'strain') && ~isempty(opts.strain);
sx_ = zeros(0,1); sy_ = zeros(0,1); ss_ = zeros(0,1);
if has_strain
  % The joint fit is refused without shape sigmas, and the reason is worth
  % stating: the two observables share one amplitude, and the strain data
  % constrain the modulus THROUGH that amplitude - dh scales as D^(-1/2)
  % only once amp is pinned by the shape normalisation. With unweighted
  % shape rows against sigma-weighted strain rows, the normal equations
  % let the strain set amp for itself, the amplitude information cancels
  % exactly, and the joint fit comes back WIDER than shape alone while
  % looking perfectly healthy. Found by unit test, kept as a guard.
  assert(have_sigma, ['vdef.invertElasticModulus: a joint fit needs ' ...
    'opts.sigma on the shape data, or the shared amplitude decouples ' ...
    'and the strain constraint silently cancels']);
  st = opts.strain;
  assert(all(isfield(st, {'x','y','sigma'})), ...
    'opts.strain needs fields x, y, sigma');
  sx_ = st.x(:); sy_ = st.y(:); ss_ = st.sigma(:);
  assert(numel(sy_) == numel(sx_) && numel(ss_) == numel(sx_), ...
    'opts.strain fields disagree on length');
  oks = isfinite(sx_) & isfinite(sy_) & isfinite(ss_) & ss_ > 0;
  sx_ = sx_(oks); sy_ = sy_(oks); ss_ = ss_(oks);
  assert(numel(sx_) >= 3, ...
    'opts.strain has %d usable points; need at least 3', numel(sx_));
  if ~isfield(st,'ref_depth') || isempty(st.ref_depth), st.ref_depth = 100; end
  if ~isfield(st,'avg_width') || isempty(st.avg_width)
    st.avg_width = opts.avg_width;
    assert(isscalar(st.avg_width), ...
      'give opts.strain.avg_width explicitly when opts.avg_width is per-observation');
  end
  st.x = sx_;
  opts.strain = st;
  hmin = min(thickness_on(opts.h, sx_));
  assert(st.ref_depth < hmin, ...
    'strain ref_depth %.0f m is below the thinnest ice (%.0f m)', st.ref_depth, hmin);
end
Ns = numel(sx_);

E_grid  = sort(opts.E_grid(:)).';
x0_grid = sort(opts.x0_grid(:)).';
nE  = numel(E_grid);
nx0 = numel(x0_grid);
assert(nE >= 3 && nx0 >= 1, 'need at least 3 trial moduli and 1 trial boundary');
assert(all(x0_grid < x_obs(end)), ...
  'every trial landward boundary must sit landward of the last observation');

% Both datasets stacked once; every misfit evaluation uses the stack. The
% shape rows come first, and Na is what tells fit_linear which rows the
% optional offset column applies to.
u  = 1 ./ sigma.^2;
Na = M;
yy = [w_obs; sy_];
uu = [u; 1 ./ ss_.^2];

%% Misfit surface
J = inf(nx0, nE);
for a = 1:nx0
  for e = 1:nE
    Wm = model_shape(x0_grid(a), E_grid(e), x_obs, opts);
    if isempty(Wm), continue; end
    J(a,e) = fit_linear(Wm, yy, uu, opts.fit_offset, Na);
  end
end
assert(any(isfinite(J(:))), 'no trial model could be solved - check opts.h and the grids');

dlogE = mean(diff(log10(E_grid)));
dx0   = 0; if nx0 > 1, dx0 = mean(diff(x0_grid)); end

%% Profile over the clamp position at every modulus
% Taking min(J) down the x0 grid is NOT a profile likelihood, because a
% discrete x0 cannot follow E along the trade-off valley. The misfit away
% from the optimum then comes out too high, the profile too steep, and the
% interval read off its curvature several times too narrow. Refining x0
% continuously at each E is what makes R.J_profile mean what it says.
Jprof  = inf(1, nE);
X0prof = nan(1, nE);
for e = 1:nE
  [j0, a0] = min(J(:,e));
  if ~isfinite(j0), continue; end
  [Jprof(e), X0prof(e)] = refine_x0(E_grid(e), x0_grid(a0), dx0, j0, ...
                                    x_obs, yy, uu, Na, opts);
end
[Jgrid, ep] = min(Jprof);

%% Pattern search from the grid minimum
% Elgart and others use MATLAB's patternsearch at this step. The same idea
% is written out here so that it needs no toolbox and behaves identically in
% Octave: poll the four axis neighbours at the current step, move to the
% best improvement, halve the step when none of them improves.
%
% ONE ROUND OF LOCAL SAMPLING IS NOT ENOUGH, which is what this replaces.
% E* and x0 trade off along a narrow diagonal valley, so quantising x0 to
% the search grid pushes the coarse minimum's E* off by more than one E
% step - and a refinement that only spans one step then converges to a
% solution that fits ~20 times worse than the truth while looking perfectly
% well behaved. Walking the valley is the whole job.
Jmin   = Jgrid;
E_best = E_grid(ep); x0_best = X0prof(ep);
sE = dlogE; sx = dx0;
for it = 1:400
  poll = [log10(E_best)+sE, x0_best; log10(E_best)-sE, x0_best];
  if sx > 0
    poll = [poll; log10(E_best), x0_best+sx; log10(E_best), x0_best-sx]; %#ok<AGROW>
  end
  Jbest = Jmin; qbest = 0;
  for q = 1:size(poll,1)
    if poll(q,2) >= x_obs(end), continue; end
    Wm = model_shape(poll(q,2), 10^poll(q,1), x_obs, opts);
    if isempty(Wm), continue; end
    Jv = fit_linear(Wm, yy, uu, opts.fit_offset, Na);
    if Jv < Jbest, Jbest = Jv; qbest = q; end
  end
  if qbest > 0
    Jmin = Jbest; E_best = 10^poll(qbest,1); x0_best = poll(qbest,2);
  else
    sE = sE/2; sx = sx/2;
    if sE < dlogE/512, break; end
  end
end

[Wbest, info] = model_shape(x0_best, E_best, x_obs, opts);
[~, amp, off] = fit_linear(Wbest, yy, uu, opts.fit_offset, Na);
resid   = w_obs - (amp*Wbest.W + off);
resid_s = sy_  - amp*Wbest.S;

npar = 2 + 1 + double(opts.fit_offset);          % E*, x0, amplitude, offset
dof  = max(M + Ns - npar, 1);
s2   = Jmin / dof;                               % variance implied by the fit

%% Interval from the profiled misfit
% The misfit is quadratic in log10 E near its minimum, and the interval is
% often NARROWER THAN ONE GRID STEP - a well-sampled profile pins E* to a
% few percent while a search grid spanning three decades cannot step finer
% than ten. Interpolating linearly between grid points in that regime puts
% the bounds in the wrong place, and can even put them on one side of the
% refined minimum. So the bounds come from a parabola through the three
% points bracketing the profile minimum, with its floor pinned to the
% refined minimum. Once the implied half-width exceeds a grid step the
% quadratic approximation has stopped being safe and the profile is by then
% well enough sampled to walk directly, interpolating linearly in log10 E.
uE   = log10(E_grid);
step = dlogE;
lo = NaN; hi = NaN; curv = NaN;
if isfinite(s2) && s2 > 0
  if ep > 1 && ep < nE && all(isfinite(Jprof(ep-1:ep+1)))
    p = polyfit(uE(ep-1:ep+1) - uE(ep), Jprof(ep-1:ep+1), 2);
    curv = p(1);
  end
  if isfinite(curv) && curv > 0 && sqrt(s2/curv) <= step
    du = sqrt(s2/curv);
    lo = 10^(log10(E_best) - du);
    hi = 10^(log10(E_best) + du);
  else
    thr = Jprof(ep) + s2;
    lo  = crossing(E_grid, Jprof, thr, ep, -1);
    hi  = crossing(E_grid, Jprof, thr, ep, +1);
  end
  % An interval that excludes its own point estimate is never the honest
  % answer, whichever branch produced it.
  if isfinite(lo), lo = min(lo, E_best); end
  if isfinite(hi), hi = max(hi, E_best); end
end

% Distinct minima in the profiled misfit. Only RIVALS are counted, on two
% conditions, because a bare local-minimum count says "ambiguous" about
% every result:
%   separated  the barrier back to the best minimum exceeds one data
%              variance. Profiling over a DISCRETE x0 grid takes the lower
%              envelope of a family of curves and leaves a shallow dimple
%              wherever the winning x0 changes hands; this drops those.
%   competitive  it sits within 9 data variances of the best, the
%              delta-chi-squared = 3 sigma level. A clean profile has local
%              minima hundreds of variances up the walls of the valley, and
%              nothing about them is a rival explanation of the data.
% Endpoints are excluded - a run to the edge is what R.interior reports.
nloc = 1;
for e = 2:nE-1
  if e == ep, continue; end
  if ~(isfinite(Jprof(e)) && Jprof(e) < Jprof(e-1) && Jprof(e) < Jprof(e+1))
    continue;
  end
  if e < ep, barrier = max(Jprof(e:ep)) - Jprof(e);
  else,      barrier = max(Jprof(ep:e)) - Jprof(e);
  end
  if barrier > s2 && (Jprof(e) - Jprof(ep)) < 9*s2
    nloc = nloc + 1;
  end
end

%% Full solution on its own grid, for plotting
[xg, hg, ~, E_solve] = model_grid(x0_best, E_best, x_obs, opts);
bopts = opts.beam; bopts.A0 = 1;
Fb = vdef.beamFlexure(xg, hg, E_solve, bopts);
% Grounded ice landward of the clamp: flat, at the fitted offset.
if x_obs(1) < x0_best
  x_pre = linspace(x_obs(1), x0_best, 20).'; x_pre(end) = [];
else
  x_pre = zeros(0,1);
end

R = [];
R.E             = E_best;
R.E_lo          = lo;
R.E_hi          = hi;
R.x0            = x0_best;
R.amplitude     = amp;
R.offset        = off;
R.h_mean        = mean(hg);
R.h_spec        = opts.h;      % echoed so a caller can refit identically
R.D             = E_best * R.h_mean^3 / (12*(1 - beam_nu(opts)^2));
R.lambda        = info.lambda;
R.dlogE_dlogh   = -3;
R.w_model       = amp*Wbest.W + off;
R.x_grid        = [x_pre; xg];
R.w_grid        = [repmat(off, numel(x_pre), 1); amp*Fb.w + off];
R.rms           = sqrt(mean(resid.^2));
R.s2            = s2;
R.chi2red       = NaN;
if have_sigma, R.chi2red = Jmin / dof; end

% The strain dataset, echoed with its model and its own goodness of fit.
% chi2_strain uses the INPUT sigmas, not s2, so a strain misfit cannot
% hide inside a shape-dominated variance estimate - the two datasets have
% very different point counts and precisions, and this is the number that
% says whether one E* really explains both.
R.has_strain = has_strain;
R.n_strain   = Ns;
R.chi2_shape  = NaN;
R.chi2_strain = NaN;
if has_strain
  R.strain_x      = sx_;
  R.strain_obs    = sy_;
  R.strain_sigma  = ss_;
  R.strain_model  = amp*Wbest.S;
  R.chi2_strain   = sum((resid_s ./ ss_).^2) / max(Ns - 1, 1);
  if have_sigma
    R.chi2_shape  = sum((resid ./ sigma).^2) / max(M - npar, 1);
  end
  % The strain admittance on the plotting grid, for drawing the predicted
  % profile rather than only its block averages.
  K2g = strain_lever(opts, Fb.x);
  R.strain_grid = amp * K2g .* Fb.d2wdx2;
end
R.J             = J;
R.J_profile     = Jprof;
R.x0_profile    = X0prof;
R.E_grid        = E_grid;
R.x0_grid       = x0_grid;
R.interior      = (ep > 1) && (ep < nE) && (nx0 == 1 || ...
  (X0prof(ep) > x0_grid(1) && X0prof(ep) < x0_grid(end)));
R.n_local_min   = nloc;
R.curvature     = curv;
R.n_lambda      = info.n_lambda;
R.n_obs         = M;
R.dof           = dof;
R.resid         = resid;

end

%% ========================================================================
function nu = beam_nu(opts)
nu = 0.3;
if isfield(opts,'beam') && isfield(opts.beam,'nu') && ~isempty(opts.beam.nu)
  nu = opts.beam.nu;
end
end

%% ========================================================================
function [xg, hg, lambda, E_solve] = model_grid(x0, E, x_obs, opts)
%MODEL_GRID Uniform grid from the clamp to a far-field boundary.
%   Both the spacing and the seaward extent scale with the flexural length,
%   which is itself a function of E and h, so the grid is rebuilt for every
%   candidate rather than fixed once: a stiff beam needs a longer domain to
%   reach the far field, and a soft one needs finer steps to resolve the
%   hinge.
%
%   In E_patch mode the searched E applies only inside the patch, through
%   the D-equivalent thickness h*(E/E_ref)^(1/3), and the solver runs at
%   E_ref (returned as E_solve). Grid sizing uses E_ref too: the patch is
%   a local perturbation, not what sets the flexural length.
C  = vdef.constants();
nu = beam_nu(opts);
rho_w = C.rho_sea; g = C.g;
if isfield(opts,'beam')
  if isfield(opts.beam,'rho_w') && ~isempty(opts.beam.rho_w), rho_w = opts.beam.rho_w; end
  if isfield(opts.beam,'g')     && ~isempty(opts.beam.g),     g     = opts.beam.g;     end
end

E_solve = E;
patch   = isfield(opts,'E_patch') && ~isempty(opts.E_patch);
if patch, E_solve = opts.E_patch.E_ref; end

h_ref  = mean(thickness_on(opts.h, x_obs));
D_ref  = E_solve * h_ref^3 / (12*(1 - nu^2));
lambda = (4*D_ref / (rho_w*g))^(1/4);

L  = x_obs(end) + opts.seaward_pad*lambda;
dx = min(opts.dx_max, lambda/30);
dx = max(dx, (L - x0)/opts.max_nodes);
N  = max(round((L - x0)/dx) + 1, 7);
xg = x0 + dx*(0:N-1).';
hg = thickness_on(opts.h, xg);
if patch
  in = xg >= opts.E_patch.lo & xg <= opts.E_patch.hi;
  hg(in) = hg(in) * (E / opts.E_patch.E_ref)^(1/3);
end
end

%% ========================================================================
function [Mdl, info] = model_shape(x0, E, x_obs, opts)
%MODEL_SHAPE Unit-amplitude model rows for both observables.
%   Mdl.W is the predicted deflection at the shape observation points;
%   Mdl.S, when opts.strain is present, is the predicted strain admittance
%   at the strain points UP TO the shared amplitude:
%
%       dh(zr) = amp * K2(x) * w''(x),   K2 = nu/(1-nu)*(z_n*zr - zr^2/2)
%
%   with z_n = h(x)/2 the neutral plane and zr the reference depth. This
%   is the thin-plate bending strain eps_zz = nu/(1-nu)*(z_n - z)*w''
%   integrated from the surface to zr - the same expression the project
%   validated against the stacked radar profile via the GPS a''(x), with
%   the ill-conditioned polynomial curvature replaced by the beam's.
%
%   Everything is zero landward of the clamp: grounded ice neither
%   deflects nor bends, so observations there are predictions of zero
%   rather than missing data, and they are what pins x0 down.
%
%   Returns empty for a candidate whose domain is too short to hold a far
%   field - fewer than two flexural lengths between the clamp and the
%   seaward boundary. Such a beam is not badly resolved, it is a different
%   problem, and letting it compete in the misfit would let a stiff
%   solution win by being clipped rather than by fitting.
Mdl = []; info = struct('n_lambda', NaN, 'lambda', NaN);
[xg, hg, lambda, E_solve] = model_grid(x0, E, x_obs, opts);
if (xg(end) - xg(1)) < 2*lambda || numel(xg) < 7
  return;
end
bopts = opts.beam; bopts.A0 = 1;
F = vdef.beamFlexure(xg, hg, E_solve, bopts);

Mdl = struct();
Mdl.W = sample_field(F.x, F.w, x_obs, x0, opts.avg_width);
if isfield(opts,'strain') && ~isempty(opts.strain)
  st  = opts.strain;
  K2  = strain_lever(opts, st.x);
  Mdl.S = K2 .* sample_field(F.x, F.d2wdx2, st.x, x0, st.avg_width);
else
  Mdl.S = zeros(0,1);
end
info.n_lambda = F.n_lambda;
info.lambda   = lambda;
end

%% ========================================================================
function K2 = strain_lever(opts, xq)
%STRAIN_LEVER The depth-integrated bending lever arm at each position.
%   Integral of nu/(1-nu)*(z_n - z) from the surface to the reference
%   depth, with the neutral plane z_n = h(x)/2 following the LOCAL
%   thickness. Positive while zr < z_n everywhere, so on a rising tide the
%   column thickens where the beam is concave-up (at the clamp) and thins
%   where it is concave-down (mid-line) - the sign pattern the radar
%   admittance map shows.
nu = beam_nu(opts);
zr = opts.strain.ref_depth;
zn = thickness_on(opts.h, xq) / 2;
K2 = (nu/(1-nu)) * (zn*zr - zr^2/2);
end

%% ========================================================================
function V = sample_field(Fx, Fv, x_obs, x0, wid)
%SAMPLE_FIELD A model field as the observations actually sample it.
%   Point samples when avg_width is zero; otherwise the mean of the field
%   across each observation's own footprint, taken on 21 points so a 500 m
%   block is sampled every 25 m against a flexural length of a kilometre
%   or more. The part of a footprint landward of the clamp contributes
%   zero, which matters for the one block that straddles it.
M = numel(x_obs);
if isscalar(wid) && wid <= 0
  V  = zeros(M, 1);
  in = x_obs >= x0;
  V(in) = interp1(Fx, Fv, x_obs(in), 'linear');
  V(~isfinite(V)) = 0;
  return;
end
if isscalar(wid), wid = repmat(wid, M, 1); end
ns = 21;                                   % odd, as Simpson's rule needs
XS = repmat(x_obs(:), 1, ns) + wid(:) * linspace(-0.5, 0.5, ns);
WS = zeros(size(XS));
in = XS >= x0;
WS(in) = interp1(Fx, Fv, XS(in), 'linear');
WS(~isfinite(WS)) = 0;
% Simpson rather than a plain mean of the samples. A plain mean carries an
% O(h^2) error that scales with the curvature of the deflection - the very
% quantity being fitted - so at 25 m sampling it leaves a residual several
% times the formal error on a(x) and simply moves the block-centre bias
% somewhere less obvious. Simpson is O(h^4) for the same 21 solves.
sw = ones(1, ns); sw(2:2:ns-1) = 4; sw(3:2:ns-2) = 2;
V  = (WS * sw.') / sum(sw);
end

%% ========================================================================
function [Jb, xb] = refine_x0(E, x0, step, J0, x_obs, yy, uu, Na, opts)
%REFINE_X0 Best clamp position at a fixed modulus, by 1-D pattern search.
%   Started from the best point on the x0 grid and stepped down by halves
%   to a 64th of the grid spacing, which is metres for a grid stepped in
%   hundreds of them.
Jb = J0; xb = x0;
if step <= 0, return; end
s = step;
for it = 1:40
  moved = false;
  for d = [-1 1]
    xc = xb + d*s;
    if xc >= x_obs(end), continue; end
    Wm = model_shape(xc, E, x_obs, opts);
    if isempty(Wm), continue; end
    Jv = fit_linear(Wm, yy, uu, opts.fit_offset, Na);
    if Jv < Jb, Jb = Jv; xb = xc; moved = true; end
  end
  if ~moved
    s = s/2;
    if s < step/64, break; end
  end
end
end

%% ========================================================================
function [J, amp, off] = fit_linear(Mdl, yy, uu, fit_offset, Na)
%FIT_LINEAR Weighted least squares for the parameters the model is linear in.
%   ONE amplitude is shared by both observables: the deflection and the
%   bending strain of the same beam scale together with the far-field
%   tide, so a = amp*W and dh = amp*S with the same amp. That coupling is
%   the entire point of the joint fit - the shape data alone cannot see
%   the amplitude (it is normalised away), while the strain data are
%   absolute, so they contribute the D^(-1/2) curvature-amplitude
%   constraint that shape alone lacks. The optional offset is a shape-only
%   nuisance column: the strain admittance is reference-invariant by
%   construction and gets no such freedom.
amp = NaN; off = 0; J = inf;
rows = [Mdl.W; Mdl.S];
if fit_offset
  X = [rows, [ones(Na,1); zeros(numel(Mdl.S),1)]];
else
  X = rows;
end
Xw = X .* repmat(uu, 1, size(X,2));
A  = X.' * Xw;
if rcond(A) < 1e-12, return; end
p   = A \ (Xw.' * yy);
amp = p(1);
if fit_offset, off = p(2); end
r = yy - X*p;
J = sum(uu .* r.^2);
end

%% ========================================================================
function Ec = crossing(E, Jp, thr, k0, dirn)
%CROSSING Where the profiled misfit first rises through thr, walking out
%   from the minimum. Linear in log E between grid points; NaN if the walk
%   reaches the end of the grid still below the threshold, which says the
%   data do not bound the modulus on that side. Grid entries whose profile
%   could not be evaluated are stepped over without becoming the
%   interpolation anchor, so the bound always sits between two FINITE
%   points of the profile.
Ec = NaN;
n  = numel(E);
k  = k0;
kn = k0;
while true
  kn = kn + dirn;
  if kn < 1 || kn > n, return; end
  if ~isfinite(Jp(kn)), continue; end
  if Jp(kn) >= thr
    f  = (thr - Jp(k)) / (Jp(kn) - Jp(k));
    f  = min(max(f, 0), 1);
    Ec = 10^(log10(E(k)) + f*(log10(E(kn)) - log10(E(k))));
    return;
  end
  k = kn;
end
end

%% ========================================================================
function h = thickness_on(hspec, xq)
%THICKNESS_ON Ice thickness at the query points.
%   A table is interpolated linearly between its samples and HELD CONSTANT
%   outside them. Thickness gradients near a grounding line are steep and
%   sign-changing, so extrapolating one off the end of the measurements
%   would put made-up structure into the rigidity exactly where the fit is
%   most sensitive to it.
xq = xq(:);
if isa(hspec, 'function_handle')
  h = hspec(xq); h = h(:);
elseif isscalar(hspec)
  h = repmat(hspec, numel(xq), 1);
else
  T = hspec;
  assert(size(T,2) == 2, 'opts.h as a table must be Nh x 2 of [x h]');
  [xs, is] = sort(T(:,1)); hs = T(is,2);
  keep = isfinite(xs) & isfinite(hs);
  xs = xs(keep); hs = hs(keep);
  assert(numel(xs) >= 2, 'opts.h table needs at least 2 finite rows');
  h = interp1(xs, hs, min(max(xq, xs(1)), xs(end)), 'linear');
end
assert(all(isfinite(h)) && all(h > 0), 'thickness must be finite and positive');
end
