function S = invertStrainRate(V, blk, opts)
%INVERTSTRAINRATE Legendre inversion for the vertical strain-rate profile.
%   S = INVERTSTRAINRATE(V, blk, opts) fits, independently in each
%   along-track block, the relative vertical velocity profile v(d) with a
%   Legendre expansion in normalised depth and differentiates it to get the
%   vertical strain rate.
%
%   MODEL. With x = 2*d/H - 1 mapping the fitted depth range onto [-1, 1],
%
%     v(d)      = sum_{k=0..K} c_k P_k(x)                     [m/yr]
%     eps_zz(d) = dv/dd = (2/H) sum_{k=0..K} c_k P_k'(x)      [1/yr]
%
%   For the default K = 2 this reduces to a strain rate that is linear in
%   depth,
%
%     eps_zz(x) = S1 + S2*x,   S1 = 2*c_1/H,   S2 = 6*c_2/H
%
%   so S1 is the depth-averaged vertical strain rate over the fitted range
%   and S2 is its top-to-bottom gradient. c_0 is a nuisance offset: it
%   absorbs any residual constant in the surface referencing and does not
%   affect eps_zz. This is the same parameterisation and the same reported
%   quantities (c0/c1/c2, S1, S2, epszz_mean, p_quad) as the earlier EGIG
%   2011-2012 vertical-strain products, so the two are directly comparable.
%
%   WEIGHTING AND DEGREES OF FREEDOM. Samples are weighted by inverse
%   variance where a per-sample sigma is available (V.v_std), otherwise by
%   coherence. The reported uncertainty and p_quad use an EFFECTIVE sample
%   count n_eff = n_used / opts.bins_per_look, because adjacent fast-time
%   bins in a multilooked interferogram are not independent - the fast-time
%   multilook window and the range resolution correlate them. Treating
%   every range bin as independent would inflate the degrees of freedom by
%   the look length and make p_quad meaninglessly small.
%
%   opts fields:
%     .order          Legendre order K (default 2)
%     .fit_top_depth  shallowest depth included in the fit [m]
%     .fit_bot_depth  deepest depth included [m]; [] = deepest valid
%     .norm_depth     H used to normalise depth [m]; [] = fitted bottom
%     .reg            ridge weight on the coefficients (default 0)
%     .bins_per_look  fast-time correlation length in bins (default 1)
%     .min_samples    minimum raw samples for a block to be inverted
%
%   Returns per-block row vectors unless noted:
%     S.coef        (K+1) x Nblk Legendre coefficients [m/yr]
%     S.coef_mmyr   the same in mm/yr (EGIG product units)
%     S.coef_std    (K+1) x Nblk 1-sigma of the coefficients [m/yr]
%     S.H           normalisation depth used [m]
%     S.top_depth, S.bot_depth   fitted depth range [m]
%     S.S1, S.S2    strain-rate mean and gradient [1/yr]
%     S.epszz_mean  depth-averaged vertical strain rate [1/yr] (= S1)
%     S.p_quad      two-sided p-value for the quadratic term c_2
%     S.depth_grid  Nz x Nblk reporting depth grid [m], per block
%     S.eps_zz      Nz x Nblk strain-rate profile [1/yr]
%     S.v_fit       Nz x Nblk fitted velocity profile [m/yr]
%     S.rms         weighted RMS residual [m/yr]
%     S.n_used      raw samples used
%     S.n_eff       effective independent samples
%
%   See also vdef.verticalDisplacement, vdef.legendreBasis.

K = opts.order;
Nblk = size(V.v, 2);
Nz = 201;

bins_per_look = max(1, opts.bins_per_look);

S = [];
S.coef       = nan(K+1, Nblk);
S.coef_std   = nan(K+1, Nblk);
S.H          = nan(1, Nblk);
S.top_depth  = nan(1, Nblk);
S.bot_depth  = nan(1, Nblk);
S.S1         = nan(1, Nblk);
S.S2         = nan(1, Nblk);
S.epszz_mean = nan(1, Nblk);
S.p_quad     = nan(1, Nblk);
S.rms        = nan(1, Nblk);
S.n_used     = zeros(1, Nblk);
S.n_eff      = zeros(1, Nblk);
S.eps_zz     = nan(Nz, Nblk);
S.v_fit      = nan(Nz, Nblk);
S.depth_grid = nan(Nz, Nblk);

for b = 1:Nblk
  d = V.depth(:,b);
  v = V.v(:,b);

  sel = isfinite(d) & isfinite(v) & d >= opts.fit_top_depth;
  if ~isempty(opts.fit_bot_depth)
    sel = sel & d <= opts.fit_bot_depth;
  end
  if nnz(sel) < opts.min_samples || nnz(sel) < (K+2)
    continue;
  end

  d = d(sel);
  v = v(sel);

  % Inverse-variance weights when available, coherence otherwise
  sd = V.v_std(sel,b);
  if any(isfinite(sd) & sd > 0)
    w = 1 ./ max(sd, eps).^2;
    w(~isfinite(w)) = 0;
  else
    w = blk.coh(sel,b);
    w(~isfinite(w)) = 0;
  end
  if ~any(w > 0)
    continue;
  end
  w = w / max(w);

  top = min(d);
  bot = max(d);
  if ~isempty(opts.norm_depth)
    H = opts.norm_depth;
  else
    H = bot;
  end
  if ~(H > 0)
    continue;
  end

  x = 2*d/H - 1;
  A = vdef.legendreBasis(x, K);

  W    = w(:);
  AtWA = A' * bsxfun(@times, W, A);
  AtWb = A' * (W .* v);

  if opts.reg > 0
    % Ridge on the coefficients, scaled to the problem's own conditioning
    lam = opts.reg * mean(diag(AtWA));
    AtWA = AtWA + lam * eye(K+1);
  end

  c = AtWA \ AtWb;

  resid = v - A*c;
  n_used = numel(v);
  n_eff  = max(n_used / bins_per_look, K+2);
  dof    = max(n_eff - (K+1), 1);

  % Variance scale factor: weighted SSE per effective degree of freedom.
  % Down-weighting W by bins_per_look to reflect the true information
  % content cancels against inv(A'WA), so the correlation correction acts
  % entirely through dof - which is the point.
  Cov  = (sum(W .* resid.^2) / dof) * (AtWA \ eye(K+1));
  cstd = sqrt(abs(diag(Cov)));

  % Weighted RMS residual, for reporting
  rms_w = sqrt(sum(W .* resid.^2) / sum(W));

  S.coef(:,b)     = c;
  S.coef_std(:,b) = cstd;
  S.H(b)          = H;
  S.top_depth(b)  = top;
  S.bot_depth(b)  = bot;
  S.rms(b)        = rms_w;
  S.n_used(b)     = n_used;
  S.n_eff(b)      = n_eff;

  if K >= 1
    S.S1(b) = 2*c(2)/H;
  end
  if K >= 2
    S.S2(b) = 6*c(3)/H;
    % Two-sided Student-t p-value for c_2, via the exact incomplete-beta
    % form (betainc is base MATLAB and Octave; tcdf would need a toolbox)
    if cstd(3) > 0
      t = c(3)/cstd(3);
      S.p_quad(b) = betainc(dof/(dof + t^2), dof/2, 0.5);
    end
  end

  % Depth-averaged strain rate over the fitted range: (v(bot)-v(top))/(bot-top)
  xg = linspace(2*top/H - 1, 2*bot/H - 1, Nz).';
  [Pg, dPg] = vdef.legendreBasis(xg, K);
  S.depth_grid(:,b) = (xg + 1) * H / 2;
  S.v_fit(:,b)      = Pg * c;
  S.eps_zz(:,b)     = (2/H) * (dPg * c);
  S.epszz_mean(b)   = (S.v_fit(end,b) - S.v_fit(1,b)) / (bot - top);
end

end
