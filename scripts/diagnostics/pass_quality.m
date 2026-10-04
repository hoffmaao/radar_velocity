function OUT = pass_quality(opts)
%PASS_QUALITY Which PASSES are bad, not just which pairs.
%
%   THE GAP THIS FILLS. The chain already gates PAIRS: coalignPair refuses
%   a pair whose surface-profile correlation falls below
%   coalign_min_quality, the baseline gate drops pairs wider than
%   max_baseline, and vdef.surfaceAdmittance rejects a pass whose GPS line
%   mean departs from the tide. What nothing does is ask whether the
%   rejected pairs SHARE A PASS. Nine of GL4's pairs failed coalignment at
%   correlations of 0.09 to 0.14 while the rest sat above 0.87; if those
%   nine all involve one pass, the honest description is one bad pass, not
%   nine bad pairs, and the remedy is to drop it rather than to keep
%   discovering it pair by pair. This function builds the full pass-by-pass
%   matrix and takes its marginals so that question has an answer.
%
%   WHAT IS MEASURED, PER PASS, WITHOUT PAIRING
%     surf_dB     median surface-return power - a weak return means little
%                 came back, from attitude, burial, or a bad transmit
%     width_bin   median -6 dB width of the surface return in fast-time
%                 bins. This is the CHIRP/compression proxy: a degraded or
%                 mismatched chirp compresses to a broader peak at the same
%                 energy, so width rises while surf_dB need not fall
%     snr_dB      median ice-column power over the noise estimate taken
%                 above the surface, in the reported depth window
%     cover       fraction of columns with a usable surface reference
%     gps_resid   line-mean height residual about the CATS2008 fit [m],
%                 the quantity vdef.surfaceAdmittance gates on
%     dy_med      median cross-track offset from the main pass [m], and
%                 dy_rng its along-track range - a pass that wanders
%                 off-track images different ice
%
%   WHAT IS MEASURED, PER PAIR, THEN MARGINALISED
%     quality     coalignPair's surface-profile correlation
%     applied     whether the shift survived the quality floor
%     coh         median interferometric coherence in the depth window
%     baseline    max |cross-track separation| [m]
%   For each pass the marginals are the MEDIAN over its partners, so one
%   bad partner cannot condemn a good pass. A pass that is bad shows a low
%   median against everything.
%
%   CALIBRATION AGAINST THE MAIN PASS. multipass calibrates every pass
%   relative to the main pass - a fast-time shift (coregistration) and a
%   complex scalar (equalization) - so whether that worked is measured
%   against the main pass, per pass, with FIXED targets rather than robust
%   z, because the right answer is known: zero.
%     coh_ns      multipass's coherent criterion re-run on the product: the
%                 shift maximising |sum(shifted .* conj(main))| over the
%                 coreg window, surface included [ns, positive = pass
%                 later]. Exactly what the coreg stage optimised, recovers
%                 0.98-1.04 of an injected shift; the alignment check.
%                 Target 0
%     coh_col_ns  the same criterion over the ice column only. RECORDED
%                 ONLY: it is biased by the deformation signal itself. A
%                 coherent sum over a column whose phase ramps with depth
%                 moves its optimum, so on GL3 it read +0.26..+0.40 ns on
%                 passes days from the main pass, made of -0.11..-0.24 ns in
%                 the upper column and +0.51..+0.78 ns in the deep column -
%                 opposite signs, which no misalignment can produce, growing
%                 with time from the main pass as the deep phase does
%     phase_deg   surface-referenced interferometric phase of the upper
%                 column (10-150 bins below the surface): a geometric shift
%                 tau of the layers relative to the surface would put
%                 2*pi*fc*tau here (1 deg = 3.7 ps at 750 MHz), so this is
%                 the independent, picosecond check on registration
%     layer_ns    power-profile ("envelope") delay of the layers: log
%                 trace-averaged power, minus a 41-bin running mean,
%                 Hann-tapered, 64x upsampled. RECORDED ONLY. It recovers
%                 0.94-1.12 of an injected shift, but between DIFFERENT
%                 passes it also responds to profile shape: on the
%                 surface-coupled GL2/GL3/GL4 builds it read up to +/-0.5 ns
%                 in a date-ordered pattern that neither the column coherent
%                 criterion (<= 0.1 ns) nor the phase (<= 11 deg where 40-130
%                 would be implied) shows
%     coalign_ns  coalignPair's surface-envelope estimate, RECORDED ONLY.
%                 On the surface-coupled GL1 build it read -0.16 to -0.55 ns
%                 where both estimators above read |x| < 0.09 ns, and it
%                 recovered 85% of an injected shift from a -0.29 ns start:
%                 it is biased by surface-return shape, so it cannot judge
%                 an alignment finer than its own bias
%     drift_ns    along-track spread of the layer envelope delay over 10
%                 windows, p95 - p5; target 0. A scalar shift cannot remove
%                 a drift. The envelope estimator is phase-independent and
%                 its shape bias is common to a pass's windows, so it
%                 cancels in a spread; across GL1 passes its windows share a
%                 0.17 ns pattern with 0.16 ns per-pass noise, where the
%                 coherent criterion in the same windows shares 0.68 ns with
%                 0.54 ns (a bias from the common main pass's local phase)
%     leak_r, leak_slope
%                 correlation and regression slope of the windowed layer
%                 delay against -(ref_z_k(x) - ref_z_main(x))/(c/2), the
%                 shape the z-motion compensation imposes. On a floating
%                 shelf ref_z is the tide, so a product built WITH the
%                 compensation shows slope ~1 here; a correctly built one
%                 shows none
%     selftest    per product, a known +0.2-bin shift (multipass's own sinc
%                 convention) is injected into the pass most coherent with
%                 the main pass; [layer, coherent full, coherent column]
%                 recovery, and the two that are flagged on (layer for
%                 drift/leak, coherent full for alignment) must be within
%                 15% - the
%                 audit checks its own instruments before believing them.
%                 selftest_worst repeats it on the least coherent pass, as
%                 the estimators' worst case, without gating
%     gain_dB     ice-column power relative to the main pass (10 bins below
%                 the surface down to ~3 us), the quantity the
%                 surface-coupled build equalizes on; target 0
%     srf_dB      the same for the surface return, for the record - it also
%                 carries surface conditions, not only system gain
%     gam_dB      20 log10 of the pass's global coherence with the main pass
%                 over the ice column. multipass's own equalization estimate
%                 is gain_dB + gam_dB, which is why it inflates a
%                 decorrelated pass
%     eq_deg      phase of the mean interferogram with the main pass over
%                 every bin, exactly as multipass computes its equalization:
%                 in an equalized product this is the phase closure, target 0
%     wander_deg  along-track circular spread of the surface interferometric
%                 phase. A constant equalization phase cannot flatten a
%                 phase that wanders; target small
%
%   HOW THE FLAGS WORK. Each per-pass metric is scored by robust z, the
%   deviation from the median in units of 1.4826 times the median absolute
%   deviation, computed within a product because absolute levels differ
%   between legs. A pass is FLAGGED when any metric is worse than
%   opts.z_flag robust sigmas AND that metric is also worse than a fixed
%   floor, so a tight distribution cannot manufacture an outlier out of a
%   pass that is fine in absolute terms. Both conditions are reported, so
%   a marginal call can be overridden by eye.
%
%   THIS DOES NOT DROP ANYTHING. It reports. Excluding a pass changes every
%   downstream number, so that decision belongs in the driver that quotes
%   the result, not in a diagnostic.
%
%   The calibration metrics are flagged separately (cal_flag, cal_why)
%   against the fixed limits in opts.cal_limits.
%
%   opts, all optional: .products .mp_dir .zsel .block .z_flag .out_dir
%   .tag (suffix for pass_quality<tag>.mat, so auditing another build does
%   not overwrite the file the tidal-stack pass dropping reads)
%   .cal_limits
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "addpath('<code>/scripts/diagnostics'); pass_quality"
%
%   See also vdef.coalignPair, vdef.surfaceAdmittance, pass_tide.

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
addpath(root); addpath(here); addpath(fullfile(root, 'opr_vvel'));
def = struct('products', {vdef.surveyLines()}, ...
  'mp_dir', '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass', ...
  'zsel', 20:10:250, 'z_flag', 3.5, 'out_dir', vdef.figureDir(), 'tag', '', ...
  'cal_limits', struct('align_ns', 0.25, 'phase_deg', 30, 'drift_ns', 1.0, 'leak_r', 0.5, ...
                       'leak_slope', 0.3, 'eq_deg', 5, 'gain_dB', 1.0, 'wander_deg', 30));
fn = fieldnames(def);
for i = 1:numel(fn)
  if ~isfield(opts, fn{i}) || (isempty(opts.(fn{i})) && ~ischar(def.(fn{i}))), opts.(fn{i}) = def.(fn{i}); end
end
C = vdef.constants();
P = vdef.firnColumn(vdef.defaultParams());

OUT = struct('name', {}, 'seg', {}, 'tmid', {}, 'main_idx', {}, 'surf_dB', {}, ...
  'width_bin', {}, 'snr_dB', {}, 'cover', {}, 'gps_resid', {}, 'dy_med', {}, ...
  'dy_rng', {}, 'Q', {}, 'A', {}, 'COH', {}, 'BL', {}, 'q_med', {}, 'frac_applied', {}, ...
  'coh_med', {}, 'flag', {}, 'why', {}, 'zmotion', {}, 'coh_col_ns', {}, 'coh_ns', {}, ...
  'phase_deg', {}, 'layer_ns', {}, 'coalign_ns', {}, 'selftest', {}, 'selftest_worst', {}, 'drift_ns', {}, ...
  'leak_r', {}, 'leak_slope', {}, 'gain_dB', {}, 'srf_dB', {}, 'gam_dB', {}, 'eq_deg', {}, ...
  'wander_deg', {}, 'cal_flag', {}, 'cal_why', {});
for n = 1:numel(opts.products)
  pn = opts.products{n}; nm = strrep(pn, 'EAGER_2022_', '');
  fn_mp = fullfile(opts.mp_dir, sprintf('%s_multipass03.mat', pn));
  if ~exist(fn_mp, 'file'), fprintf('%s: missing product\n', nm); continue; end
  param = vvel_defaults(struct('vvel', struct('pass_name', pn))); o = param.vvel;
  fprintf('\n===== %s =====\n', nm); t0 = tic;
  D = load(fn_mp, 'data', 'pass', 'param_multipass');
  pm = D.param_multipass.multipass;
  if isfield(pm, 'pass_en_mask') && ~isempty(pm.pass_en_mask), en = find(pm.pass_en_mask);
  else, en = 1:numel(D.pass); end
  main_idx = pm.baseline_master_idx; mp_ = D.pass(main_idx);
  Time = mp_.time(:); Nt = numel(Time); Nx = size(D.data, 2); Np = numel(D.pass);
  if isfield(mp_, 'surface') && ~isempty(mp_.surface), Surface = mp_.surface(:).';
  else, Surface = mp_.layers(1).twtt_ref(:).'; end
  Surface = Surface(1:Nx);
  if isfield(mp_, 'wfs') && isfield(mp_.wfs(1), 'fc') && ~isempty(mp_.wfs(1).fc), fc = mp_.wfs(1).fc;
  else, fc = 750e6; end
  depth = vdef.depthFromTwtt(P, Time - mean(Surface, 'omitnan'));
  zi = arrayfun(@(z) find(depth >= z, 1, 'first'), opts.zsel, 'uni', 0); zi = [zi{:}];
  ref_bin = round(interp1(Time, 1:Nt, Surface, 'linear', NaN));
  [~, pass_ok_gps, info] = pass_tide(pn, opts.mp_dir);

  %% per-pass, single-pass metrics
  seg = repmat({''},1,Np); tmid = nan(1,Np);
  surf_dB = nan(1,Np); width_bin = nan(1,Np); snr_dB = nan(1,Np); cover = nan(1,Np);
  dy_med = nan(1,Np); dy_rng = nan(1,Np);
  ry_main = D.pass(main_idx).ref_y(:).';
  for k = 1:Np
    if isfield(D.pass(k),'param_pass') && isfield(D.pass(k).param_pass,'day_seg')
      seg{k} = char(D.pass(k).param_pass.day_seg);
    end
    tmid(k) = mean(D.pass(k).gps_time, 'omitnan');
    ry = D.pass(k).ref_y(:).'; m = min(numel(ry), numel(ry_main));
    dy = ry(1:m) - ry_main(1:m);
    dy_med(k) = median(dy, 'omitnan'); dy_rng(k) = range(dy(isfinite(dy)));
    kk = find(en == k, 1);
    if isempty(kk), continue; end
    s = D.data(:,:,kk);
    pw = abs(s).^2;
    good = isfinite(ref_bin) & ref_bin > 25 & ref_bin < Nt-5;
    cover(k) = mean(good);
    cols = find(good); if isempty(cols), continue; end
    cols = cols(1:max(1,round(numel(cols)/400)):end);      % subsample for speed
    pk = nan(1,numel(cols)); wd = nan(1,numel(cols)); nf = nan(1,numel(cols));
    for c = 1:numel(cols)
      b = ref_bin(cols(c));
      seg_pw = pw(max(1,b-20):min(Nt,b+20), cols(c));
      [pmax, imax] = max(seg_pw);
      if ~isfinite(pmax) || pmax <= 0, continue; end
      pk(c) = pmax;
      % -6 dB width in bins, contiguous about the peak
      thr = pmax/4; lo = imax; hi = imax;
      while lo > 1 && seg_pw(lo-1) >= thr, lo = lo-1; end
      while hi < numel(seg_pw) && seg_pw(hi+1) >= thr, hi = hi+1; end
      wd(c) = hi - lo + 1;
      nf(c) = median(pw(1:max(2,b-40), cols(c)), 'omitnan');   % above the surface
    end
    surf_dB(k)   = 10*log10(median(pk, 'omitnan'));
    width_bin(k) = median(wd, 'omitnan');
    colpw = median(pw(zi, cols), 1, 'omitnan');
    snr_dB(k) = 10*log10(median(colpw,'omitnan') / max(median(nf,'omitnan'), eps));
    clear s pw
  end
  gps_resid = info.resid;

  %% calibration against the main pass: gain, coherence, phase closure, wander
  km = find(en == main_idx, 1);
  sm = double(D.data(:,:,km));   % double: the coherent sums below refine a peak flatter than single precision resolves
  dt = Time(2) - Time(1);
  sb = round(median(ref_bin, 'omitnan'));
  col_rows = (sb + 10) : min(Nt, sb + round(3.0e-6/dt));   % as the equalize stage
  srf_rows = max(1, sb - 3) : min(Nt, sb + 3);
  gain_dB = nan(1,Np); srf_dB = nan(1,Np); gam_dB = nan(1,Np); eq_deg = nan(1,Np);
  wander_deg = nan(1,Np);
  Pm_col = mean(abs(sm(col_rows,:)).^2, 'all', 'omitnan');
  Pm_srf = mean(abs(sm(srf_rows,:)).^2, 'all', 'omitnan');
  for k = 1:Np
    kk = find(en == k, 1);
    if isempty(kk), continue; end
    s = D.data(:,:,kk);
    x = s .* conj(sm);
    eq_deg(k) = angle(mean(x(isfinite(x)))) * 180/pi;      % multipass's estimate, all bins
    xc = x(col_rows,:); pc = abs(s(col_rows,:)).^2; pmc = abs(sm(col_rows,:)).^2;
    g = isfinite(xc) & isfinite(pc) & isfinite(pmc);
    gain_dB(k) = 10*log10(mean(pc(g)) / Pm_col);
    gam_dB(k) = 20*log10(abs(mean(xc(g))) / sqrt(mean(pc(g)) * mean(pmc(g))));
    srf_dB(k) = 10*log10(mean(abs(s(srf_rows,:)).^2, 'all', 'omitnan') / Pm_srf);
    xs = sum(x(srf_rows,:), 1); xs = xs(isfinite(xs) & abs(xs) > 0);
    if ~isempty(xs)
      Rl = abs(mean(xs ./ abs(xs)));
      wander_deg(k) = sqrt(-2*log(max(Rl, eps))) * 180/pi;
    end
    clear s x xc pc pmc
  end
  % fast-time alignment against the main pass: coherent criterion over the
  % ice column (decisive), over the coreg window (closure), the upper
  % column's surface-referenced phase (geometry), and the layer envelope
  % (record)
  coh_col_ns = nan(1,Np); coh_ns = nan(1,Np); phase_deg = nan(1,Np); layer_ns = nan(1,Np);
  coalign_ns = nan(1,Np); drift_ns = nan(1,Np); leak_r = nan(1,Np); leak_slope = nan(1,Np);
  coh_col_ns(main_idx) = 0; coh_ns(main_idx) = 0; phase_deg(main_idx) = 0; layer_ns(main_idx) = 0;
  coalign_ns(main_idx) = 0; drift_ns(main_idx) = 0;
  coreg_rows = max(1, sb - 20) : min(Nt, sb + round(3.0e-6/dt));
  lo = max(1, min(coreg_rows) - 8); hi = min(Nt, max(coreg_rows) + 8);   % only these rows are ever shifted
  rf = coreg_rows - lo + 1; rc = col_rows - lo + 1;
  ccols = 1:4:Nx; sgrid = -0.4:0.01:0.4;
  nw = 10; we = round(linspace(1, Nx+1, nw+1));
  smb = sm(lo:hi, :);
  sbc = round(interp1(Time, 1:Nt, Surface, 'linear', NaN));             % per-column surface bin
  rz_m = D.pass(main_idx).ref_z(:);
  for k = 1:Np
    kk = find(en == k, 1);
    if isempty(kk) || k == main_idx, continue; end
    s = double(D.data(:,:,kk)); sbk = s(lo:hi, :);
    coh_col_ns(k) = coherent_shift(sbk(:, ccols), smb(:, ccols), rc, sgrid) * dt*1e9;
    coh_ns(k) = coherent_shift(sbk(:, ccols), smb(:, ccols), rf, sgrid) * dt*1e9;
    layer_ns(k) = layer_delay(s, sm, col_rows, 1:Nx) * dt*1e9;
    acc = 0;
    for c = 1:Nx
      b = sbc(c); if ~isfinite(b) || b < 3 || b + 150 > Nt, continue; end
      r = sum(s(b-2:b+2, c) .* conj(sm(b-2:b+2, c))); if abs(r) == 0, continue; end
      acc = acc + sum(s(b+10:b+150, c) .* conj(sm(b+10:b+150, c))) * conj(r)/abs(r);
    end
    phase_deg(k) = angle(acc) * 180/pi;
    wl = nan(1,nw); tw = nan(1,nw);
    rz_k = D.pass(k).ref_z(:); nn = min([numel(rz_k) numel(rz_m) Nx]);
    tpl = -(rz_k(1:nn) - rz_m(1:nn)) / (C.c/2) * 1e9;     % ns, positive = pass later
    for w = 1:nw
      cw = we(w):min(we(w+1)-1, nn);
      if numel(cw) < 20, continue; end
      wl(w) = layer_delay(s, sm, col_rows, cw) * dt*1e9;
      tw(w) = mean(tpl(cw), 'omitnan');
    end
    gw = isfinite(wl) & isfinite(tw);
    if nnz(gw) >= 4, drift_ns(k) = diff(prctile(wl(gw), [5 95])); end
    if nnz(gw) >= 5 && std(tw(gw)) > 0
      cc = corrcoef(tw(gw), wl(gw)); leak_r(k) = cc(1,2);
      pf = polyfit(tw(gw), wl(gw), 1); leak_slope(k) = pf(1);
    end
  end
  % the audit's own instruments: a known shift must come back. The GATING
  % test uses the pass most coherent with the main pass, because it tests
  % the estimator, not the pass - on a decorrelated pass a column-only
  % coherent sum is weak (GL2's 20221207_03, coherence 0.24, returned 0.83
  % of the injection). The least coherent pass is run too and reported as
  % the worst case, without gating.
  cand = en(en ~= main_idx); [~, ib] = max(gam_dB(cand)); [~, iw] = min(gam_dB(cand));
  selftest = injection_test(D.data(:,:,find(en == cand(ib), 1)), sm, smb, lo, hi, ccols, rf, rc, col_rows, sgrid);
  selftest_worst = injection_test(D.data(:,:,find(en == cand(iw), 1)), sm, smb, lo, hi, ccols, rf, rc, col_rows, sgrid);
  clear s sbk sm smb
  fprintf(['  estimator self-test (inject +0.20 bin) on the most coherent pass %d: layer %.2f, coherent full %.2f, ' ...
           'coherent column %.2f%s\n  worst case, least coherent pass %d: layer %.2f, full %.2f, column %.2f\n'], ...
    cand(ib), selftest, repmat(' *** AUDIT ESTIMATORS FAIL - TIME FLAGS NOT TRUSTWORTHY ***', 1, ...
    any(abs(selftest(1:2) - 1) > 0.15)), cand(iw), selftest_worst);

  %% pair matrices
  Q = nan(Np); A = false(Np); COH = nan(Np); BL = nan(Np);
  t1 = tic; npair = 0;
  for i = 1:Np-1
    for j = i+1:Np
      ki = find(en == i, 1); kj = find(en == j, 1);
      if isempty(ki) || isempty(kj), continue; end
      by = D.pass(j).ref_y(:) - D.pass(i).ref_y(:);
      BL(i,j) = max(abs(by(isfinite(by)))); BL(j,i) = BL(i,j);
      s_ref = D.data(:,:,ki); s_sec = D.data(:,:,kj);
      [s_sec, ca] = vdef.coalignPair(s_ref, s_sec, ...
        struct('Time', Time, 'Surface', Surface, 'fc', fc), o);
      Q(i,j) = ca.quality; Q(j,i) = Q(i,j);
      A(i,j) = ca.applied;  A(j,i) = A(i,j);
      if i == main_idx || j == main_idx
        % coalignPair reports the secondary (j) relative to the reference
        % (i); recorded as "pass k later than main", not flagged on
        k = i + j - main_idx;
        coalign_ns(k) = (1 - 2*(j == main_idx)) * ca.dtau_bulk * 1e9;
      end
      [~, coh] = vdef.multilook(s_ref, s_sec, o.mlook_window);
      cv = coh(zi, :); COH(i,j) = median(cv(:), 'omitnan'); COH(j,i) = COH(i,j);
      clear s_ref s_sec coh cv
      npair = npair + 1;
    end
  end
  clear D
  fprintf('  %d passes, %d pairs formed in %.0f s\n', Np, npair, toc(t1));

  %% marginals and flags
  q_med = nan(1,Np); frac_applied = nan(1,Np); coh_med = nan(1,Np);
  for k = 1:Np
    q_med(k) = median(Q(k,[1:k-1 k+1:Np]), 'omitnan');
    a = A(k,[1:k-1 k+1:Np]); qq = Q(k,[1:k-1 k+1:Np]);
    frac_applied(k) = mean(a(isfinite(qq)));
    coh_med(k) = median(COH(k,[1:k-1 k+1:Np]), 'omitnan');
  end
  % robust z against the product's own distribution, paired with an
  % absolute floor so a tight spread cannot invent an outlier
  tests = { 'coalign',  q_med,        -1, o.coalign_min_quality; ...
            'coherence',coh_med,      -1, 0.30; ...
            'applied',  frac_applied, -1, 0.70; ...
            'surf_dB',  surf_dB,      -1, -Inf; ...
            'width',    width_bin,    +1,  Inf; ...
            'snr_dB',   snr_dB,       -1, -Inf; ...
            'coverage', cover,        -1, 0.80; ...
            'gps',      abs(gps_resid),+1, 0.25 };
  flag = false(1,Np); why = repmat({''},1,Np);
  for t = 1:size(tests,1)
    v = tests{t,2}; sgn = tests{t,3}; floor_v = tests{t,4};
    g = isfinite(v); if nnz(g) < 4, continue; end
    md = median(v(g)); sd = 1.4826*median(abs(v(g)-md));
    if sd <= 0, continue; end
    z = sgn*(v - md)/sd;                       % positive z = worse
    bad = g & z > opts.z_flag;
    if isfinite(floor_v)
      if sgn < 0, bad = bad & (v < floor_v); else, bad = bad & (v > floor_v); end
    end
    for k = find(bad)
      flag(k) = true;
      why{k} = strtrim([why{k} ' ' sprintf('%s(z%+.1f)', tests{t,1}, z(k))]);
    end
  end

  fprintf('  %3s %-14s %7s %7s %7s %6s %7s %7s %7s %6s  %s\n', 'k', 'day_seg', ...
    'q_med', 'coh', 'applied', 'cover', 'surf_dB', 'width', 'snr_dB', 'gps_m', 'flags');
  for k = 1:Np
    star = ' '; if k == main_idx, star = '*'; end
    fprintf('  %2d%s %-14s %7.3f %7.3f %7.2f %6.2f %7.1f %7.1f %7.1f %6.2f  %s\n', ...
      k, star, seg{k}, q_med(k), coh_med(k), frac_applied(k), cover(k), ...
      surf_dB(k), width_bin(k), snr_dB(k), gps_resid(k), why{k});
  end
  if any(flag)
    fprintf('  FLAGGED: %s\n', strjoin(arrayfun(@(k) sprintf('%d (%s)', k, seg{k}), ...
      find(flag), 'uni', 0), ', '));
  else
    fprintf('  no pass flagged at %.1f robust sigma with an absolute floor\n', opts.z_flag);
  end
  if ~all(pass_ok_gps)
    fprintf('  (the GPS gate in vdef.surfaceAdmittance already rejects: %s)\n', ...
      mat2str(find(~pass_ok_gps)));
  end

  % calibration flags: fixed targets, not robust z (the answer is zero)
  L = opts.cal_limits; cal_flag = false(1,Np); cal_why = repmat({''},1,Np);
  ctests = { 'align',  abs(coh_ns) > L.align_ns; ...
             'geom',   abs(phase_deg) > L.phase_deg; ...
             'drift',  drift_ns > L.drift_ns; ...
             'leak',   abs(leak_r) > L.leak_r & abs(leak_slope) > L.leak_slope; ...
             'eqphase',abs(eq_deg) > L.eq_deg; ...
             'gain',   abs(gain_dB) > L.gain_dB; ...
             'wander', wander_deg > L.wander_deg };
  for t = 1:size(ctests,1)
    for k = find(ctests{t,2})
      if k == main_idx, continue; end
      cal_flag(k) = true; cal_why{k} = strtrim([cal_why{k} ' ' ctests{t,1}]);
    end
  end
  zmotion = vdef.zmotionApplied(pm);
  fprintf('\n  calibration against main pass %d (z-motion compensation %s):\n', main_idx, ...
    char(string(zmotion)));
  fprintf('  %3s %-14s %7s %7s %7s %8s %8s %7s %7s %7s %7s %7s %7s %7s %7s  %s\n', 'k', 'day_seg', ...
    '(col)', 'coh_ns', 'phase', '(layer)', '(coalgn)', 'drift', 'leak_r', 'slope', 'gain_dB', ...
    'srf_dB', 'gam_dB', 'eq_deg', 'wander', 'cal flags');
  for k = 1:Np
    star = ' '; if k == main_idx, star = '*'; end
    fprintf('  %2d%s %-14s %7.3f %7.3f %7.1f %8.3f %8.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.1f %7.1f  %s\n', ...
      k, star, seg{k}, coh_col_ns(k), coh_ns(k), phase_deg(k), layer_ns(k), coalign_ns(k), ...
      drift_ns(k), leak_r(k), leak_slope(k), gain_dB(k), srf_dB(k), gam_dB(k), eq_deg(k), ...
      wander_deg(k), cal_why{k});
  end
  if any(cal_flag)
    fprintf('  CALIBRATION FLAGGED: %s\n', strjoin(arrayfun(@(k) sprintf('%d (%s: %s)', k, ...
      seg{k}, cal_why{k}), find(cal_flag), 'uni', 0), ', '));
  else
    fprintf('  every pass meets the calibration limits\n');
  end

  OUT(end+1) = struct('name', pn, 'seg', {seg}, 'tmid', tmid, 'main_idx', main_idx, ...
    'surf_dB', surf_dB, 'width_bin', width_bin, 'snr_dB', snr_dB, 'cover', cover, ...
    'gps_resid', gps_resid, 'dy_med', dy_med, 'dy_rng', dy_rng, 'Q', Q, 'A', A, ...
    'COH', COH, 'BL', BL, 'q_med', q_med, 'frac_applied', frac_applied, ...
    'coh_med', coh_med, 'flag', flag, 'why', {why}, 'zmotion', zmotion, 'coh_col_ns', coh_col_ns, ...
    'coh_ns', coh_ns, 'phase_deg', phase_deg, 'layer_ns', layer_ns, 'coalign_ns', coalign_ns, ...
    'selftest', selftest, 'selftest_worst', selftest_worst, ...
    'drift_ns', drift_ns, 'leak_r', leak_r, 'leak_slope', leak_slope, 'gain_dB', gain_dB, ...
    'srf_dB', srf_dB, 'gam_dB', gam_dB, 'eq_deg', eq_deg, 'wander_deg', wander_deg, ...
    'cal_flag', cal_flag, 'cal_why', {cal_why}); %#ok<AGROW>
  fprintf('  %s done in %.0f s\n', nm, toc(t0));
end
if ~exist(opts.out_dir,'dir'), mkdir(opts.out_dir); end
fn_out = fullfile(opts.out_dir, sprintf('pass_quality%s.mat', opts.tag));
save(fn_out, 'OUT');
fprintf('\nsaved %s\n', fn_out);
fprintf(['q_med/coh/applied are medians over a pass''s partners, so one bad partner cannot\n' ...
         'condemn a good pass. width is the -6 dB width of the compressed surface return in\n' ...
         'bins - the chirp proxy, which rises when compression degrades even if power does not\n' ...
         'fall. Nothing here is dropped: excluding a pass is the driver''s decision.\n']);
end

%% ========================================================================
function y = mp_shift(x, s)
%MP_SHIFT multipass's 13-tap sinc fractional shift (multipass.m, "Pass: 3"),
%   positive s moves the image toward the radar; NaN counts as zero, as the
%   nansum there makes it. Vectorised over rows, identical arithmetic.
x(~isfinite(x)) = 0; wh = round(s); fr = wh - s; N = size(x, 1);
y = zeros(size(x), 'like', x);
for m = -6:6
  src = (1:N) + wh + m; ok = src >= 1 & src <= N;
  y(ok,:) = y(ok,:) + sinc(m + fr) * x(src(ok),:);
end
end

%% ========================================================================
function d = layer_delay(a, b, rows, cols)
%LAYER_DELAY Delay of image a relative to b in bins (positive = a later),
%   from the internal layers: log trace-averaged power over `rows`, minus
%   a 41-bin running mean so the layers and not the depth decay drive the
%   correlation (without it a raw profile reads ~60% of a known shift),
%   Hann-tapered, 64x FFT-upsampled, parabolic peak. Recovers 95-100% of
%   injected shifts on the GL1 product, precision ~0.02 bin.
up = 64;
pa = log10(mean(abs(a(rows,cols)).^2, 2, 'omitnan'));
pb = log10(mean(abs(b(rows,cols)).^2, 2, 'omitnan'));
pa = pa - movmean(pa, 41); pb = pb - movmean(pb, 41);
n = numel(pa); w = 0.5 - 0.5*cos(2*pi*(0:n-1).'/(n-1));
pa = (pa - mean(pa)).*w; pb = (pb - mean(pb)).*w;
ua = interpft(pa, up*n); ub = interpft(pb, up*n);
L = 3*up; c = zeros(2*L+1, 1);
for q = -L:L
  if q >= 0, c(q+L+1) = sum(ua(1+q:end) .* ub(1:end-q));
  else,      c(q+L+1) = sum(ua(1:end+q) .* ub(1-q:end)); end
end
[~, i] = max(c); i = min(max(i, 2), numel(c)-1);
dd = 0.5*(c(i-1) - c(i+1)) / (c(i-1) - 2*c(i) + c(i+1));
d = ((i - L - 1) + dd) / up;
end

%% ========================================================================
function s = coherent_shift(a, b, rows, grid)
%COHERENT_SHIFT multipass's coregistration criterion re-run: the shift
%   toward the radar that maximises |sum(shift(a) .* conj(b))| over `rows`
%   (positive = a is later and must move earlier), sub-grid by parabola.
%   a and b are already cut to the rows that can be reached (a margin
%   around `rows`) and to the columns wanted. NaN at the grid edge, which
%   means the search window was too small.
bb = conj(b(rows, :)); v = zeros(size(grid));
for g = 1:numel(grid)
  y = mp_shift(a, grid(g)); x = y(rows,:) .* bb; v(g) = abs(sum(x(isfinite(x))));
end
[~, i] = max(v);
if i == 1 || i == numel(v), s = NaN; return; end
dd = 0.5*(v(i-1) - v(i+1)) / (v(i-1) - 2*v(i) + v(i+1));
s = grid(i) + dd*(grid(2) - grid(1));
end

%% ========================================================================
function r = injection_test(s0, sm, smb, lo, hi, ccols, rf, rc, col_rows, sgrid)
%INJECTION_TEST Fraction of a known +0.2-bin shift (toward the radar) that
%   each estimator recovers: [layer envelope, coherent full, coherent column].
s0 = double(s0); s1 = mp_shift(s0, 0.2); Nx = size(s0, 2);
b0 = s0(lo:hi, ccols); b1 = s1(lo:hi, ccols); bm = smb(:, ccols);
r = [ (layer_delay(s1, sm, col_rows, 1:Nx) - layer_delay(s0, sm, col_rows, 1:Nx)) / -0.2, ...
      (coherent_shift(b1, bm, rf, sgrid) - coherent_shift(b0, bm, rf, sgrid)) / -0.2, ...
      (coherent_shift(b1, bm, rc, sgrid) - coherent_shift(b0, bm, rc, sgrid)) / -0.2 ];
end
