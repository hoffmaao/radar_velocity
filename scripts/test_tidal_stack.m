%TEST_TIDAL_STACK Unit test for vdef.tidalStack and vdef.complexBlocks.
%   The claim the estimator exists to make good: the tidal column response
%   can be read off the surface-referenced COMPLEX interferograms of all
%   pairs at once, by counter-rotating each pair by the model phase and
%   summing, with no unwrapping anywhere - and the result is honest about
%   what it cannot separate.
%
%   Covers, on synthetic pair fields shaped like the EAGER lines (60 pairs,
%   2 m tide range, k(100 m) = -54 rad/m, 0.02 rad of phase noise on each
%   block-mean field - the noise of a 200-column coherence-weighted mean,
%   not of a pixel; at 0.3 rad the tide's 0.06 rad swing is invisible and
%   the scan's sd is 3 mm/m, which is the precision statement, not a bug):
%     1. a known a(z, x) in -4..+4 mm/m is recovered within a grid step or
%        two at every cell, with a high coherent peak
%     2. the self-test: a response injected as a phase rotation returns as
%        an EXACT shift of the estimate (this is what pins the phase sign)
%     3. the bootstrap sd covers the actual error
%     4. an artefact term g*d_res with d_res UNCORRELATED with the tide is
%        separated by the 2-D scan: both a and g recovered, rcol ~ 0
%     5. the same artefact with d_res COLLINEAR with the tide is reported
%        as such (rcol > 0.9) rather than silently attributed to a
%     6. vdef.complexBlocks references each column by phase only and gives
%        zero weight to columns with no reference or below-threshold
%        coherence
%     7. an error that belongs to a PASS enters every pair of that pass:
%        the delete-one-pass jackknife tracks the resulting scatter of a
%        and the pair bootstrap understates it
%     8. a secular trend whose pair intervals correlate with the tide
%        difference leaks into the 1-D estimate by v*cov(dt,dtide)/var(dtide);
%        the trend-controlled scan (opts.dt) recovers a and the rate v
%     9. vdef.trackBase follows a base that ramps out of any fixed window
%    10. a response that LAGS a diurnal tide is recovered by the quadrature
%        fit (opts.dq) as the right b/a, i.e. the right lag; an in-phase
%        response gives b = 0 and the same a as the grid scan
%
%   Runs in MATLAB or Octave:
%     docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work/scripts \
%       gnuoctave/octave:latest octave --no-gui test_tidal_stack.m

addpath(fileparts(fileparts(mfilename('fullpath'))));   % +vdef
rand('seed', 9); randn('seed', 9);   %#ok<RAND>

nz = 5; nb = 4; np = 60;
dtide = -1 + 2*rand(1, np);                      % m, pair tide differences
kz = -54.4 * ones(nz, 1);                        % rad per metre of dh, ~100 m in ice
a_true = 1e-3 * (-4 + 8*rand(nz, nb));           % m per m of tide
sig_ph = 0.02;                                    % rad of phase noise per block-mean field
% expected 1-sigma of the scan: sig_ph / (|k| sqrt(np) std(dtide)) ~ 0.08 mm/m
I = zeros(nz, nb, np); W = ones(nz, nb, np);
for p = 1:np
  I(:,:,p) = exp(1i*(kz .* a_true * dtide(p) + sig_ph*randn(nz, nb)));
end
o = struct('nboot', 60, 'self_test', 1.2e-3);

%% 1. Recovery
R = vdef.tidalStack(I, W, dtide, kz, o);
err = R.a - a_true; step = mean(diff(R.agrid));
fprintf('1. recovery: max |err| %.3f mm/m (grid step %.2f), min F %.2f\n', ...
  1e3*max(abs(err(:))), 1e3*step, min(R.F(:)));
assert(max(abs(err(:))) < 0.35e-3, 'a(z,x) not recovered: max error %.3f mm/m', 1e3*max(abs(err(:))));
assert(min(R.F(:)) > 0.8, 'coherent peak too low (%.2f) for 0.3 rad phase noise', min(R.F(:)));

%% 2. Self-test is exact
fprintf('2. self-test: injected 1.20 mm/m, worst shift error %.2e mm/m, %d edge cells\n', 1e3*R.inj_err, R.n_edge);
assert(R.inj_err < 1e-9, 'injected response did not return as an exact shift (%.3g)', R.inj_err);

%% 3. Bootstrap sd covers the error
z = abs(err) ./ max(R.a_sd, step/2);
fprintf('3. bootstrap: rms |err|/sd %.2f, max %.2f\n', sqrt(mean(z(:).^2)), max(z(:)));
assert(sqrt(mean(z(:).^2)) < 2.0, 'bootstrap sd understates the error (rms z %.2f)', sqrt(mean(z(:).^2)));

%% 4. Separable artefact
g_true = 1e-3 * (-3 + 6*rand(nb, 1));            % m/m per ns, per block
d_ind = 0.4*randn(nb, np);                       % ns, uncorrelated with the tide
Ia = I;
for p = 1:np, Ia(:,:,p) = I(:,:,p) .* exp(1i * kz * (g_true .* d_ind(:,p)).'); end
o4 = o; o4.dres = d_ind;
R4 = vdef.tidalStack(Ia, W, dtide, kz, o4);
ea = R4.a2 - a_true; eg = R4.g2 - repmat(g_true.', nz, 1);
fprintf('4. separable artefact: max |rcol| %.2f; a2 max err %.3f mm/m; g2 max err %.3f mm/m/ns\n', ...
  max(abs(R4.rcol)), 1e3*max(abs(ea(:))), 1e3*max(abs(eg(:))));
assert(max(abs(R4.rcol)) < 0.5, 'test setup: d_res should be uncorrelated with the tide');
assert(max(abs(ea(:))) < 0.5e-3, 'a not separated from an independent artefact (%.3f mm/m)', 1e3*max(abs(ea(:))));
assert(max(abs(eg(:))) < 0.6e-3, 'artefact gain not recovered (%.3f)', 1e3*max(abs(eg(:))));
% and the 1-D estimate is now biased by the artefact, as it must be
e1 = R4.a - a_true;
fprintf('   (1-D estimate, artefact uncontrolled: max err %.3f mm/m)\n', 1e3*max(abs(e1(:))));

%% 5. Collinear artefact is REPORTED, not absorbed
d_col = repmat(dtide, nb, 1) * 5 + 0.15*randn(nb, np);   % ~0.96 correlated with the tide
Ic = I;
for p = 1:np, Ic(:,:,p) = I(:,:,p) .* exp(1i * kz * (g_true .* d_col(:,p)).'); end
o5 = o; o5.dres = d_col;
R5 = vdef.tidalStack(Ic, W, dtide, kz, o5);
fprintf('5. collinear artefact: min |rcol| %.2f (must be reported high); a2 max err %.2f mm/m\n', ...
  min(abs(R5.rcol)), 1e3*max(abs(R5.a2(:) - a_true(:))));
assert(min(abs(R5.rcol)) > 0.9, 'collinearity not reported (min |rcol| %.2f)', min(abs(R5.rcol)));

%% 6. complexBlocks
Nt = 60; Nx = 40; rows = [10 20 30];
ig = exp(1i*0.7) * ones(Nt, Nx); coh = 0.9*ones(Nt, Nx);
ref = 5*ones(1, Nx); ref(3) = NaN;               % column 3 has no reference
coh(:, 7) = 0.1;                                  % column 7 below threshold
[Ib, Wb] = vdef.complexBlocks(ig, coh, ref, [1 21], 20, 0.3, rows);
fprintf('6. complexBlocks: |I| %.3f, angle %.2e, block-1 weight %.1f (18 usable cols x 0.9 = %.1f)\n', ...
  abs(Ib(1,1)), angle(Ib(1,1)), Wb(1,1), 18*0.9);
assert(all(abs(abs(Ib(:)) - 1) < 1e-12) && all(abs(angle(Ib(:))) < 1e-12), ...
  'phase referencing should give unit modulus and zero phase for a constant field');
assert(abs(Wb(1,1) - 18*0.9) < 1e-9, 'unreferenced and low-coherence columns must carry zero weight');

%% 7. Pass-level errors: the jackknife sees them, the pair bootstrap does not
% Each pass gets an error of its own that enters every pair it is in -
% phase_p += delta_j - delta_i - as a registration or surface error does.
% Over independent draws of delta the empirical scatter of a is the truth.
npass = 15; [pj, pi_] = meshgrid(1:npass, 1:npass); sel = pi_ < pj;
pairs7 = [pi_(sel) pj(sel)]; np7 = size(pairs7, 1);          % 105 pairs
Tp = -0.5 + rand(1, npass);                                  % pass tides
dtide7 = Tp(pairs7(:,2)) - Tp(pairs7(:,1));
sig_pass = 0.05;                                              % rad per pass
nrep = 12; est = nan(nz, nb, nrep); sd_jk = nan(nz, nb, nrep); sd_bt = nan(nz, nb, nrep);
o7 = struct('nboot', 40, 'self_test', 0, 'pairs', pairs7);
for r = 1:nrep
  delta = sig_pass * randn(1, npass);
  I7 = zeros(nz, nb, np7); W7 = ones(nz, nb, np7);
  for p = 1:np7
    I7(:,:,p) = exp(1i*(kz .* a_true * dtide7(p) + (delta(pairs7(p,2)) - delta(pairs7(p,1))) ...
      + 0.01*randn(nz, nb)));
  end
  R7 = vdef.tidalStack(I7, W7, dtide7, kz, o7);
  est(:,:,r) = R7.a; sd_jk(:,:,r) = R7.a_sd_jk; sd_bt(:,:,r) = R7.a_sd_boot;
end
emp = std(est, 0, 3); jk = mean(sd_jk, 3); bt = mean(sd_bt, 3);
q_jk = median(jk(:) ./ emp(:)); q_bt = median(bt(:) ./ emp(:));
fprintf('7. pass-level errors: empirical sd %.3f mm/m; jackknife %.3f (ratio %.2f), pair bootstrap %.3f (ratio %.2f)\n', ...
  1e3*median(emp(:)), 1e3*median(jk(:)), q_jk, 1e3*median(bt(:)), q_bt);
assert(q_jk > 0.6 && q_jk < 1.7, 'jackknife sd does not track the empirical scatter (ratio %.2f)', q_jk);
assert(q_bt < 0.7*q_jk, 'the pair bootstrap should understate pass-level errors (ratios %.2f vs %.2f)', q_bt, q_jk);
assert(isequal(R7.a_sd, R7.a_sd_jk), 'with pairs given, a_sd must be the jackknife');
% the trend-controlled estimate gets a jackknife too
dt7 = 0.5 + 3*rand(1, np7);
R7t = vdef.tidalStack(I7, W7, dtide7, kz, setfield(o7, 'dt', dt7)); %#ok<SFLD>
assert(isequal(size(R7t.a_t_sd), [nz nb]) && all(isfinite(R7t.a_t_sd(:))), 'a_t_sd missing or not finite');

%% 8. A tide-correlated trend leaks into a; opts.dt removes it
v_true = 1.0e-3;                                              % m per day, firn-compaction scale
dt8 = abs(1.5 + 0.8*dtide + 0.6*randn(1, np));                % days, r ~ 0.6 with dtide, as on the lines
c8 = corrcoef(dt8, dtide); r8 = c8(1,2);
cxy = mean((dt8 - mean(dt8)) .* (dtide - mean(dtide))); vxx = mean((dtide - mean(dtide)).^2);
bias_pred = v_true * cxy / vxx;
I8 = I;
for p = 1:np, I8(:,:,p) = I(:,:,p) .* exp(1i * kz * v_true * dt8(p)); end
o8 = o; o8.dt = dt8;
R8 = vdef.tidalStack(I8, W, dtide, kz, o8);
bias_1d = median(R8.a(:) - a_true(:));
e_t = R8.a_t - a_true; e_v = R8.v_t - v_true;
fprintf(['8. trend leak: corr(dt, dtide) %.2f (reported %.2f); 1-D bias %.2f mm/m against %.2f predicted; ' ...
         'trend-controlled: a max err %.3f mm/m, v median err %.3f mm/day\n'], ...
  r8, R8.r_tt, 1e3*bias_1d, 1e3*bias_pred, 1e3*max(abs(e_t(:))), 1e3*median(abs(e_v(:))));
assert(abs(r8) > 0.5 && abs(bias_pred) > 0.5e-3, 'test setup: the trend must be correlated enough to leak');
assert(abs(R8.r_tt - r8) < 1e-9, 'r_tt not reported');
assert(abs(bias_1d - bias_pred) < 0.35*abs(bias_pred), '1-D bias %.2f mm/m, predicted %.2f', 1e3*bias_1d, 1e3*bias_pred);
assert(max(abs(e_t(:))) < 0.4e-3 && max(abs(e_t(:))) < 0.5*abs(bias_pred), ...
  'trend-controlled a not recovered (%.3f mm/m against a %.2f mm/m leak)', 1e3*max(abs(e_t(:))), 1e3*bias_pred);
assert(median(abs(e_v(:))) < 0.3e-3, 'rate not recovered (%.3f mm/day)', 1e3*median(abs(e_v(:))));

%% 9. The base tracker follows a base that leaves any fixed search window
% A flat base at 285 m that ramps up to 210 m over the last fifth of the
% line, under a brighter flat stripe at 280 m (the system artefact a
% windowed pick snapped onto on the real lines).
zt = (0:1:400).'; nct = 200; truth = 285*ones(1, nct); rr = 161:nct; truth(rr) = linspace(285, 210, numel(rr));
imgt = -120 + 2*randn(numel(zt), nct);
for c = 1:nct, imgt(abs(zt - truth(c)) <= 1, c) = -80; end
imgt(abs(zt - 280) <= 0.5, rr) = -75;
bt = vdef.trackBase(imgt, zt);
% Where the base passes within one search step (6 m) of the brighter stripe
% the tracker may sit on the stripe briefly; it must recover once they part.
near = abs(truth - 280) <= 6 & (1:nct) >= rr(1); err = abs(bt - truth);
fprintf('9. base tracker: max |error| %.1f m away from the stripe, %.1f m where the base crosses it (ramp 285 -> 210 m)\n', ...
  max(err(~near)), max(err(near)));
assert(max(err(~near)) <= 3, 'base not tracked (max error %.1f m away from the stripe)', max(err(~near)));
assert(max(err(near)) <= 6, 'tracker strayed more than one step at the stripe crossing (%.1f m)', max(err(near)));

%% 10. Quadrature tide: a lagging response is seen, an elastic one is not
% 13 passes at random times over 3 days under a pure diurnal tide, every
% pair stacked, with a trend; the column responds to eta(t - tau).
Tq = 24.84/24; wq = 2*pi/Tq;                                 % days, rad per day
tp = sort(3*rand(1, 13)); [pj, pi_] = meshgrid(1:13, 1:13); sel = pi_ < pj;
pr = [pi_(sel) pj(sel)]; np10 = size(pr, 1);
eta = @(t) 0.55*sin(wq*t); etq = @(t) 0.55*cos(wq*t);         % etq = (T/2pi) deta/dt
dT = eta(tp(pr(:,2))) - eta(tp(pr(:,1))); dQ = etq(tp(pr(:,2))) - etq(tp(pr(:,1)));
dt10 = tp(pr(:,2)) - tp(pr(:,1)); v10 = 0.8e-3;
c10 = corrcoef(dT, dQ);
o10 = struct('nboot', 2, 'self_test', 0, 'pairs', pr, 'dt', dt10, 'dq', dQ);
for tau_h = [0 3]
  tau = tau_h/24; I10 = zeros(nz, nb, np10); W10 = ones(nz, nb, np10);
  resp = eta(tp(pr(:,2)) - tau) - eta(tp(pr(:,1)) - tau);    % the lagged tide difference
  for p = 1:np10
    I10(:,:,p) = exp(1i*(kz .* a_true * resp(p) + kz * v10 * dt10(p) + 0.01*randn(nz, nb)));
  end
  R10 = vdef.tidalStack(I10, W10, dT, kz, o10);
  b_true = -a_true * sin(wq*tau); a_amp = a_true * cos(wq*tau);
  e_a = R10.a_q - a_amp; e_b = R10.b_q - b_true; e_v = R10.v_q - v10;
  big = abs(a_true) > 1e-3;                                   % lag is defined where there is a response
  lag = atan2(-R10.b_q(big) .* sign(a_true(big)), R10.a_q(big) .* sign(a_true(big))) / wq * 24;
  fprintf(['10. quadrature, lag %d h: corr(dq, dtide) %.2f; max err a %.3f, b %.3f mm/m, v %.3f mm/day; ' ...
           'lag recovered %.2f to %.2f h; F %.2f; a_q - a_t (grid) max %.3f mm/m\n'], tau_h, c10(1,2), ...
    1e3*max(abs(e_a(:))), 1e3*max(abs(e_b(:))), 1e3*max(abs(e_v(:))), min(lag), max(lag), min(R10.F_q(:)), ...
    1e3*max(abs(R10.a_q(:) - R10.a_t(:))));
  % bounds from the realised design: least-squares sigma of each parameter
  % for 0.01 rad of phase noise, with the common phase free
  Xc = [dT; dQ; dt10] - mean([dT; dQ; dt10], 2); sg = 0.01/abs(kz(1)) * sqrt(diag(inv(Xc*Xc.')));
  assert(max(abs(e_a(:))) < 4.5*sg(1) && max(abs(e_b(:))) < 4.5*sg(2) && max(abs(e_v(:))) < 4.5*sg(3), ...
    'quadrature fit did not recover a, b, v (lag %d h; errors %.1f, %.1f, %.1f sigma)', tau_h, ...
    max(abs(e_a(:)))/sg(1), max(abs(e_b(:)))/sg(2), max(abs(e_v(:)))/sg(3));
  lag_sd = max(sg(1:2)) ./ abs(a_true(big)) / wq * 24;       % hours, small-angle
  assert(all(abs(lag - tau_h) < 4.5*lag_sd), 'lag not recovered (%d h)', tau_h);
  assert(all(isfinite(R10.a_q_sd(:))) && all(isfinite(R10.b_q_sd(:))), 'quadrature jackknife missing');
  if tau_h == 0
    assert(max(abs(R10.a_q(:) - R10.a_t(:))) < 1.5*step, 'in phase: the fit must stay at the grid peak');
  end
end

fprintf('\nall tidal-stack checks passed\n');
