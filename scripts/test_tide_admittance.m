%TEST_TIDE_ADMITTANCE Unit test for vdef.fitTideAdmittance.
%   The claim that function exists to make good: the tide admittance of a
%   joint strain = a + b*t + c*tide fit is invariant to the choice of
%   reference pass, while the plain correlation of strain with tide is
%   not, because the secular trend aliases into it through the sample
%   covariance of t with tide. That aliasing is what left the two builds
%   of leg 1 (references three days apart) with different absolute r
%   profiles after the coalignment fix, and it is why the acceptance test
%   moved from r to the admittance.
%
%   Covers, on synthetic series shaped like the EAGER 2022 line (13
%   pairs over 3.2 days, diurnal tide of ~1 m range, an admittance
%   profile that changes sign along track like the flexure hinge):
%     1. recovery of trend and admittance from noisy series, with the
%        quoted 1-sigma consistent with the actual errors
%     2. exact reference invariance: re-referencing strain, t and tide to
%        another pass changes intercept only - trend, admittance and
%        r_partial reproduce to machine precision
%     3. the aliasing the fit removes: two noiseless builds of the same
%        truth referencing different epochs give bit-equal admittance but
%        visibly different plain correlations (and the exact fit drives
%        r_partial to +/-1, exercising that branch)
%     4. NaN pairs are dropped per block; blocks below min_pairs and
%        all-NaN blocks come back NaN with n reported
%     5. a tide collinear with time is refused, not fitted: r_tt is
%        reported and the coefficients stay NaN
%
%   Runs in MATLAB or Octave:
%     docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work/scripts \
%       gnuoctave/octave:latest octave --no-gui test_tide_admittance.m

addpath(fileparts(fileparts(mfilename('fullpath'))));   % +vdef

rand('seed', 7); randn('seed', 7);   %#ok<RAND> % Octave-compatible seeding

%% Common geometry: the real line's sampling, in miniature
% =====================================================================
Npair = 13;
Nblk  = 24;
% Pass times [days] over the 12-09..12-12 window - FIXED, not drawn, so
% that test 3's aliasing magnitude is identical in Octave and MATLAB
% (their legacy seeded generators produce different streams). The clumps
% mimic the real out-and-back walking pattern.
t     = [0 0.02 0.45 0.48 1.02 1.05 1.55 2.02 2.05 2.55 3.05 3.15 3.20];
M2    = 12.42/24;                              % [days] principal lunar semidiurnal
tide  = 0.55 * sin(2*pi*t/M2 + 0.8) ...
      + 0.25 * sin(2*pi*t/1.0758 + 2.1);       % [m] two-constituent tide, ~1.4 m range

a_true = 1e-4 * randn(1, Nblk);                          % reference constants
b_true = -80e-6 + 10e-6 * randn(1, Nblk);                % [strain/day] secular trend
c_true = 120e-6 * cos(pi * (0:Nblk-1) / (Nblk-1));       % [strain/m] hinge: + to -

noise_sigma = 5e-6;

truth = @(a, b, c) a(:) * ones(1,Npair) + b(:) * t + c(:) * tide;
strain = truth(a_true, b_true, c_true) + noise_sigma * randn(Nblk, Npair);

%% 1. Recovery from noisy series
% =====================================================================
A = vdef.fitTideAdmittance(strain, t, tide);

assert(all(A.n == Npair), 'expected every pair used in every block');
err_c = A.admittance - c_true;
err_b = A.trend - b_true;
assert(all(abs(err_c) < 4 * A.admittance_std), ...
  'admittance error exceeds 4 sigma (max %.3g at %.3g sigma)', ...
  max(abs(err_c)), max(abs(err_c ./ A.admittance_std)));
assert(all(abs(err_b) < 4 * A.trend_std), 'trend error exceeds 4 sigma');
% The quoted sigma should describe the actual scatter: the rms error over
% 24 blocks should be within a factor 2 of the mean quoted sigma
rms_c = sqrt(mean(err_c.^2));
assert(rms_c > 0.5*mean(A.admittance_std) && rms_c < 2*mean(A.admittance_std), ...
  'quoted admittance sigma (%.3g) inconsistent with actual rms error (%.3g)', ...
  mean(A.admittance_std), rms_c);
% The hinge survives: r_partial changes sign along the profile
assert(A.r_partial(1) > 0.5 && A.r_partial(end) < -0.5, ...
  'r_partial does not resolve the sign change of the admittance profile');
fprintf('1. recovery: rms admittance error %.2f ue/m against truth range %.0f ue/m\n', ...
  1e6*rms_c, 1e6*(max(c_true)-min(c_true)));

%% 2. Exact reference invariance under a common re-referencing
% =====================================================================
k0 = 9;   % adopt pair 9's epoch as the new reference
A2 = vdef.fitTideAdmittance(strain - strain(:,k0)*ones(1,Npair), ...
  t - t(k0), tide - tide(k0));

tol = 1e-9;   % relative, against the truth scale
assert(max(abs(A2.admittance - A.admittance)) < tol * max(abs(c_true)), ...
  're-referencing changed the admittance');
assert(max(abs(A2.trend - A.trend)) < tol * max(abs(b_true)), ...
  're-referencing changed the trend');
assert(max(abs(A2.r_partial - A.r_partial)) < 1e-9, ...
  're-referencing changed r_partial');
fprintf('2. invariance: re-referenced admittance agrees to %.2e ue/m\n', ...
  1e6*max(abs(A2.admittance - A.admittance)));

%% 3. What the fit removes: reference-epoch aliasing of the plain r
% =====================================================================
% Two noiseless builds of the SAME ice: build A references the first
% pass, build B the last, and each build measures only its own pairs
% (the reference pass itself contributes no pair). Same truth, no noise -
% any difference between the builds is pure analysis artefact.
iA = 2:Npair;                    % pairs of build A (ref = pass 1)
iB = 1:(Npair-1);                % pairs of build B (ref = pass Npair)
sA = truth(a_true, b_true, c_true);   % already relative to an arbitrary origin
sB = sA - sA(:,Npair)*ones(1,Npair);  % the same series re-referenced to pass Npair

FA = vdef.fitTideAdmittance(sA(:,iA), t(iA), tide(iA));
FB = vdef.fitTideAdmittance(sB(:,iB), t(iB), tide(iB));

assert(max(abs(FA.admittance - c_true)) < 1e-12, 'noiseless build A not exact');
assert(max(abs(FB.admittance - c_true)) < 1e-12, 'noiseless build B not exact');
assert(all(abs(FA.r_partial) > 1 - 1e-6) && all(abs(FB.r_partial) > 1 - 1e-6), ...
  'exact fit should saturate r_partial at +/-1');

rA = nan(1, Nblk); rB = nan(1, Nblk);
for bk = 1:Nblk
  rm = corrcoef(tide(iA), sA(bk,iA)); rA(bk) = rm(1,2);
  rm = corrcoef(tide(iB), sB(bk,iB)); rB(bk) = rm(1,2);
end
dr = max(abs(rA - rB));
assert(dr > 0.02, ...
  'expected the plain correlation to differ between the builds (got %.3g): the aliasing this test demonstrates has vanished', dr);
fprintf('3. aliasing: plain r differs by up to %.2f between builds; admittance exact in both\n', dr);

%% 4. NaN pairs and short blocks
% =====================================================================
s4 = strain;
s4(3, [1 4 8]) = NaN;            % a block with three pairs knocked out
s4(5, 1:Npair-3) = NaN;          % a block below min_pairs
s4(7, :) = NaN;                  % a dead block
A4 = vdef.fitTideAdmittance(s4, t, tide);

assert(A4.n(3) == Npair-3 && isfinite(A4.admittance(3)), ...
  'partially-NaN block should still fit on the surviving pairs');
assert(abs(A4.admittance(3) - c_true(3)) < 4 * A4.admittance_std(3), ...
  'partially-NaN block no longer recovers its admittance');
assert(A4.n(5) == 3 && ~isfinite(A4.admittance(5)), ...
  'block below min_pairs must come back NaN');
assert(A4.n(7) == 0 && ~isfinite(A4.admittance(7)) && ~isfinite(A4.r_partial(7)), ...
  'all-NaN block must come back NaN with n = 0');
untouched = setdiff(1:Nblk, [3 5 7]);
assert(isequal(A4.admittance(untouched), A.admittance(untouched)), ...
  'NaNs in one block must not perturb the others');
fprintf('4. NaN handling: %d/%d/%d pairs used in the damaged blocks\n', ...
  A4.n(3), A4.n(5), A4.n(7));

%% 5. Collinear tide is refused
% =====================================================================
tide5 = 0.3 * (t - mean(t));     % tide exactly linear in time
s5 = a_true(:)*ones(1,Npair) + b_true(:)*t + c_true(1)*ones(Nblk,1)*tide5;
A5 = vdef.fitTideAdmittance(s5, t, tide5);
assert(all(abs(A5.r_tt - 1) < 1e-12), 'r_tt should report the collinearity');
assert(all(~isfinite(A5.admittance)), ...
  'a tide collinear with time must be refused, not fitted');
fprintf('5. collinearity: refused with r_tt = %.6f\n', A5.r_tt(1));

fprintf('\ntest_tide_admittance: all checks passed\n');
