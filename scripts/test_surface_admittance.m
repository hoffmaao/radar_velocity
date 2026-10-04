%TEST_SURFACE_ADMITTANCE Unit test for vdef.surfaceAdmittance.
%   The claim the function exists to make good: the surface tidal
%   admittance a(x) regressed against an EXTERNAL tide is unbiased in the
%   presence of per-pass height offsets that are uniform along the line,
%   while the line-mean regressor the chain first used is compressed
%   toward one by them - and that a pass with a gross height error is
%   rejected rather than fitted.
%
%   Covers, on synthetic passes shaped like the EAGER 2022 line (13
%   passes, ~1 m tide range, 10 blocks of 200 samples, a(x) rising from
%   0.2 to 1.4 seaward):
%     1. clean data: a(x) recovered within its quoted sigma, and the
%        legacy line-mean regressor equals the direct slope formula
%     2. uniform per-pass offsets of 10 cm rms: the line-mean estimator
%        compresses the profile (slope of estimate on truth well below 1),
%        the external-tide estimator does not
%     3. the common-mode nuisance column leaves a(x) unchanged and
%        brings a_std down to the block's own noise - which, by design,
%        excludes the uniform shift of the whole profile that the
%        offsets leave behind
%     4. the pass gate rejects an injected 1 m pass and only that pass,
%        and a(x) then matches the fit that never saw it
%
%   Runs in MATLAB or Octave:
%     docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work/scripts \
%       gnuoctave/octave:latest octave --no-gui test_surface_admittance.m

addpath(fileparts(fileparts(mfilename('fullpath'))));   % +vdef
rand('seed', 5); randn('seed', 5);   %#ok<RAND>

Np = 13; BLOCK = 200; nb = 10; Nx = nb*BLOCK;
a_true = linspace(0.2, 1.4, nb).';
a_x    = kron(a_true, ones(BLOCK,1));          % per sample
tide   = -0.5 + rand(1, Np);                   % ~1 m range
sig_s  = 0.15;                                 % per-sample height noise [m]
sig_b  = sig_s / sqrt(BLOCK);                  % block-mean noise, ~1 cm
opts   = struct('block', BLOCK, 'min_pass', 5);

%% 1. Clean recovery, and the legacy estimator reproduces the direct slope
Z = a_x * tide + sig_s * randn(Nx, Np);
S = vdef.surfaceAdmittance(Z, tide, opts);
z1 = (S.a - a_true) ./ S.a_std;
fprintf('1. clean: a recovered, |z| max %.2f, rms %.2f (expect ~1), all %d passes ok\n', ...
  max(abs(z1)), sqrt(mean(z1.^2)), nnz(S.pass_ok));
assert(all(S.pass_ok), 'the gate rejected a clean pass');
assert(max(abs(z1)) < 3.5 && sqrt(mean(z1.^2)) < 1.8, 'clean a(x) not recovered within sigma');
assert(all(S.n == Np), 'not every pass was used');

Slm = vdef.surfaceAdmittance(Z, [], opts);
lm  = mean(Z, 1);
b   = 3; zb = mean(Z((b-1)*BLOCK+1:b*BLOCK, :), 1);
p   = polyfit(lm, zb, 1);
assert(abs(Slm.a(b) - p(1)) < 1e-12, 'legacy line-mean estimator differs from the direct slope');
assert(~Slm.common_mode && isnan(Slm.abar), 'legacy mode should carry no gate or nuisance column');
fprintf('   legacy line-mean estimator matches polyfit slope (%.4f)\n', p(1));

%% 2. Uniform per-pass offsets compress the line-mean estimator only
u  = 0.10 * randn(1, Np);
Zu = a_x * tide + repmat(u, Nx, 1) + sig_s * randn(Nx, Np);
Se = vdef.surfaceAdmittance(Zu, tide, opts);
Sl = vdef.surfaceAdmittance(Zu, [], opts);
% The external-tide estimate is affine in the truth: slope one, plus a
% uniform shift from the offsets' sample covariance with the tide (the
% part no regressor can remove). The line-mean estimate is compressed,
% which shows as a slope below one once both are put in line-mean units.
pe = polyfit(a_true, Se.a, 1);
al = Sl.a / mean(Sl.a); at = a_true / mean(a_true);
pl = polyfit(at, al, 1);
fprintf('2. 10 cm uniform offsets: external tide a = %.3f*truth %+.3f; line-mean shape slope %.3f\n', ...
  pe(1), pe(2), pl(1));
assert(abs(pe(1) - 1) < 0.05, 'external-tide a(x) is distorted (slope %.3f)', pe(1));
assert(pl(1) < 0.93, 'line-mean a(x) should be compressed here (slope %.3f)', pl(1));
assert(all(Se.pass_ok), 'the gate rejected a pass with ordinary offsets');
fprintf('   line-mean residual rms %.3f m against injected %.3f\n', ...
  sqrt(mean(Se.pass_resid.^2)), sqrt(mean(u.^2)));

%% 3. The nuisance column changes a_std, not a
Sn = vdef.surfaceAdmittance(Zu, tide, setfield(opts, 'common_mode', false)); %#ok<SFLD>
assert(max(abs(Sn.a - Se.a)) < 1e-10, 'the common-mode column moved a(x)');
ratio = median(Se.a_std ./ Sn.a_std);
fprintf('3. a_std with/without the common-mode column: median ratio %.2f; with it %.4f vs block noise/(sd(tide) sqrt(Np)) %.4f\n', ...
  ratio, median(Se.a_std), sig_b / std(tide) / sqrt(Np));
assert(ratio < 0.6, 'the nuisance column did not absorb the shared error (ratio %.2f)', ratio);
% What a_std then measures is the block's OWN scatter. The error it does
% not carry is the uniform shift of the whole profile (pe(2) above), which
% is shared by every block - so the test removes the common shift before
% asking whether the quoted sigmas describe what is left.
shift = mean(Se.a - a_true);
zn = (Se.a - a_true - shift) ./ Se.a_std;
fprintf('   common shift %+.3f (not in a_std, by design); after removing it rms z = %.2f\n', shift, sqrt(mean(zn.^2)));
assert(sqrt(mean(zn.^2)) < 2.5, 'a_std with the column understates the block scatter (rms z %.2f)', sqrt(mean(zn.^2)));

%% 4. The pass gate
Zg = Zu; kbad = 4; Zg(:, kbad) = Zg(:, kbad) + 1.1;
Sg = vdef.surfaceAdmittance(Zg, tide, opts);
fprintf('4. pass %d offset by 1.1 m: rejected = %s; residual there %+.2f m, others rms %.3f\n', ...
  kbad, mat2str(find(~Sg.pass_ok)), Sg.pass_resid(kbad), ...
  sqrt(mean(Sg.pass_resid(Sg.pass_ok).^2)));
assert(~Sg.pass_ok(kbad) && nnz(~Sg.pass_ok) == 1, 'the gate did not reject exactly the bad pass');
Zd = Zu; Zd(:, kbad) = NaN; td = tide; td(kbad) = NaN;
Sd = vdef.surfaceAdmittance(Zd, td, opts);
assert(max(abs(Sg.a - Sd.a)) < 1e-10, 'gated a(x) differs from the fit without the pass');
assert(all(Sg.n == Np - 1), 'block pass counts should drop by one');
Sng = vdef.surfaceAdmittance(Zg, tide, setfield(opts, 'max_pass_sigma', Inf)); %#ok<SFLD>
assert(all(Sng.pass_ok), 'max_pass_sigma = Inf should disable the gate');
fprintf('   without the gate the bad pass distorts a(x) by up to %.3f\n', max(abs(Sng.a - Sd.a)));

fprintf('\nall surface-admittance checks passed\n');
