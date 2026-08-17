%TEST_INVERT_NETWORK Unit test for vdef.invertNetwork.
%
%   What the network inversion has to make good: recover per-epoch
%   displacement from a redundant set of pair differences, without a
%   privileged reference, while identifying the pairs that do not fit.
%
%   Covers:
%     1. exact recovery on noiseless data, up to the additive constant the
%        data cannot determine
%     2. redundancy actually buys precision - the full network beats a
%        single-reference star network on the same noisy data
%     3. a corrupted pair is rejected, and its rejection restores the
%        answer rather than merely flagging it
%     4. the residuals ARE closure errors: on a noiseless network they
%        vanish, and a broken triangle shows up in them
%     5. an epoch left unreachable by rejection comes back NaN, not pinned
%        to the datum
%     6. many rows at once (a depth profile) invert independently
%
%   Runs in MATLAB or Octave:
%     docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work/scripts \
%       gnuoctave/octave:latest octave --no-gui test_invert_network.m

addpath(fileparts(fileparts(mfilename('fullpath'))));   % +vdef
rand('seed', 17); randn('seed', 17);   %#ok<RAND>

Nep = 10;
x_true = cumsum([0 0.6 -0.2 0.9 0.4 -0.5 1.1 0.3 -0.7 0.8]);   % arbitrary units
x_true = x_true - mean(x_true);        % the datum invertNetwork adopts

% every unordered pair
pairs = [];
for i = 1:Nep
  for j = i+1:Nep
    pairs(end+1,:) = [i j]; %#ok<AGROW>
  end
end
Npair = size(pairs,1);
d_true = x_true(pairs(:,2)) - x_true(pairs(:,1));
fprintf('%d epochs, %d pairs (a single-reference star would use %d)\n', ...
  Nep, Npair, Nep-1);

%% 1. Noiseless recovery
N = vdef.invertNetwork(pairs, d_true, []);
err = N.x - x_true;
assert(max(abs(err)) < 1e-9, 'noiseless recovery is off by %.3g', max(abs(err)));
assert(max(abs(N.resid)) < 1e-9, 'noiseless residuals should vanish (max %.3g)', ...
  max(abs(N.resid)));
assert(abs(sum(N.x)) < 1e-9, 'the sum(x) = 0 datum was not applied');
fprintf('1. noiseless: max |error| %.2e, max |residual| %.2e\n', ...
  max(abs(err)), max(abs(N.resid)));

%% 2. Redundancy buys precision
% Same noise level on every pair. The star network uses only the pairs
% touching epoch 1; the full network uses all of them.
sig = 0.05; NTRIAL = 200;
e_full = zeros(1,NTRIAL); e_star = zeros(1,NTRIAL);
star = pairs(pairs(:,1) == 1, :);
for t = 1:NTRIAL
  dn = d_true + sig*randn(1,Npair);
  Nf = vdef.invertNetwork(pairs, dn, struct('n_sigma', Inf));
  e_full(t) = sqrt(mean((Nf.x - x_true).^2, 'omitnan'));

  ds = dn(pairs(:,1) == 1);
  Ns = vdef.invertNetwork(star, ds, struct('n_sigma', Inf, 'min_pairs', 3));
  e_star(t) = sqrt(mean((Ns.x - x_true).^2, 'omitnan'));
end
gain = mean(e_star)/mean(e_full);
fprintf('2. rms error: full network %.4f, star %.4f -> %.2fx better\n', ...
  mean(e_full), mean(e_star), gain);
assert(gain > 1.5, ...
  'the full network should beat a star by a clear margin, got %.2fx', gain);

%% 3. A corrupted pair is rejected, and rejecting it restores the answer
dn = d_true + 0.02*randn(1,Npair);
bad_idx = 17;
dn(bad_idx) = dn(bad_idx) + 3.0;            % one badly wrong pair
Nnorej = vdef.invertNetwork(pairs, dn, struct('n_sigma', Inf));
Nrej   = vdef.invertNetwork(pairs, dn, struct('n_sigma', 3));
e_norej = sqrt(mean((Nnorej.x - x_true).^2,'omitnan'));
e_rej   = sqrt(mean((Nrej.x   - x_true).^2,'omitnan'));
fprintf('3. one corrupted pair: rms %.4f unrejected -> %.4f rejected\n', e_norej, e_rej);
assert(~Nrej.used(bad_idx), 'the corrupted pair was not rejected');
assert(e_rej < 0.5*e_norej, ...
  'rejection should restore the answer: %.4f vs %.4f', e_rej, e_norej);

%% 4. Residuals are closure errors
% Break ONE triangle consistently and check the residual pattern points at
% the pairs of that triangle rather than smearing over the whole network.
d4 = d_true;
tri = find((pairs(:,1)==2 & pairs(:,2)==5));
d4(tri) = d4(tri) + 1.0;
N4 = vdef.invertNetwork(pairs, d4, struct('n_sigma', Inf));
[~, worst] = max(abs(N4.resid));
fprintf('4. broken pair %d -> largest residual at pair %d\n', tri, worst);
assert(worst == tri, ...
  'the broken pair should carry the largest residual (got pair %d)', worst);

%% 5. Unreachable epochs come back NaN
% Pairs among epochs 1-5, plus one isolated pair (9,10) so the epoch count
% is still 10 and epochs 6-8 are touched by nothing at all. Only the
% LARGEST component shares a datum with the constraint, so 9 and 10 must
% come back NaN too even though they are measured against each other.
sub = [pairs(all(pairs <= 5, 2), :); 9 10];
d5  = x_true(sub(:,2)) - x_true(sub(:,1));
N5  = vdef.invertNetwork(sub, d5, []);
assert(all(isfinite(N5.x(1:5))), 'the connected epochs should be solved');
assert(all(isnan(N5.x(6:10))), ...
  'epochs no pair touches must be NaN, not pinned to the datum');
fprintf('5. disconnected: %d of %d epochs solved, rest NaN\n', ...
  nnz(isfinite(N5.x)), Nep);

%% 6. Many rows invert independently
Nz = 40;
prof = (1:Nz).'/Nz;                       % depth-dependent scaling
D = prof * d_true;                        % Nz x Npair
Nm = vdef.invertNetwork(pairs, D, []);
X_expect = prof * x_true;
assert(max(abs(Nm.x(:) - X_expect(:))) < 1e-9, 'multi-row inversion is wrong');
assert(size(Nm.x,1) == Nz, 'expected %d rows out', Nz);
fprintf('6. %d rows inverted independently, max |error| %.2e\n', ...
  Nz, max(abs(Nm.x(:) - X_expect(:))));

fprintf('\ntest_invert_network: all checks passed\n');
