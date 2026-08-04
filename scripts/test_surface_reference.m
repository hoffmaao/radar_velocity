%TEST_SURFACE_REFERENCE Regression test for the surface reference bin.
%   vdef.differentialRange references every trace to zero just below that
%   trace's own surface return. A trace whose Surface is NaN, or whose
%   reference twtt falls outside the fast-time axis, has NO reference bin,
%   and the only defensible output for it is "no data".
%
%   The failure this guards against is silent: two-argument min/max ignore
%   NaN, so clamping an unconditioned ref_bin with min(max(r,1),Nt) leaves
%   the NaN alone but a bounds violation snaps to bin 1 or bin Nt. The
%   trace then gets referenced to the top of the record instead of to its
%   surface, and comes back as a confident, fully finite, entirely wrong
%   dtau column that the downstream block average happily folds in.
%
%   Real products make this reachable rather than hypothetical: Surface
%   comes from pass.surface, which carries NaN wherever the surface
%   tracker failed, and the reference offset pushes traces near the end of
%   the record past the last bin.
%
%   Covers, on one synthetic pair with a known dtau:
%     1. NaN Surface            -> no reference bin, whole trace invalid
%     2. Surface past the axis  -> same, not clamped to the last bin
%     3. Surface before the axis-> same, not clamped to the first bin
%     4. good traces            -> unaffected, dtau still recovered
%     5. the already-unwrapped branch honours the same rule
%     6. vdef.blockAverage does not let an unreferenced trace into a mean
%
%   Runs in MATLAB or Octave:
%     docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work/scripts \
%       gnuoctave/octave:latest octave --no-gui test_surface_reference.m

addpath(fileparts(fileparts(mfilename('fullpath'))));   % +vdef

rand('seed', 3); randn('seed', 3);   %#ok<RAND> % Octave-compatible seeding

fc      = 750e6;
fs      = 300e6;
Nt      = 400;
Nx      = 60;
Surface = 0.20e-6;

Time = (0:Nt-1).' / fs;

opts = [];
opts.phase_sign          = -1;
opts.ref_twtt_offset     = 50e-9;
opts.coherence_threshold = 0.3;
opts.max_gap_bins        = 20;
opts.min_coverage        = 0.3;
opts.block_size          = Nx;
opts.mlook_window        = [5 15];

% A noiseless, perfectly coherent pair with dtau linear in fast time, so
% any surviving column can be checked against the truth exactly.
dtau_true = (30e-12) * (Time - Surface) / (Time(end) - Surface);
dtau_true(Time < Surface) = 0;

map = [];
map.Time      = Time;
map.fc        = fc;
map.phase     = angle(exp(-1i*2*pi*fc*repmat(dtau_true, 1, Nx)));
map.coherence = ones(Nt, Nx);
map.phase_is_unwrapped = false;

% Column layout: three columns that cannot be referenced, the rest good.
COL_NAN   = 5;                       % surface tracker failed here
COL_LATE  = 12;                      % reference offset pushes past the axis
COL_EARLY = 20;                      % surface sits before the record starts
bad_cols  = [COL_NAN COL_LATE COL_EARLY];
good_cols = setdiff(1:Nx, bad_cols);

map.Surface            = repmat(Surface, 1, Nx);
map.Surface(COL_NAN)   = NaN;
map.Surface(COL_LATE)  = Time(end);        % + offset lands beyond Time(end)
map.Surface(COL_EARLY) = Time(1) - 1e-6;   % + offset still before Time(1)

[dtau, info] = vdef.differentialRange(map, opts);

fprintf('Surface reference regression, %d traces (%d unreferencable)\n', ...
  Nx, numel(bad_cols));
fprintf('%-8s %-12s %-10s %-14s %-12s\n', ...
  'col', 'Surface', 'ref_bin', 'max_valid_bin', 'finite dtau');
for x = [bad_cols good_cols(1)]
  fprintf('%-8d %-12.4g %-10.4g %-14d %-12d\n', ...
    x, map.Surface(x), info.ref_bin(x), info.max_valid_bin(x), nnz(isfinite(dtau(:,x))));
end

%% 1-3. A trace with no reference bin must produce nothing at all
for x = bad_cols
  assert(~isfinite(info.ref_bin(x)), ...
    'col %d has no valid surface, so ref_bin must be NaN (got %g) - a clamped bin silently references the wrong depth', ...
    x, info.ref_bin(x));
  assert(info.max_valid_bin(x) == 0, ...
    'col %d must report max_valid_bin 0, got %d', x, info.max_valid_bin(x));
  assert(~any(isfinite(dtau(:,x))), ...
    'col %d has no surface reference but returned %d finite dtau samples', ...
    x, nnz(isfinite(dtau(:,x))));
  assert(~any(info.valid(:,x)), 'col %d must be entirely invalid', x);
end

%% 4. The good traces are untouched, and still recover the truth
ref_bin_good = info.ref_bin(good_cols(1));
assert(all(info.ref_bin(good_cols) == ref_bin_good), ...
  'the referencable traces must all share the same reference bin here');
assert(all(info.max_valid_bin(good_cols) == Nt), ...
  'a fully coherent referencable trace must stay valid to the last bin');

expect = dtau_true - dtau_true(ref_bin_good);
cmp    = ref_bin_good:Nt;
err    = dtau(cmp, good_cols) - repmat(expect(cmp), 1, numel(good_cols));
assert(max(abs(err(:))) < 1e-15, ...
  'the referencable traces must still recover dtau (max error %.3g s)', max(abs(err(:))));
fprintf('\ngood traces recover dtau to %.2e s; %d of %d traces dropped\n', ...
  max(abs(err(:))), Nx - numel(good_cols), Nx);

%% 5. The already-unwrapped branch obeys the same rule
map_uw = map;
map_uw.phase_is_unwrapped = true;
map_uw.phase = -2*pi*fc*repmat(dtau_true, 1, Nx);
[dtau_uw, info_uw] = vdef.differentialRange(map_uw, opts);
for x = bad_cols
  assert(~any(isfinite(dtau_uw(:,x))), ...
    'col %d leaked through the already-unwrapped branch', x);
  assert(~isfinite(info_uw.ref_bin(x)), ...
    'col %d must have no ref_bin in the already-unwrapped branch', x);
end
assert(any(isfinite(dtau_uw(:, good_cols(1)))), ...
  'the already-unwrapped branch dropped a referencable trace too');

%% 6. An unreferencable trace must not reach a block mean
blk = vdef.blockAverage(dtau, map, info, opts);
assert(all(blk.coverage(:) <= (numel(good_cols)/Nx) + 1e-12), ...
  'block coverage %.4f exceeds the fraction of referencable traces %.4f - an unreferenced trace was averaged in', ...
  max(blk.coverage(:)), numel(good_cols)/Nx);

deep = isfinite(blk.dtau);
assert(any(deep(:)), 'the block average produced nothing at all');
blk_err = blk.dtau(cmp) - expect(cmp);
assert(max(abs(blk_err(isfinite(blk_err)))) < 1e-15, ...
  'the block mean was contaminated (max error %.3g s)', ...
  max(abs(blk_err(isfinite(blk_err)))));
fprintf('block coverage %.4f = %d/%d referencable traces, mean clean to %.2e s\n', ...
  max(blk.coverage(:)), numel(good_cols), Nx, ...
  max(abs(blk_err(isfinite(blk_err)))));

fprintf('\nPASS: traces without a surface reference are dropped, not clamped.\n');
