function OUT = elastic_modulus(opts)
%ELASTIC_MODULUS Effective Young's modulus of the grounding zone, from GPS flexure.
%
%   WHAT THIS DOES. Fits an elastic beam of varying thickness on a
%   hydrostatic foundation to the observed along-track profile of surface
%   tidal admittance, and reports the effective Young's modulus E* that
%   minimises the misfit. The method is that of Elgart, Minchew and Meyer
%   (2025), who apply it to ICESat-2 repeat-track flexure on the Ross Ice
%   Shelf with ROSETTA-Ice radar thickness; the substitutions here are the
%   survey's own GPS for the altimetry and this survey's own tracked bed
%   for the airborne radar thickness. Numerics live in vdef.beamFlexure,
%   vdef.surfaceAdmittance and vdef.invertElasticModulus.
%
%   WHY THE GPS AND NOT THE RADAR. The observable a beam model wants is
%   vertical DEFLECTION of the surface, and this survey measures that
%   directly and radar-independently: each pass's ref_z(x) is the platform
%   height on a floating shelf, so regressing it on the tide gives a(x),
%   the local admittance of the surface to the tide. That profile is the
%   one result of this project that never depended on the interferometry -
%   a(x) falls monotonically along the line by a factor of about six, and
%   the legs agree on it despite different main passes and pass sets. The
%   radar column-strain admittance is a different quantity - it goes as the
%   CURVATURE of this profile, not the profile - and it sits at the
%   measurement's systematic floor. So the deflection is what gets inverted
%   here.
%
%   THE TIDE IS CATS2008, NOT THE LINE MEAN. The first version of this
%   driver regressed each pass's heights on that pass's own line-mean
%   height. Measured against CATS2008 the line means carry a 13-15 cm rms
%   non-tidal residual, about half of it a height offset UNIFORM along the
%   line, and a uniform per-pass offset in the regressor compresses a(x)
%   toward a flat profile - which the beam fit reads as a stiffer beam
%   (+58% in E* for 14 cm rms on synthetic data with this geometry). So
%   the regressor is the external prediction (scripts/diagnostics/
%   cats2008_tide.m), the line mean is demoted to a pass gate and a
%   common-mode nuisance term (vdef.surfaceAdmittance), and the same tide
%   and gate are handed to the strain chain. What the external tide does
%   NOT remove is the offsets' random projection onto it, which shifts the
%   whole profile by a constant per line; that is what the pass JACKKNIFE
%   below measures, and it is the error to quote. opts.tide = 'linemean'
%   runs the legacy estimator for comparison.
%
%   THE NORMALISATION IS NOT A PROBLEM. The inversion eliminates the
%   amplitude in closed form and fits the SHAPE of the profile, so an
%   arbitrary rescaling leaves E* bit-identical (asserted in
%   scripts/test_beam_flexure.m, check 5). With the external tide a(x) is
%   in absolute metres per metre, and the fitted amplitude then reads
%   directly as the far-field admittance the beam implies.
%
%   READ THIS BEFORE QUOTING A NUMBER. This survey is a PARTIAL WINDOW. The
%   line spans about 4.8 km inside the flexure zone and reaches neither the
%   flat grounded end nor the flat floating one, so the clamp position sits
%   outside the data and trades off against E*. Check 7 of the unit test is
%   that exact geometry on synthetic data with a known answer: the modulus
%   comes back tens of percent off, with an interval wide enough to cover
%   the truth but seven times wider than the same method achieves on a full
%   profile. The formal intervals carry that ill-posedness. What they do
%   not carry is the error shared by every block of a line - a per-pass
%   height offset moves the whole profile and the amplitude absorbs it -
%   which is why the jackknife over passes is printed beside them, and why
%   the legs agreeing with each other is not a test of it: the legs are
%   the same passes, walked minutes apart, and share those offsets.
%
%   AND E* IS ONLY EVER E*h^3. Flexure constrains the rigidity, so the
%   modulus is no better known than the thickness cubed. Several thickness
%   cases are run rather than one, and the spread between them is part of
%   the answer rather than a sensitivity test bolted on afterwards.
%
%   A FUNCTION, NOT A SCRIPT, unlike the other diagnostics here. Two
%   reasons: local functions in a script are only defined once execution
%   reaches them in Octave, so a script version cannot be run outside
%   MATLAB; and a function cannot leak workspace state between invocations,
%   which is a trap this project has already been bitten by.
%
%   IT PRINTS AND RETURNS; IT DOES NOT DRAW. scripts/figures/flexure_inversion.m
%   calls this and does the drawing, so there is one implementation of the
%   inversion and one of the plot, and neither can drift from the other.
%
%   opts, all optional:
%     .mp_dir     multipass product directory
%     .products   product names to invert
%     .block      traces per along-track block (200 = 500 m)
%     .net_dir    all-pairs vvel products; source of the englacial strain
%                 admittance that joins the primary fit. The strain is the
%                 beam's own bending, dh = amp*K2*w''(x), sharing the
%                 shape fit's amplitude. Missing dir = all cases run
%                 shape-only.
%     .strain_source  'network' (default): the strain admittance at 100 m
%                 from the vvel networks in .net_dir, surface-referenced.
%                 'stack': the coherent tidal stack (scripts/diagnostics/
%                 tidal_stack.m) on the surface-coupled rebuild, read from
%                 .stack_file in vdef.figureDir - the in-phase response of
%                 the phase-lag fit (a_q, with its pass jackknife a_q_sd)
%                 at .stack_depth (100 m). Its phase reference is
%                 .stack_ref m below the surface (5: the chain's 50 ns; 60
%                 for tidal_stack_nozc_ref60.mat), so the lever is
%                 integrated from there, and with .stack_offset (default
%                 true) a line-uniform constant joins the strain rows to
%                 absorb the shallow-firn response a 5 m reference carries.
%                 The jackknife below then varies a(x) only: the stack's
%                 own error is already a pass jackknife.
%     .tide       'cats' (default) or 'linemean' - the regressor for a(x)
%                 and for the strain admittance
%     .max_pass_sigma  pass-gate threshold in robust sigmas of the
%                 line-mean residual against the tide (default 4)
%     .jackknife  leave-one-pass-out on the primary case (default true;
%                 one inversion per pass per line)
%     .profiles   inject a(x) profiles directly and skip the file load,
%                 for testing the inversion path without the products
%
%   Returns OUT.lines (the a(x) profiles, with the per-pass heights, tide
%   and gate they came from), OUT.fits (one inversion per thickness case
%   per line), OUT.jackknife (per line, the primary case's leave-one-pass-
%   out moduli and the standard error they imply), so the figures reuse
%   the numbers rather than refitting them. OUT.block echoes the block
%   size the fits used. OUT.constrained(c,i) is the averaging gate the
%   result figures share: E* bounded from above AND the shape data
%   described (reduced chi-squared under 3) - computed here, once, so the
%   figures cannot drift apart on the criterion.
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "addpath('<code>/scripts/diagnostics'); elastic_modulus"

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));                       % +vdef
addpath(here);                                             % cats2008_tide

if ~isfield(opts,'mp_dir') || isempty(opts.mp_dir)
  opts.mp_dir = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
end
if ~isfield(opts,'products') || isempty(opts.products)
  % EAGER_2022 is excluded throughout this project: it carries neither a
  % coregistration time shift nor an equalization vector, its
  % combine_passes input was replaced after the product was built so it
  % cannot be regenerated, and it is the same leg as GL1 anyway.
  opts.products = {'EAGER_2022_GL1','EAGER_2022_GL2', ...
                   'EAGER_2022_GL3','EAGER_2022_GL4'};
end
if ~isfield(opts,'block') || isempty(opts.block), opts.block = 200; end
if ~isfield(opts,'profiles'), opts.profiles = []; end
if ~isfield(opts,'layer_dir') || isempty(opts.layer_dir)
  opts.layer_dir = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_layer';
end
if ~isfield(opts,'bed_dir') || isempty(opts.bed_dir)
  opts.bed_dir = '/kucresis/scratch/hoffmana_sta/vvel/opr_out/accum/2022_Antarctica_Ground/CSARP_layer_bed';
end
if ~isfield(opts,'net_dir') || isempty(opts.net_dir)
  opts.net_dir = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel_net';
end
if ~isfield(opts,'tide') || isempty(opts.tide), opts.tide = 'cats'; end
assert(any(strcmp(opts.tide, {'cats','linemean'})), 'opts.tide must be ''cats'' or ''linemean''');
if ~isfield(opts,'max_pass_sigma') || isempty(opts.max_pass_sigma), opts.max_pass_sigma = 4; end
if ~isfield(opts,'jackknife') || isempty(opts.jackknife), opts.jackknife = true; end
if ~isfield(opts,'strain_source') || isempty(opts.strain_source), opts.strain_source = 'network'; end
assert(any(strcmp(opts.strain_source, {'network','stack'})), 'opts.strain_source must be ''network'' or ''stack''');
if ~isfield(opts,'stack_file') || isempty(opts.stack_file), opts.stack_file = 'tidal_stack_nozc.mat'; end
if ~isfield(opts,'stack_depth') || isempty(opts.stack_depth), opts.stack_depth = 100; end
if ~isfield(opts,'stack_ref') || isempty(opts.stack_ref), opts.stack_ref = 5; end
if ~isfield(opts,'stack_offset') || isempty(opts.stack_offset), opts.stack_offset = true; end

MIN_PASS = 5;                  % passes needed to regress a block on the tide
E_GRID   = logspace(log10(0.05e9), log10(30e9), 61);
X0_GRID  = -6000:250:1000;     % clamp position, seaward coordinates
% The strain observable's interval and nuisance term, echoed to callers
% via R.strain_spec: the network admittance is read at 100 m from the
% surface; the stack from its own phase reference down to stack_depth
if strcmp(opts.strain_source, 'stack')
  SSPEC = struct('ref_depth', opts.stack_depth, 'top_depth', opts.stack_ref, 'fit_offset', opts.stack_offset, ...
    'avg_width', []);                  % set from the stack's own blocks when it is loaded
else
  SSPEC = struct('ref_depth', 100, 'top_depth', 0, 'fit_offset', false, 'avg_width', []);
end

% THICKNESS. The primary case is the TRACKED BED - the OPR layer_tracker
% (Viterbi, the toolbox's own ACCUM bottom configuration, multiple
% suppression off because this is a ground platform) run over the four
% main-pass segments into opts.bed_dir. It is the only source that both
% resolves h(x) along the line and is validated: at GA10, 64-65 m off the
% line, it lands within 0.026 us (~2 m) of the ApRES bed on all four
% segments.
%
% 'BedMachine' is NOT a radar pick and is kept only as a comparison: the
% archive's bottom_mc is the BedMachine v3 thickness grid inserted as a
% pseudo-layer (runOpsInsertLayer.m; _top/_bot are the grid -/+ errbed).
% At GA10 it is 0.54 us (~45 m) shallower than both ApRES and the tracked
% return, and fitting against it would make E* partly a function of a
% model grid.
%
% The two constants bracket the plausible range and keep the h^-3
% exponent check meaningful: 265 is the mean of the two stable ApRES bed
% picks, 295 the depth where the coherent column ends in the multipass
% data.
% The third column says whether the ENGLACIAL STRAIN ADMITTANCE joins the
% fit. The primary case is joint - the strain is an absolute observable of
% the same beam, and it supplies the amplitude constraint the shape alone
% cannot (see the opts.strain block in vdef.invertElasticModulus). The
% joint fit re-weights the two datasets from their own residuals
% (opts.rescale): at face value the strain rows misfit at reduced
% chi-squared 2-4 against 0.1-1 for the shape, and left at that weight
% they set the answer while contributing far less information than their
% sigmas claim. The 'a(x) only' row is the same thickness without the
% strain, so the strain's contribution is a printed difference rather than
% an assertion. The two 'strain only' rows are the converse: the beam
% fitted to the englacial strain ALONE, no surface deflection, with the
% amplitude fixed at the freely floating far field (1 m per m of tide) or
% fitted - the radar's own estimate of E*. On one line the thin-plate
% bending signal at 100 m peaks near 3 mm/m against 1.6-3.8 mm/m per block,
% and noise reads as curvature, so expect these to scatter LOW (median
% 1.0 GPa for a 2.4 truth at 2 mm/m on this window, test_beam_flexure
% check 17) with intervals that still cover.
BED_TR = struct('kind','bed','dir',opts.bed_dir,   'layer','bottom');
BED_BM = struct('kind','bed','dir',opts.layer_dir, 'layer','bottom_mc');
H_CASES = { ...
  'tracked bed', BED_TR, true;  ...
  'a(x) only',   BED_TR, false; ...
  'BedMachine',  BED_BM, true;  ...
  'ApRES 265 m', 265,    true;  ...
  'radar 295 m', 295,    true;  ...
  'strain only', BED_TR, 'strain'; ...
  'strain, amp free', BED_TR, 'strain_free' };

%% Flexure profiles
topts = struct('block', opts.block, 'min_pass', MIN_PASS, ...
               'max_pass_sigma', opts.max_pass_sigma, 'tide', opts.tide);
if ~isempty(opts.profiles)
  L = opts.profiles;
else
  L = [];
  fprintf('=== surface tidal admittance: passes (tide regressor: %s) ===\n', opts.tide);
  for n = 1:numel(opts.products)
    S = gps_profile(opts.products{n}, opts.mp_dir, topts);
    if isempty(S), continue; end
    if isempty(L), L = S; else, L(end+1) = S; end %#ok<AGROW>
  end
end
assert(~isempty(L), 'no lines loaded from %s', opts.mp_dir);

%% Englacial strain admittance per line
% Loaded once and converted into each line's own seaward frame with the
% SAME flip parameters gps_profile used, so the two observables of a line
% cannot disagree about where its blocks are - and fitted against the
% SAME tide, with the same passes gated out. A line without usable
% products just runs shape-only.
for f = {'s_x','s_y','s_sig','s_epoch_x','s_rows','s_tday','s_tide','pass_ok','tide'}
  if ~isfield(L, f{1}), L(1).(f{1}) = []; end
end
if strcmp(opts.strain_source, 'stack')
  fs = fullfile(vdef.figureDir(), opts.stack_file);
  assert(exist(fs, 'file') == 2, 'stack file %s not found', fs);
  S = load(fs); ST = S.OUT;
  assert(isfield(ST, 'a_q'), '%s has no phase-lag fit (a_q): rerun tidal_stack', fs);
  % the stack has its own blocks (500 m), whatever block a(x) is built on
  SSPEC.avg_width = median(diff(sort(ST(1).along(:))));
  fprintf('\n=== englacial strain from the coherent stack (%s), %d m to %d m ===\n', ...
    opts.stack_file, SSPEC.top_depth, SSPEC.ref_depth);
  for i = 1:numel(L)
    O = ST(strcmp({ST.name}, L(i).name));
    if isempty(O), continue; end
    [~, iz] = min(abs(O.zsel - SSPEC.ref_depth));
    xs = L(i).x_flip_sign * (O.along(:) - L(i).x_flip_ref);
    sy = O.a_q(iz, :).'; sg = O.a_q_sd(iz, :).';
    okg = isfinite(xs) & isfinite(sy) & isfinite(sg) & sg > 0;
    [L(i).s_x, is] = sort(xs(okg));
    sy = sy(okg); L(i).s_y = sy(is); sg = sg(okg); L(i).s_sig = sg(is);
    fprintf('%-5s %2d blocks, %2d pairs, a_q(%.0f m) %+.1f to %+.1f mm/m (median sigma %.1f)\n', ...
      short(L(i).name), nnz(okg), O.npair, O.zsel(iz), 1e3*min(L(i).s_y), 1e3*max(L(i).s_y), ...
      1e3*median(L(i).s_sig));
  end
elseif exist(opts.net_dir,'dir')
  fprintf('\n=== englacial strain admittance (from %s) ===\n', opts.net_dir);
  for i = 1:numel(L)
    sopts = struct('ref_depth', SSPEC.ref_depth);
    if isfield(L(i),'tide') && ~isempty(L(i).tide) && strcmp(opts.tide,'cats')
      sopts.tide = L(i).tide; sopts.pass_ok = L(i).pass_ok;
    end
    G = load_strain_admittance(L(i).name, opts.net_dir, opts.mp_dir, sopts);
    if isempty(G), continue; end
    xs = L(i).x_flip_sign * (G.along - L(i).x_flip_ref);
    okg = isfinite(xs) & isfinite(G.adm) & isfinite(G.adm_std) & G.adm_std > 0;
    rows = find(okg);
    [L(i).s_x, is] = sort(xs(okg));
    sy = G.adm(okg);     L(i).s_y   = sy(is);
    sg = G.adm_std(okg); L(i).s_sig = sg(is);
    % Enough to refit the admittance on a subset of passes (the jackknife)
    % without touching the products again: the per-epoch column change of
    % every block, in the sorted order, and the axes it was fitted on.
    L(i).s_epoch_x = G.epoch_x(rows(is), :);
    L(i).s_rows    = rows(is);
    L(i).s_tday    = G.tday;
    L(i).s_tide    = G.tide;
    fprintf('%-5s %2d blocks, %2d pairs, dh(100 m) %+.1f to %+.1f mm/m (median sigma %.1f)\n', ...
      short(L(i).name), nnz(okg), G.n_pair, 1e3*min(L(i).s_y), ...
      1e3*max(L(i).s_y), 1e3*median(L(i).s_sig));
  end
else
  fprintf('\nno vvel network products at %s - all cases run shape-only\n', opts.net_dir);
end

fprintf('\n=== surface tidal admittance a(x), %.0f m blocks ===\n', opts.block*2.5);
fprintf('%-5s %6s %10s %10s %10s %9s %6s\n', ...
  'line','nblk','a north','a south','span km','tide m','passes');
for i = 1:numel(L)
  ok = isfinite(L(i).a);
  fprintf('%-5s %6d %10.3f %10.3f %10.2f %9.3f %6d\n', ...
    short(L(i).name), nnz(ok), L(i).a(find(ok,1,'first')), ...
    L(i).a(find(ok,1,'last')), ...
    (max(L(i).x_sea(ok))-min(L(i).x_sea(ok)))/1e3, L(i).tide_range, nnz(L(i).pass_ok));
end
fprintf(['blocks are ordered from the north (grounding) end, so x_sea and ' ...
         'a both rise\ndown each line. a is in metres of surface per metre of %s.\n'], ...
  tide_units(opts.tide));

%% Inversion, one fit per line per thickness case
fprintf('\n=== effective Young''s modulus ===\n');
fprintf('%-13s %-5s %7s %7s %7s %8s %7s %6s %6s %5s %5s\n', ...
  'thickness','line','E GPa','lo','hi','clamp km','a_ff','X2a','X2s','scale','flag');
nC   = size(H_CASES,1);
Efit = nan(nC, numel(L));
Rall = cell(nC, numel(L));
for c = 1:nC
  for i = 1:numel(L)
    ok = isfinite(L(i).a) & isfinite(L(i).a_std) & L(i).a_std > 0;
    if nnz(ok) < 6
      fprintf('%-13s %-5s   too few usable blocks (%d)\n', ...
        H_CASES{c,1}, short(L(i).name), nnz(ok));
      continue;
    end
    hspec = thickness_spec(H_CASES{c,2}, L(i), opts);
    iopts = fit_opts(L(i), hspec, ok, H_CASES{c,3}, opts.block, E_GRID, X0_GRID, SSPEC);
    if ischar(H_CASES{c,3})                       % strain only: no deflection rows
      if ~isfield(iopts, 'strain'), continue; end
      R = vdef.invertElasticModulus([], [], iopts);
    else
      R = vdef.invertElasticModulus(L(i).x_sea(ok), L(i).a(ok), iopts);
    end
    Efit(c,i) = R.E;
    Rall{c,i} = R;
    fprintf('%-13s %-5s %7.2f %7.2f %7.2f %8.2f %7.2f %6.2f %s %5s %5s\n', ...
      H_CASES{c,1}, short(L(i).name), R.E/1e9, R.E_lo/1e9, R.E_hi/1e9, ...
      R.x0/1e3, R.amplitude, chi2a_of(R), x2s_of(R), scale_of(R), flags_of(R));
  end
end
% The strain-only cases combined over the lines. The lines share no
% passes, so their profiled misfits (each in units of its own variance, the
% clamp profiled per line) add: one E* for the shelf from the radar alone.
COMB = struct('case', {}, 'E', {}, 'E_lo', {}, 'E_hi', {}, 'n', {});
for c = 1:nC
  if ~ischar(H_CASES{c,3}), continue; end
  Jsum = zeros(size(E_GRID)); n = 0;
  for i = 1:numel(L)
    Ri = Rall{c,i};
    if isempty(Ri) || ~isfinite(Ri.s2), continue; end
    Jsum = Jsum + (Ri.J_profile - min(Ri.J_profile)) / Ri.s2; n = n + 1;
  end
  if n < 2, continue; end
  [jm, im] = min(Jsum); lE = log10(E_GRID);
  lo = NaN; hi = NaN;
  k = find(Jsum(1:im) > jm + 1, 1, 'last');
  if ~isempty(k), lo = 10^interp1(Jsum([k k+1]), lE([k k+1]), jm + 1); end
  k = im - 1 + find(Jsum(im:end) > jm + 1, 1, 'first');
  if ~isempty(k), hi = 10^interp1(Jsum([k-1 k]), lE([k-1 k]), jm + 1); end
  COMB(end+1) = struct('case', H_CASES{c,1}, 'E', E_GRID(im), 'E_lo', lo, 'E_hi', hi, 'n', n); %#ok<AGROW>
  fprintf('%-16s %d lines combined: E* = %.2f GPa [%.2f, %.2f] (grid-resolved, %.2f decade steps)\n', ...
    H_CASES{c,1}, n, E_GRID(im)/1e9, lo/1e9, hi/1e9, mean(diff(lE)));
end
fprintf(['flags: E minimum on a grid edge, M rival minima in the misfit, ' ...
         'P under 3 flexural lengths of seaward pad\n']);
fprintf(['a_ff is the far-field admittance the fit implies, per metre of %s. ' ...
         'X2a and X2s\nare the reduced chi-squareds of the shape and strain ' ...
         'datasets against their INPUT\nsigmas; scale is the factor the joint ' ...
         'fit applied to the strain sigmas to bring its\nweight in line with ' ...
         'its scatter. Intervals are delta-chi-squared = 1 in the sigmas\n' ...
         'used, never narrower than the input sigmas imply.\n'], tide_units(opts.tide));

%% Jackknife over passes, primary case
% The error the formal interval cannot see. Every block of a line shares
% its passes, so a per-pass height error shifts the whole profile at once;
% the fit absorbs most of that in the amplitude and the residuals come out
% SMALLER than the sigmas (shape chi2 0.1-0.2 on GL3/GL4), which is not
% better data. Leaving one pass out at a time, recomputing a(x) AND the
% strain admittance without it, and refitting, measures what the passes
% actually pin down. The refit keeps the primary fit's dataset weights
% (no re-scaling per pass), so the jackknife varies the data alone.
JK = struct('name', {}, 'E', {}, 'se_dec', {}, 'E_lo', {}, 'E_hi', {}, 'n', {});
if opts.jackknife && isempty(opts.profiles)
  fprintf('\n=== leave-one-pass-out jackknife, %s ===\n', H_CASES{1,1});
  fprintf('%-5s %7s %15s %15s %7s %3s  %s\n', 'line','E GPa','formal','jackknife','SE dec','n','per-pass E*');
  for i = 1:numel(L)
    R0 = Rall{1,i};
    if isempty(R0), continue; end
    hspec = thickness_spec(H_CASES{1,2}, L(i), opts);
    Np = numel(L(i).pass_ok); Ek = nan(1, Np);
    for k = find(L(i).pass_ok)
      Lk = drop_pass(L(i), k, topts);
      if isempty(Lk), continue; end
      ok = isfinite(Lk.a) & isfinite(Lk.a_std) & Lk.a_std > 0;
      if nnz(ok) < 6, continue; end
      iopts = fit_opts(Lk, hspec, ok, H_CASES{1,3}, opts.block, E_GRID, X0_GRID, SSPEC);
      iopts.rescale = false;
      iopts.sigma = iopts.sigma * R0.sigma_scale_shape;
      if isfield(iopts,'strain'), iopts.strain.sigma = iopts.strain.sigma * R0.sigma_scale_strain; end
      Rk = vdef.invertElasticModulus(Lk.x_sea(ok), Lk.a(ok), iopts);
      if Rk.interior
        Ek(k) = Rk.E;
      else
        fprintf('%-5s without pass %d: fit ran to a grid edge, left out of the jackknife\n', ...
          short(L(i).name), k);
      end
    end
    g = isfinite(Ek); m = nnz(g);
    if m < 3, continue; end
    le = log10(Ek(g));
    se = sqrt((m-1)/m * sum((le - mean(le)).^2));
    JK(end+1) = struct('name', L(i).name, 'E', Ek, 'se_dec', se, ...
      'E_lo', R0.E*10^(-se), 'E_hi', R0.E*10^se, 'n', m); %#ok<AGROW>
    Rall{1,i}.E_lo_jack = JK(end).E_lo;
    Rall{1,i}.E_hi_jack = JK(end).E_hi;
    fprintf('%-5s %7.2f [%5.2f, %6.2f] [%5.2f, %6.2f] %7.3f %3d  %s\n', ...
      short(L(i).name), R0.E/1e9, R0.E_lo/1e9, R0.E_hi/1e9, JK(end).E_lo/1e9, ...
      JK(end).E_hi/1e9, se, m, sprintf('%.2f ', Ek(g)/1e9));
  end
  fprintf(['the jackknife interval is E* x/ 10^SE; where it is wider than the ' ...
           'formal one the\ndifference is what the passes share.\n']);
end

%% What the numbers actually support
fprintf('\n=== spread ===\n');
for c = 1:nC
  e = Efit(c,:); e = e(isfinite(e));
  if numel(e) < 2, continue; end
  fprintf('%-13s %d lines: %.2f to %.2f GPa, mean %.2f, sd %.2f (%.0f%%)\n', ...
    H_CASES{c,1}, numel(e), min(e)/1e9, max(e)/1e9, mean(e)/1e9, ...
    std(e)/1e9, 100*std(e)/mean(e));
end
eall = Efit(isfinite(Efit));
if ~isempty(eall)
  fprintf(['every line and thickness case: %.2f to %.2f GPa, mean %.2f, ' ...
           'sd %.2f\n'], min(eall)/1e9, max(eall)/1e9, mean(eall)/1e9, std(eall)/1e9);
end

% The thickness lever, measured rather than asserted. E* should scale as
% h^-3 between the two constant-thickness cases; printing the achieved
% exponent checks that the two fits found the same rigidity rather than
% wandering off to a different clamp position.
% Cases 4 and 5 are the two CONSTANT thicknesses; the exponent only means
% h^-3 between fits that differ in nothing but a scalar h.
if all(isfinite(Efit(4,:))) && all(isfinite(Efit(5,:)))
  ratio = mean(Efit(5,:)) / mean(Efit(4,:));
  fprintf(['h 265 -> 295 m moves mean E* by x%.3f, an exponent of %.2f ' ...
           '(h^-3 predicts -3.00)\n'], ratio, log(ratio)/log(295/265));
end

fprintf(['\nThe intervals above are per-line. On a window this short they are ' ...
         'wide because\nE* trades off against a clamp the data never see ' ...
         '(scripts/test_beam_flexure.m,\ncheck 7). The legs share their passes, ' ...
         'so agreement between them does not test\nthe per-pass height error; ' ...
         'the jackknife does. Quote it alongside the spread.\n']);

% The averaging gate the result figures share (elasticity_results,
% elasticity_map): the data bound E* from above AND the beam describes the
% line's own shape data (reduced chi-squared below 3). Criterion, not
% name, so the gate tracks the data across rebuilds.
con = false(nC, numel(L));
for c = 1:nC
  for i = 1:numel(L)
    Ri = Rall{c,i};
    if isempty(Ri), continue; end
    x2a = chi2a_of(Ri);
    con(c,i) = isfinite(Ri.E_hi) && isfinite(x2a) && x2a < 3;
  end
end

OUT = struct('lines', L, 'fits', {Rall}, 'E', Efit, ...
             'cases', {H_CASES(:,1)}, 'block', opts.block, ...
             'constrained', con, 'jackknife', JK, 'tide', opts.tide, ...
             'strain_source', opts.strain_source, 'strain_spec', SSPEC, 'strain_only', {COMB});   % braces: an empty COMB would empty OUT

end

%% ========================================================================
function iopts = fit_opts(Li, hspec, ok, mode, block, E_GRID, X0_GRID, SSPEC)
%FIT_OPTS The inversion options for one line, built in one place so the
%   main fits and the jackknife refits cannot differ in anything but data.
% a(x) is a BLOCK MEAN, so the model has to be block-averaged too. The
% deflection is curved on the scale of a flexural length and a 500 m
% block mean is not the value at the block centre; left uncorrected that
% bias runs an order of magnitude above the formal error on a(x).
% mode: true (joint), false (a(x) only), 'strain' (strain alone, amplitude
% fixed at the floating far field) or 'strain_free' (amplitude fitted)
iopts = struct('h', hspec, 'sigma', Li.a_std(ok), 'avg_width', block*2.5, ...
               'E_grid', E_GRID, 'x0_grid', X0_GRID);
strain_only = ischar(mode);
if strain_only
  iopts = rmfield(iopts, 'sigma');
  % the clamp is searched against the strain blocks alone; the grid must sit
  % landward of the last of them
  iopts.x0_grid = X0_GRID(X0_GRID < max(Li.s_x));
  if strcmp(mode, 'strain'), iopts.amplitude = 1; end
end
if (strain_only || mode) && numel(Li.s_x) >= 3
  iopts.strain = struct('x', Li.s_x, 'y', Li.s_y, 'sigma', Li.s_sig, ...
                        'ref_depth', SSPEC.ref_depth, 'top_depth', SSPEC.top_depth, ...
                        'fit_offset', SSPEC.fit_offset, 'avg_width', block*2.5);
  if ~isempty(SSPEC.avg_width), iopts.strain.avg_width = SSPEC.avg_width; end
  iopts.rescale = ~strain_only;
end
end

%% ========================================================================
function Lk = drop_pass(Li, k, topts)
%DROP_PASS The line's two observables recomputed without pass k.
%   a(x) from the stored heights with the pass's tide set NaN (the gate is
%   not re-run: the accepted set is the full fit's, minus this pass), and
%   the strain admittance refitted on the stored per-epoch column change.
Lk = Li;
tide = Li.tide; tide(~Li.pass_ok) = NaN; tide(k) = NaN;
if strcmp(topts.tide, 'linemean'), tide = []; end
Zk = Li.Z; Zk(:, k) = NaN;
if isempty(tide) && nnz(any(isfinite(Zk),1)) < topts.min_pass, Lk = []; return; end
% The per-SAMPLE geometry, not the per-block one the line struct exposes
% as along/lat/lon: the blocks are rebuilt from the samples here.
P = profile_from_Z(Zk, tide, Li.along_full, Li.lat_full, Li.lon_full, ...
      setfield(topts, 'max_pass_sigma', Inf), Li.name); %#ok<SFLD>
Lk.a = P.a; Lk.a_std = P.a_std; Lk.x_sea = P.x_sea;
Lk.pass_ok = P.pass_ok;
if ~isempty(Li.s_epoch_x)
  st = Li.s_tide; st(k) = NaN;
  A  = vdef.fitTideAdmittance(Li.s_epoch_x, Li.s_tday, st);
  Lk.s_y = A.admittance(:); Lk.s_sig = A.admittance_std(:);
  okk = isfinite(Lk.s_y) & isfinite(Lk.s_sig) & Lk.s_sig > 0;
  Lk.s_x = Li.s_x(okk); Lk.s_y = Lk.s_y(okk); Lk.s_sig = Lk.s_sig(okk);
end
end

%% ========================================================================
function s = tide_units(kind)
if strcmp(kind, 'cats'), s = 'CATS2008 tide'; else, s = 'line-mean height'; end
end

function v = x2s_of(R)
v = '     -';
if R.has_strain, v = sprintf('%6.2f', R.chi2_strain); end
end

function v = scale_of(R)
v = '    -';
if R.has_strain && isfield(R,'sigma_scale_strain'), v = sprintf('%5.2f', R.sigma_scale_strain); end
end

function flag = flags_of(R)
flag = '';
if ~R.interior,       flag = [flag 'E']; end
if R.n_local_min > 1, flag = [flag 'M']; end
if R.n_lambda < 3,    flag = [flag 'P']; end
end

%% ========================================================================
function S = gps_profile(pn, mp_dir, topts)
%GPS_PROFILE Surface tidal admittance a(x) and its standard error.
%   Loads the per-pass platform heights of one multipass product and hands
%   them to vdef.surfaceAdmittance against the CATS2008 tide at each
%   pass's mid-time (or the line mean, for the legacy estimator), then
%   orients the blocks seaward. Prints the per-pass table: tide, line
%   mean, residual, and whether the gate kept the pass.
S = [];
fn = fullfile(mp_dir, sprintf('%s_multipass03.mat', pn));
if ~exist(fn,'file')
  fprintf('missing %s\n', fn);
  return;
end
D = load(fn, 'pass','param_multipass');
main_idx = D.param_multipass.multipass.baseline_master_idx;
Np = numel(D.pass);
Nx = numel(D.pass(main_idx).ref_z);

Z = nan(Nx, Np); tmid = nan(1, Np); segs = repmat({''}, 1, Np);
for k = 1:Np
  z = D.pass(k).ref_z(:);
  if numel(z) ~= Nx, continue; end
  Z(:,k)  = z;
  tmid(k) = mean(D.pass(k).gps_time, 'omitnan');
  if isfield(D.pass(k),'param_pass') && isfield(D.pass(k).param_pass,'day_seg')
    segs{k} = char(D.pass(k).param_pass.day_seg);
  end
end
along = D.pass(main_idx).along_track(:);
lat   = D.pass(main_idx).lat(:);
lon   = D.pass(main_idx).lon(:);
day_seg = segs{main_idx};
clear D;

if strcmp(topts.tide, 'cats')
  tide = cats2008_tide(tmid);
  if any(isfinite(tmid) & ~isfinite(tide))
    fprintf('%-5s %d pass(es) fall outside the CATS2008 window and are dropped\n', ...
      short(pn), nnz(isfinite(tmid) & ~isfinite(tide)));
  end
else
  tide = [];
end

P = profile_from_Z(Z, tide, along, lat, lon, topts, pn);

fprintf('%-5s %2d passes, %2d kept, %2d blocks, corr(a, seaward x) = %+.2f', ...
  short(pn), Np, nnz(P.pass_ok), nnz(isfinite(P.a)), P.corr);
if ~isempty(tide)
  fprintf(', line mean = %.2f x tide, residual rms %.3f m\n', P.abar, ...
    sqrt(mean(P.pass_resid(P.pass_ok).^2)));
  for k = 1:Np
    if ~isfinite(tmid(k)), continue; end
    tag = ''; if ~P.pass_ok(k), tag = '  REJECTED by the pass gate'; end
    if k == main_idx, tag = [tag '  (main)']; end %#ok<AGROW>
    when = datestr(datenum(1970,1,1) + tmid(k)/86400, 'mm-dd HH:MM');   %#ok<DATST,DATNM> datetime is not in Octave
    fprintf('      %2d %-12s %s  tide %+.3f  line mean %+.3f  resid %+.3f%s\n', ...
      k, segs{k}, when, ...
      tide(k), mean(Z(:,k), 'omitnan'), P.pass_resid(k), tag);
  end
else
  fprintf('\n');
end

S = struct('name', pn, 'a', P.a, 'a_std', P.a_std, 'x_sea', P.x_sea, ...
  'along', P.along, 'lat', P.lat, 'lon', P.lon, 'day_seg', day_seg, ...
  'x_flip_sign', P.sgn, 'x_flip_ref', P.ref, ...
  'tide_range', P.tide_range, 'n_pass', Np, ...
  'Z', Z, 'tide', P.tide, 'pass_ok', P.pass_ok, 'pass_resid', P.pass_resid, ...
  'abar', P.abar, 'tmid', tmid, 'segs', {segs}, ...
  'along_full', along, 'lat_full', lat, 'lon_full', lon);
end

%% ========================================================================
function P = profile_from_Z(Z, tide, along, lat, lon, topts, pn)
%PROFILE_FROM_Z Blocks, regression, and seaward orientation.
%   Shared by the full fit and the jackknife so a dropped pass changes
%   nothing but the pass set.
A = vdef.surfaceAdmittance(Z, tide, struct('block', topts.block, ...
      'min_pass', topts.min_pass, 'max_pass_sigma', topts.max_pass_sigma));
nb = numel(A.a);
xb = nan(nb,1); latb = nan(nb,1); lonb = nan(nb,1);
for b = 1:nb
  idx = A.cols{b};
  xb(b)   = mean(along(idx), 'omitnan');
  latb(b) = mean(lat(idx),   'omitnan');
  lonb(b) = mean(lon(idx),   'omitnan');
end
a = A.a; a_std = A.a_std;

% Orientation is set by GEOGRAPHY, not by the data. Grounding is at the
% NORTH end of this line, so seaward is south and x_sea counts from the
% north end. Legs are walked out and back, so along_track runs north on
% some products and south on others and the sign has to be read off the
% latitudes. Deciding it from the sign of the a(x) trend instead would let
% a noisy line silently invert its own coordinates and still "fit".
ok = isfinite(a) & isfinite(latb);
assert(any(ok), '%s: no block yielded an admittance', pn);
first = find(ok,1,'first'); last = find(ok,1,'last');
if latb(first) > latb(last)
  sgn = +1;  ref = xb(first);       % the first block is the northern end
else
  sgn = -1;  ref = xb(last);
end
x_sea = sgn * (xb - ref);

% Half the legs are walked north and half south, so for half of them x_sea
% now runs backwards through the arrays. Sort every per-block field together
% into increasing x_sea, which is the order the inversion requires and the
% order that makes "first block" mean the grounding end on every line.
[x_sea, is] = sort(x_sea);
a = a(is); a_std = a_std(is);
xb = xb(is); latb = latb(is); lonb = lonb(is);
ok = ok(is);

% Now check the data agree with the geography. They should: a rises seaward.
sel = ok & isfinite(x_sea);
u = x_sea(sel) - mean(x_sea(sel)); v = a(sel) - mean(a(sel));
den = sqrt((u.'*u) * (v.'*v));
if den <= 0, cc = NaN; else, cc = (u.'*v)/den; end
assert(cc > 0, ['%s: a(x) FALLS seaward, which contradicts a grounding ' ...
  'line at the north end (corr %+.2f). Sort out the orientation before ' ...
  'fitting anything to it.'], pn, cc);

tk = A.tide(A.pass_ok);
P = struct('a', a, 'a_std', a_std, 'x_sea', x_sea, 'along', xb, 'lat', latb, ...
  'lon', lonb, 'sgn', sgn, 'ref', ref, 'corr', cc, 'tide', A.tide, ...
  'pass_ok', A.pass_ok, 'pass_resid', A.pass_resid, 'abar', A.abar, ...
  'tide_range', max(tk) - min(tk));
end

%% ========================================================================
function s = short(name)
s = strrep(name, 'EAGER_2022_', '');
end

%% ========================================================================
function v = chi2a_of(R)
%CHI2A_OF The shape dataset's reduced chi-squared, whichever fit produced it.
%   Joint fits report it per dataset; shape-only fits report the overall
%   chi2red, which for them is the same quantity.
if R.has_strain && isfinite(R.chi2_shape)
  v = R.chi2_shape;
else
  v = R.chi2red;
end
end

%% ========================================================================
function hspec = thickness_spec(spec, S, opts)
%THICKNESS_SPEC Turn a thickness case into something the inversion accepts.
%   A scalar passes straight through; a bed spec (struct with kind='bed')
%   reads that layer product. The earlier two-point ApRES h(x) case is
%   retired: the tracked bed supersedes it, resolving h(x) at every block
%   instead of interpolating between two sites.
if isstruct(spec) && isfield(spec,'kind') && strcmp(spec.kind,'bed')
  hspec = bed_profile(S, spec, opts);
  return;
end
assert(isnumeric(spec) && isscalar(spec), 'unrecognised thickness spec');
hspec = spec;
end

%% ========================================================================
function hspec = bed_profile(S, spec, opts)
%BED_PROFILE Ice thickness along the line, from a layer-product bed.
%   Reads layer spec.layer out of the layer product spec.dir for the
%   day_seg of this line's MAIN pass, subtracts the surface layer's
%   traveltime (from the archive CSARP_layer, since the tracked product
%   carries only the bottom), maps the result onto the line's own blocks,
%   and converts traveltime below the surface to depth through the firn
%   column. The surface term is ~3 ns here - the layer twtt axis is
%   already surface-referenced to within a bin - so this is rigor rather
%   than a correction that moves anything.
%
%   WHY NOT THE LAYERS IN THE MULTIPASS PRODUCT. They are there - two of
%   them, named surface and bottom - but bottom is identically zero in
%   every pass of every product. multipass loads the standard surface and
%   bottom pair, and this season's bed pick is not stored under that name.
%
%   MAPPING. Each day_seg walks the same line four times out and back, so
%   the layer file covers about 21 km and the multipass leg covers 4.8 km
%   of it. Points are matched by NEAREST NEIGHBOUR in a local metric frame
%   rather than by index or by gps_time, which is what picks the right leg
%   out of the four: the legs are 120-213 m apart while the layer sampling
%   is ~15 m along track, so the nearest point is always on the correct
%   leg. The block MEDIAN is then taken over an along/across-track CAPSULE
%   about the block centre - half a block (~250 m) along the line's own
%   bearing but only MAX_CROSS (~60 m) across it - rather than a Euclidean
%   disc, whose 250 m radius would reach the adjacent leg corridors and
%   mix their bed picks into this leg's median. Achieved match distances
%   are a few metres, and a block whose nearest point sits further away
%   than MAX_MATCH is dropped rather than filled.
%
%   SCREENING. About 13% of the picks fail outright, mostly at the
%   turnarounds and toward the grounding end where the basal return
%   weakens, and they fail SHALLOW - the 1st and 5th percentiles are 64 and
%   93 m against a 239 m median. A median-absolute-deviation cut removes
%   them; a block that loses most of its points is dropped, and
%   vdef.invertElasticModulus then holds thickness constant beyond the
%   outermost surviving block rather than extrapolating a trend into the
%   region where the pick just failed. That is an ASSUMPTION, and it lands
%   exactly where the clamp is, so it is reported per line.
MAX_MATCH = 60;      % m; blocks matched further away than this are dropped
MAX_CROSS = 60;      % m; cross-track half-width of the block capsule, well
                     % under the 120-213 m leg separation
N_MAD     = 5;       % robust outlier cut on the depth samples
MIN_FRAC  = 0.25;    % fraction of a block's samples that must survive

hspec = [];
fn = fullfile(spec.dir, S.day_seg, sprintf('Data_%s_001.mat', S.day_seg));
lf = fullfile(spec.dir, S.day_seg, sprintf('layer_%s.mat',  S.day_seg));
if ~exist(fn,'file') || ~exist(lf,'file')
  fprintf('  %-5s no layer product for %s in %s - bed thickness unavailable\n', ...
    short(S.name), S.day_seg, spec.dir);
  return;
end
T = load(fn); Q = load(lf);
nms = cellfun(@char, Q.lyr_name, 'uni', 0);
j = find(strcmp(nms, spec.layer), 1);
if isempty(j)
  fprintf('  %-5s layer "%s" not in %s (have: %s)\n', short(S.name), ...
    spec.layer, S.day_seg, strjoin(nms, ', '));
  return;
end
row = find(T.id == Q.lyr_id(j), 1);
tb  = T.twtt(row,:).';

% Surface twtt from the archive layer product, interpolated by gps_time so
% a tracked product with a different point count still subtracts cleanly.
ts = zeros(size(tb));
sfn = fullfile(opts.layer_dir, S.day_seg, sprintf('Data_%s_001.mat', S.day_seg));
slf = fullfile(opts.layer_dir, S.day_seg, sprintf('layer_%s.mat',  S.day_seg));
if exist(sfn,'file') && exist(slf,'file')
  Ts = load(sfn); Qs = load(slf);
  ns = cellfun(@char, Qs.lyr_name, 'uni', 0);
  js = find(strcmp(ns,'surface'), 1);
  if ~isempty(js)
    rs = find(Ts.id == Qs.lyr_id(js), 1);
    ts = interp1(Ts.gps_time(:), Ts.twtt(rs,:).', T.gps_time(:), ...
                 'linear', 'extrap');
  end
end

P = vdef.firnColumn(vdef.defaultParams());
h_all = vdef.depthFromTwtt(P, tb - ts);

% Robust cut on the pick, over the whole segment.
good = isfinite(h_all);
med  = median(h_all(good));
mad  = median(abs(h_all(good) - med));
tol  = N_MAD * 1.4826 * max(mad, eps);
good = good & abs(h_all - med) <= tol;

% Nearest layer point to each of this line's along-track samples.
lat0 = mean(S.lat, 'omitnan'); lon0 = mean(S.lon, 'omitnan');
bx = (S.lon - lon0) * cosd(lat0) * 111320;
by = (S.lat - lat0) * 110540;
lx = (T.lon(:) - lon0) * cosd(lat0) * 111320;
ly = (T.lat(:) - lat0) * 110540;

nb = numel(S.x_sea);
xs = nan(nb,1); hs = nan(nb,1); nkept = zeros(nb,1); dm = nan(nb,1);
for b = 1:nb
  if ~isfinite(bx(b)) || ~isfinite(by(b)), continue; end
  d2 = (lx - bx(b)).^2 + (ly - by(b)).^2;
  % The block spans opts.block samples of track; take every layer point
  % within half a block ALONG the line rather than the single nearest, so
  % the thickness is a block mean like the admittance it will be fitted
  % to - but only within MAX_CROSS ACROSS it, so the median stays on this
  % leg's corridor. The local bearing comes from the nearest other block
  % centre; a block with no finite neighbour cannot orient a capsule and
  % is dropped like an unmatched one.
  half = 0.5 * opts.block * 2.5;
  db2 = (bx - bx(b)).^2 + (by - by(b)).^2;
  db2(b) = inf; db2(~isfinite(db2)) = inf;
  [dn2, jn] = min(db2);
  if ~isfinite(dn2), continue; end
  tvx = (bx(jn) - bx(b))/sqrt(dn2); tvy = (by(jn) - by(b))/sqrt(dn2);
  dal = (lx - bx(b))*tvx + (ly - by(b))*tvy;
  dcr = -(lx - bx(b))*tvy + (ly - by(b))*tvx;
  sel = abs(dal) <= half & abs(dcr) <= MAX_CROSS;
  if ~any(sel), continue; end
  dm(b) = sqrt(min(d2));
  if dm(b) > MAX_MATCH, continue; end
  hb = h_all(sel & good);
  nkept(b) = numel(hb);
  if numel(hb) < MIN_FRAC * nnz(sel), continue; end
  xs(b) = S.x_sea(b);
  hs(b) = median(hb);
end

ok = isfinite(xs) & isfinite(hs);
fprintf(['  %-5s bed from %s/%s: %d of %d blocks kept, h %.0f-%.0f m ' ...
         '(match %.0f m max)\n'], short(S.name), S.day_seg, spec.layer, ...
  nnz(ok), nb, min(hs(ok)), max(hs(ok)), max(dm(ok)));
if nnz(ok) < 2
  fprintf('  %-5s too few bed blocks - falling back to the segment median %.0f m\n', ...
    short(S.name), med);
  hspec = med;
  return;
end
if ~all(ok)
  drop = find(~ok);
  fprintf(['  %-5s blocks dropped at x_sea %s km; thickness is HELD ' ...
           'CONSTANT there\n'], short(S.name), ...
    strjoin(cellstr(num2str(S.x_sea(drop)/1e3, '%.2f')).', ', '));
end
hspec = [xs(ok), hs(ok)];
end

