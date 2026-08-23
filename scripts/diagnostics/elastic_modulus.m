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
%   for the airborne radar thickness. Numerics live in vdef.beamFlexure
%   and vdef.invertElasticModulus.
%
%   WHY THE GPS AND NOT THE RADAR. The observable a beam model wants is
%   vertical DEFLECTION of the surface, and this survey measures that
%   directly and radar-independently: each pass's ref_z(x) is the platform
%   height on a floating shelf, so regressing it on that pass's own line
%   mean gives a(x), the local admittance of the surface to the tide. That
%   profile is the one result of this project that never depended on the
%   interferometry - a(x) falls monotonically along the line by a factor of
%   about four, and GL3 and GL4 agree on it to 0.02 despite different main
%   passes and pass sets. The radar column-strain admittance is a different
%   quantity - it goes as the CURVATURE of this profile, not the profile -
%   and it sits at the measurement's systematic floor. So the deflection is
%   what gets inverted here.
%
%   THE NORMALISATION IS NOT A PROBLEM. a(x) is normalised by a line mean,
%   not by the tide, because the survey never reaches freely floating ice
%   and the far-field amplitude is therefore not observed. That would matter
%   if the amplitude were being fitted for its own sake; it is not. The
%   inversion eliminates the amplitude in closed form and fits the SHAPE of
%   the profile, so an arbitrary rescaling leaves E* bit-identical (asserted
%   in scripts/test_beam_flexure.m, check 5). What the fitted amplitude does
%   give is a by-product worth reading: the far-field admittance the beam
%   implies, in the same line-mean units, which says how much further the
%   shelf has to run past the end of the line before it floats freely.
%
%   READ THIS BEFORE QUOTING A NUMBER. This survey is a PARTIAL WINDOW. The
%   line spans about 4.8 km inside the flexure zone and reaches neither the
%   flat grounded end nor the flat floating one, so the clamp position sits
%   outside the data and trades off against E*. Check 7 of the unit test is
%   that exact geometry on synthetic data with a known answer: the modulus
%   comes back tens of percent off, with an interval wide enough to cover
%   the truth but seven times wider than the same method achieves on a full
%   profile. So the formal intervals below carry the ill-posedness honestly,
%   and a tight one on this geometry would be the thing to distrust. What
%   they do not carry is anything the four lines disagree about, which is
%   reported separately as the spread.
%
%   AND E* IS ONLY EVER E*h^3. Flexure constrains the rigidity, so the
%   modulus is no better known than the thickness cubed. Three thickness
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
%                 shape fit's amplitude - an ABSOLUTE observable, so it
%                 supplies the D^(-1/2) amplitude constraint that the
%                 line-mean normalised a(x) cannot. Missing dir = all
%                 cases run shape-only.
%     .profiles   inject a(x) profiles directly and skip the file load,
%                 for testing the inversion path without the products
%
%   Returns OUT.lines (the a(x) profiles) and OUT.fits (one inversion per
%   thickness case per line), so the figure reuses the numbers rather than
%   refitting them.
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "addpath('<code>/scripts/diagnostics'); elastic_modulus"

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));                       % +vdef

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

MIN_PASS = 5;                  % passes needed to regress a block on the tide
E_GRID   = logspace(log10(0.05e9), log10(30e9), 61);
X0_GRID  = -6000:250:1000;     % clamp position, seaward coordinates

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
% the same beam, and it supplies the amplitude constraint the line-mean
% normalised a(x) cannot (see the opts.strain block in
% vdef.invertElasticModulus). The 'a(x) only' row is the same thickness
% without it, so the strain's contribution is a printed difference rather
% than an assertion.
BED_TR = struct('kind','bed','dir',opts.bed_dir,   'layer','bottom');
BED_BM = struct('kind','bed','dir',opts.layer_dir, 'layer','bottom_mc');
H_CASES = { ...
  'tracked bed', BED_TR, true;  ...
  'a(x) only',   BED_TR, false; ...
  'BedMachine',  BED_BM, true;  ...
  'ApRES 265 m', 265,    true;  ...
  'radar 295 m', 295,    true };

%% Flexure profiles
if ~isempty(opts.profiles)
  L = opts.profiles;
else
  L = [];
  for n = 1:numel(opts.products)
    S = gps_profile(opts.products{n}, opts.mp_dir, opts.block, MIN_PASS);
    if isempty(S), continue; end
    if isempty(L), L = S; else, L(end+1) = S; end %#ok<AGROW>
  end
end
assert(~isempty(L), 'no lines loaded from %s', opts.mp_dir);

%% Englacial strain admittance per line
% Loaded once and converted into each line's own seaward frame with the
% SAME flip parameters gps_profile used, so the two observables of a line
% cannot disagree about where its blocks are. A line without usable
% products just runs shape-only.
if ~isfield(L,'s_x'), L(1).s_x = []; L(1).s_y = []; L(1).s_sig = []; end
if exist(opts.net_dir,'dir')
  fprintf('\n=== englacial strain admittance (from %s) ===\n', opts.net_dir);
  for i = 1:numel(L)
    G = load_strain_admittance(L(i).name, opts.net_dir, opts.mp_dir, ...
          struct('ref_depth', 100));
    if isempty(G), continue; end
    xs = L(i).x_flip_sign * (G.along - L(i).x_flip_ref);
    okg = isfinite(xs) & isfinite(G.adm) & isfinite(G.adm_std) & G.adm_std > 0;
    [L(i).s_x, is] = sort(xs(okg));
    sy = G.adm(okg);   L(i).s_y   = sy(is);
    sg = G.adm_std(okg); L(i).s_sig = sg(is);
    fprintf('%-5s %2d blocks, %2d pairs, dh(100 m) %+.1f to %+.1f mm/m (median sigma %.1f)\n', ...
      short(L(i).name), nnz(okg), G.n_pair, 1e3*min(L(i).s_y), ...
      1e3*max(L(i).s_y), 1e3*median(L(i).s_sig));
  end
else
  fprintf('\nno vvel network products at %s - all cases run shape-only\n', opts.net_dir);
end

fprintf('\n=== surface tidal admittance a(x), %.0f m blocks ===\n', opts.block*2.5);
fprintf('%-5s %6s %10s %10s %10s %9s\n', ...
  'line','nblk','a north','a south','span km','tide m');
for i = 1:numel(L)
  ok = isfinite(L(i).a);
  fprintf('%-5s %6d %10.3f %10.3f %10.2f %9.3f\n', ...
    short(L(i).name), nnz(ok), L(i).a(find(ok,1,'first')), ...
    L(i).a(find(ok,1,'last')), ...
    (max(L(i).x_sea(ok))-min(L(i).x_sea(ok)))/1e3, L(i).tide_range);
end
fprintf(['blocks are ordered from the north (grounding) end, so x_sea and ' ...
         'a both rise\ndown each line.\n']);

%% Inversion, one fit per line per thickness case
fprintf('\n=== effective Young''s modulus ===\n');
fprintf('%-13s %-5s %7s %7s %7s %8s %7s %6s %6s %5s\n', ...
  'thickness','line','E GPa','lo','hi','clamp km','a_ff','X2a','X2s','flag');
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
    % a(x) is a BLOCK MEAN, so the model has to be block-averaged too. The
    % deflection is curved on the scale of a flexural length and a 500 m
    % block mean is not the value at the block centre; left uncorrected that
    % bias runs an order of magnitude above the formal error on a(x).
    iopts = struct('h', hspec, 'sigma', L(i).a_std(ok), ...
                   'avg_width', opts.block*2.5, ...
                   'E_grid', E_GRID, 'x0_grid', X0_GRID);
    if H_CASES{c,3} && numel(L(i).s_x) >= 3
      iopts.strain = struct('x', L(i).s_x, 'y', L(i).s_y, ...
                            'sigma', L(i).s_sig, 'ref_depth', 100, ...
                            'avg_width', opts.block*2.5);
    end
    R = vdef.invertElasticModulus(L(i).x_sea(ok), L(i).a(ok), iopts);
    Efit(c,i) = R.E;
    Rall{c,i} = R;

    flag = '';
    if ~R.interior,       flag = [flag 'E']; end   %#ok<AGROW>
    if R.n_local_min > 1, flag = [flag 'M']; end   %#ok<AGROW>
    if R.n_lambda < 3,    flag = [flag 'P']; end   %#ok<AGROW>
    x2s = '     -';
    if R.has_strain, x2s = sprintf('%6.2f', R.chi2_strain); end
    fprintf('%-13s %-5s %7.2f %7.2f %7.2f %8.2f %7.2f %6.2f %s %5s\n', ...
      H_CASES{c,1}, short(L(i).name), R.E/1e9, R.E_lo/1e9, R.E_hi/1e9, ...
      R.x0/1e3, R.amplitude, chi2a_of(R), x2s, flag);
  end
end
fprintf(['flags: E minimum on a grid edge, M rival minima in the misfit, ' ...
         'P under 3 flexural lengths of seaward pad\n']);
fprintf(['a_ff is the far-field admittance the fit implies, in line-mean ' ...
         'units. X2a and X2s\nare the reduced chi-squareds of the shape ' ...
         'and strain datasets against their own\nsigmas - the joint E* ' ...
         'means nothing unless both are order 1.\n']);

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

fprintf(['\nThe intervals above are per-line and formal. On a window this ' ...
         'short they are wide\nbecause E* trades off against a clamp the ' ...
         'data never see, and a tight one here\nwould be the thing to ' ...
         'distrust (scripts/test_beam_flexure.m, check 7). What they do\n' ...
         'not carry is anything the lines disagree about - quote the ' ...
         'spread alongside them.\n']);

OUT = struct('lines', L, 'fits', {Rall}, 'E', Efit, 'cases', {H_CASES(:,1)});

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
%   leg. Achieved match distances are a few metres, and a block whose
%   points sit further away than MAX_MATCH is dropped rather than filled.
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
  % within half a block of it rather than the single nearest, so the
  % thickness is a block mean like the admittance it will be fitted to.
  half = 0.5 * opts.block * 2.5;
  sel  = d2 <= half^2;
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

%% ========================================================================
function S = gps_profile(pn, mp_dir, BLOCK, MIN_PASS)
%GPS_PROFILE Surface tidal admittance a(x) and its standard error.
%   Per along-track block, regress each pass's mean platform height on that
%   pass's own line-mean height - its tide - across passes:
%
%       ref_z_k(block) = a(block) * tide_k + const
%
%   a is the local admittance of the ice surface to the tide, 1 where the
%   shelf floats freely and falling to 0 into the grounding zone. The
%   regression is written out rather than left to polyfit because the fit
%   uncertainty is wanted: it weights the beam inversion, and without it a
%   block carrying five noisy passes counts as much as one carrying
%   thirteen.
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

Z = nan(Nx, Np); tide = nan(1, Np);
for k = 1:Np
  z = D.pass(k).ref_z(:);
  if numel(z) ~= Nx, continue; end
  Z(:,k)  = z;
  tide(k) = mean(z, 'omitnan');
end
along = D.pass(main_idx).along_track(:);
lat   = D.pass(main_idx).lat(:);
lon   = D.pass(main_idx).lon(:);
% The day_seg of the MAIN pass, which is the segment whose layer file
% carries the bed pick for this line's own track.
day_seg = '';
if isfield(D.pass(main_idx),'param_pass') && ...
   isfield(D.pass(main_idx).param_pass,'day_seg')
  day_seg = char(D.pass(main_idx).param_pass.day_seg);
end
clear D;

nb = floor(Nx/BLOCK);
a = nan(nb,1); a_std = nan(nb,1);
xb = nan(nb,1); latb = nan(nb,1); lonb = nan(nb,1);
for b = 1:nb
  idx = (b-1)*BLOCK+1 : b*BLOCK;
  zb  = mean(Z(idx,:), 1, 'omitnan');
  ok  = isfinite(zb) & isfinite(tide);
  n   = nnz(ok);
  if n < MIN_PASS, continue; end
  X = [ones(n,1), tide(ok).'];
  if rcond(X.'*X) < 1e-12, continue; end
  beta = X \ zb(ok).';
  res  = zb(ok).' - X*beta;
  Cov  = (sum(res.^2)/(n - 2)) * ((X.'*X) \ eye(2));
  a(b)     = beta(2);
  a_std(b) = sqrt(abs(Cov(2,2)));
  xb(b)    = mean(along(idx), 'omitnan');
  latb(b)  = mean(lat(idx),   'omitnan');
  lonb(b)  = mean(lon(idx),   'omitnan');
end

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
cc  = corr_pair(x_sea(sel), a(sel));
fprintf('%-5s %2d passes, %2d blocks, corr(a, seaward x) = %+.2f\n', ...
  short(pn), Np, nnz(sel), cc);
assert(cc > 0, ['%s: a(x) FALLS seaward, which contradicts a grounding ' ...
  'line at the north end (corr %+.2f). Sort out the orientation before ' ...
  'fitting anything to it.'], pn, cc);

S = struct('name', pn, 'a', a, 'a_std', a_std, 'x_sea', x_sea, ...
  'along', xb, 'lat', latb, 'lon', lonb, 'day_seg', day_seg, ...
  'x_flip_sign', sgn, 'x_flip_ref', ref, ...
  'tide_range', max(tide) - min(tide), 'n_pass', Np);
end

%% ========================================================================
function c = corr_pair(u, v)
u = u(:) - mean(u(:)); v = v(:) - mean(v(:));
den = sqrt((u.'*u) * (v.'*v));
if den <= 0, c = NaN; else, c = (u.'*v)/den; end
end
