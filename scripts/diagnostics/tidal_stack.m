function OUT = tidal_stack(opts)
%TIDAL_STACK Coherent all-pairs tidal response of every profile, from the interferograms.
%
%   WHAT THIS DOES. For each multipass product, every pair of passes is
%   interfered, coaligned and multilooked exactly as the vvel chain does it
%   (vdef.coalignPair, vdef.multilook), referenced to the surface bin and
%   averaged into along-track blocks as a COMPLEX field (vdef.complexBlocks),
%   and the tidal column response a(z, x) is read off all pairs at once by
%   vdef.tidalStack: counter-rotate each pair by the model phase for a
%   trial response, sum, and take the response at which they add
%   coherently. Nothing is unwrapped and no per-pair dh is formed, so this
%   is an estimator independent of the unwrap -> dh -> regression chain the
%   maps and the beam inversion rest on. It also returns the response at
%   EVERY depth, not one reference depth.
%
%   WHAT IT DOES NOT DO. The mis-registration artefact - a residual
%   tide-proportional shift mis-registering the two traces - is in the
%   phase, so a(z, x) is signal plus artefact here as everywhere. Each
%   pair's residual misalignment is MEASURED (the GPS-predicted
%   compensation error minus what coalignment removed, per block) and
%   handed to the 2-D scan as a second regressor. Where that residual is
%   not collinear with the tide the two separate; where it is, the scan is
%   a ridge, and R.rcol says which is which per block.
%
%   ORDER OF TRUST, learned the hard way on this dataset: read R.inj_err
%   first (the phase-sign self-test), then R.F (how tide-locked the phase
%   is), then rcol, and only then a. A high F with |rcol| near 1 is the
%   artefact seen directly.
%
%   opts, all optional:
%     .products    (default the four calibrated legs and EAGER_2022; all
%                   are reported - nothing is dropped for disagreeing)
%     .mp_dir, .net_dir  multipass products; all-pairs vvel products for
%                   the optional unwrapped-chain comparison at 100 m
%     .drop_passes which passes to exclude ENTIRELY, on top of the gates
%                  the chain already applies. [] (default) drops nothing;
%                  'auto' reads scripts/diagnostics/pass_quality.m's saved
%                  output and drops any pass whose median coalignment
%                  correlation against its partners falls below
%                  coalign_min_quality, or that coaligned with fewer than
%                  half of them. A struct array with fields .name and
%                  .idx names them by hand.
%
%                  WHY THIS IS NOT ALREADY DONE BY THE PAIR GATES. A pair
%                  that fails coalignment is already skipped, so dropping
%                  "the failed pairs" changes nothing. What this tests is
%                  different and worth testing: whether the pairs that a
%                  bad pass DID pass are trustworthy. A pass that cannot
%                  be registered against most of its partners has no
%                  reason to be correctly registered against the rest, and
%                  a marginal alignment produces a shift that is applied
%                  rather than refused - which is worse than a gap.
%     .tag         suffix for the saved .mat, so a variant run does not
%                  overwrite the baseline
%     .block       columns per block (default 200 = 500 m, as vvel_net)
%     .zsel        depths reported [m] (default 20:10:320; each block is
%                  masked below its own ice base, so the full thickness is
%                  used wherever the ice is thick and nothing below the base
%                  is reported where it is thin)
%     .base_margin cells within this distance above the picked base [m]
%                  are masked too, clear of the basal echo (default 15)
%     .out_dir     where the .mat goes
%
%   Prints, per product, the four estimates at 100 m block by block:
%   a_stack (1-D), a_ctrl and its artefact gain (2-D), the collinearity,
%   and the unwrapped chain's slope on the same pairs. Returns OUT with
%   everything, one element per product. Draw it with
%   scripts/figures/tidal_stack_figure.m.
%
%   THE QUADRATURE FIT. Each pair's quadrature-tide difference (pass_tide's
%   info.tide_q) goes to vdef.tidalStack as opts.dq, with the trend: OUT.a_q,
%   b_q, v_q and their jackknife errors. b = 0 is an elastic column (strain
%   in phase with the tide, strain rate a quarter period ahead of it);
%   b/a = -tan(2*pi*lag/T) otherwise. The tide here is a function of the
%   hour of day, so a lag cannot be told from a diurnal process peaking at
%   another hour; b = 0 within error is the clean outcome.
%
%   opts.ref_depth moves the phase reference from the default 50 ns (~5 m)
%   to that depth [m]; the response is then relative to that level.
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "addpath('<code>/scripts/diagnostics'); tidal_stack"

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
addpath(root); addpath(here); addpath(fullfile(root, 'opr_vvel'));
def = struct( ...
  'products', {vdef.surveyLines()}, ...
  'mp_dir', '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass', ...
  'net_dir', '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel_net', ...
  'block', 200, 'zsel', 20:10:320, 'base_margin', 15, 'drop_passes', [], 'tag', '', 'ref_depth', [], ...
  'out_dir', vdef.figureDir());
fn = fieldnames(def);
for i = 1:numel(fn)
  if ~isfield(opts, fn{i}) || (isempty(opts.(fn{i})) && ~ischar(def.(fn{i})))
    opts.(fn{i}) = def.(fn{i});
  end
end
if ~isfield(opts,'drop_passes'), opts.drop_passes = []; end
if ~isfield(opts,'tag') || isempty(opts.tag), opts.tag = ''; end
C = vdef.constants();
P = vdef.firnColumn(vdef.defaultParams());

OUT = struct('name', {}, 'x_sea', {}, 'along', {}, 'lat', {}, 'lon', {}, 'zsel', {}, 'base', {}, 'base_col', {}, 'a', {}, 'a_sd', {}, ...
  'a_sd_boot', {}, 'F', {}, 'a_t', {}, 'a_t_sd', {}, 'v_t', {}, 'F_t', {}, 'r_tt', {}, 'a2', {}, 'g2', {}, ...
  'a_q', {}, 'b_q', {}, 'v_q', {}, 'F_q', {}, 'a_q_sd', {}, 'b_q_sd', {}, 'r_qt', {}, ...
  'F2', {}, 'rcol', {}, 'c_unw', {}, 'npair', {}, 'pairs', {}, 'dt_days', {}, 'kz', {}, 'inj_err', {}, 'main_idx', {});
for n = 1:numel(opts.products)
  pn = opts.products{n}; nm = strrep(pn, 'EAGER_2022_', '');
  fn_mp = fullfile(opts.mp_dir, sprintf('%s_multipass03.mat', pn));
  if ~exist(fn_mp, 'file'), fprintf('%s: missing %s\n', nm, fn_mp); continue; end
  param = vvel_defaults(struct('vvel', struct('pass_name', pn))); o = param.vvel;
  fprintf('\n===== %s =====\n', nm); t0 = tic;
  D = load(fn_mp, 'data', 'pass', 'param_multipass');
  pm = D.param_multipass.multipass;
  zc = vdef.zmotionApplied(D.param_multipass);   % false: built without ref_z compensation, predicted shift 0
  if isfield(pm, 'pass_en_mask') && ~isempty(pm.pass_en_mask), en = find(pm.pass_en_mask);
  else, en = 1:numel(D.pass); end
  assert(numel(en) == size(D.data, 3), '%s: %d enabled passes but %d image slabs', nm, numel(en), size(D.data,3));
  main_idx = pm.baseline_master_idx; mp_ = D.pass(main_idx);
  Time = mp_.time(:); Nt = numel(Time); Nx = size(D.data, 2);
  if isfield(mp_, 'surface') && ~isempty(mp_.surface), Surface = mp_.surface(:).';
  else, Surface = mp_.layers(1).twtt_ref(:).'; end
  Surface = Surface(1:Nx);
  if isfield(mp_, 'wfs') && isfield(mp_.wfs(1), 'fc') && ~isempty(mp_.wfs(1).fc), fc = mp_.wfs(1).fc;
  else, fc = 750e6; end
  [depth, nloc] = vdef.depthFromTwtt(P, Time - mean(Surface, 'omitnan'));
  zi = arrayfun(@(z) find(depth >= z, 1, 'first'), opts.zsel, 'uni', 0);
  zi = [zi{:}]; zsel = depth(zi).';
  % the ice base on the main pass (vdef.trackBase), per column and per block
  km_ = find(en == main_idx, 1); gd = isfinite(depth) & depth >= 0;
  base_col = vdef.trackBase(10*log10(abs(double(D.data(gd, :, km_))).^2), depth(gd));
  kz = o.phase_sign * (4*pi*fc*nloc(zi)/C.c);
  starts = 1:opts.block:Nx; nb = numel(starts);
  % phase reference: the chain's default is the surface pick + ref_twtt_offset
  % (50 ns, ~5 m of firn). opts.ref_depth [m] moves it down the column - the
  % test of whether a depth-constant response is the column moving or the
  % reference being contaminated by the surface return's tail
  ref_off = o.ref_twtt_offset;
  if ~isempty(opts.ref_depth)
    tsurf = Time - mean(Surface, 'omitnan');
    g = isfinite(depth) & tsurf >= 0;                 % depth is defined from the surface down
    [du, iu] = unique(depth(g)); tg = tsurf(g);
    ref_off = interp1(du, tg(iu), opts.ref_depth, 'linear');
    assert(isfinite(ref_off), 'ref_depth %.0f m is outside the column', opts.ref_depth);
    fprintf('  phase reference at %.0f m below the surface (%.0f ns), not the default %.0f ns\n', ...
      opts.ref_depth, ref_off*1e9, o.ref_twtt_offset*1e9);
  end
  ref_bin = round(interp1(Time, 1:Nt, Surface + ref_off, 'linear', NaN));
  along = mp_.along_track(:);
  xb = arrayfun(@(s) mean(along(s:min(s+opts.block-1, Nx))), starts).';
  % block centres on the map, from the main pass's own track (the geometry
  % every slab is registered to), so a figure can draw the blocks where
  % they are instead of on a per-line seaward axis
  latb = arrayfun(@(s) mean(mp_.lat(s:min(s+opts.block-1, Nx)), 'omitnan'), starts).';
  lonb = arrayfun(@(s) mean(mp_.lon(s:min(s+opts.block-1, Nx)), 'omitnan'), starts).';
  [tide, pass_ok, tinfo] = pass_tide(pn, opts.mp_dir); tide_q = tinfo.tide_q;
  Np = numel(D.pass);
  tmid = nan(1, Np);                      % pass mid-times, for the pair interval the trend term needs
  for k = 1:Np, tmid(k) = mean(D.pass(k).gps_time, 'omitnan'); end
  drop = resolve_drop(opts.drop_passes, pn, opts.out_dir, o.coalign_min_quality);
  if ~isempty(drop)
    fprintf('  dropping pass(es) %s entirely\n', mat2str(drop));
    pass_ok(drop) = false;
  end
  fprintf('  loaded in %.0f s: %d passes (%d kept), %d columns, fc %.0f MHz, main pass %d, k(100 m) = %.1f rad/m\n', ...
    toc(t0), Np, nnz(pass_ok), Nx, fc/1e6, main_idx, interp1(zsel, kz, 100));

  I = []; W = []; DT = []; DQ = []; DTD = []; DRES = []; pairs = []; t1 = tic;
  for i = 1:Np-1
    for j = i+1:Np
      if ~pass_ok(i) || ~pass_ok(j), continue; end
      ki = find(en == i, 1); kj = find(en == j, 1);
      if isempty(ki) || isempty(kj), continue; end
      by = D.pass(j).ref_y(:) - D.pass(i).ref_y(:);
      if max(abs(by(isfinite(by)))) > o.max_baseline, continue; end
      s_ref = D.data(:,:,ki); s_sec = D.data(:,:,kj);
      if zc
        % compensated build: undo the compensation's misalignment first
        [s_sec, ca] = vdef.coalignPair(s_ref, s_sec, struct('Time', Time, 'Surface', Surface, 'fc', fc), o);
        if ~ca.applied, continue; end
      end   % surface-coupled build: aligned at the product level, nothing to undo
      [ig, coh] = vdef.multilook(s_ref, s_sec, o.mlook_window);
      clear s_ref s_sec
      [Ipb, Wpb] = vdef.complexBlocks(ig, coh, ref_bin, starts, opts.block, o.coherence_threshold, zi);
      pairs(end+1,:) = [i j]; I(:,:,end+1) = Ipb; W(:,:,end+1) = Wpb; %#ok<AGROW>
      DT(end+1) = tide(j) - tide(i); DQ(end+1) = tide_q(j) - tide_q(i); DTD(end+1) = (tmid(j) - tmid(i))/86400; %#ok<AGROW>
      if zc
        % measured residual misalignment per block [ns]: GPS-predicted
        % compensation error minus what coalignment actually removed. A
        % surface-coupled build has no compensation error to predict, so
        % the 2-D control has no regressor there and is not run.
        rz_s = D.pass(j).ref_z(:).'; rz_r = D.pass(i).ref_z(:).'; nn = min([numel(rz_s) numel(rz_r) Nx]);
        pred = -(rz_s(1:nn) - rz_r(1:nn)) / (C.c/2);
        if isfield(ca, 'dtau_profile') && numel(ca.dtau_profile) >= 2
          pf = ca.dtau_profile(:).'; rmv = interp1(linspace(1, nn, numel(pf)), pf, 1:nn, 'linear', 'extrap');
        else
          rmv = repmat(ca.dtau_bulk, 1, nn);
        end
        res = (pred - rmv) * 1e9;
        DRES(:,end+1) = arrayfun(@(b) mean(res(starts(b):min(starts(b)+opts.block-1, nn)), 'omitnan'), 1:nb).'; %#ok<AGROW>
      end
    end
  end
  clear D
  I = I(:,:,2:end); W = W(:,:,2:end); npair = size(pairs, 1);
  fprintf('  %d pairs stacked in %.0f s\n', npair, toc(t1));
  if npair < 8, fprintf('  too few pairs\n'); continue; end

  R = vdef.tidalStack(I, W, DT, kz, struct('dres', DRES, 'pairs', pairs, 'dt', DTD, 'dq', DQ));
  % nothing at or below the ice base (minus the margin) is ice: mask it
  base_b = arrayfun(@(s) median(base_col(s:min(s+opts.block-1, Nx)), 'omitnan'), starts);
  below = zsel(:) > base_b(:).' - opts.base_margin;
  for fld = {'a', 'a_sd', 'a_sd_boot', 'a_sd_jk', 'F', 'a_t', 'a_t_sd', 'v_t', 'F_t', 'a_q', 'b_q', 'v_q', 'F_q', 'a_q_sd', 'b_q_sd', 'a2', 'g2', 'F2'}
    if isfield(R, fld{1}) && isequal(size(R.(fld{1})), size(below)), R.(fld{1})(below) = NaN; end
  end
  fprintf('  ice base per block %.0f-%.0f m; %d of %d cells masked below base - %d m\n', min(base_b), max(base_b), ...
    nnz(below), numel(below), opts.base_margin);
  if ~isfield(R, 'a2')   % no 2-D control (surface-coupled build)
    R.a2 = nan(size(R.a)); R.g2 = nan(size(R.a)); R.F2 = nan(size(R.a)); R.rcol = nan(nb, 1);
  end
  fprintf('  self-test: worst injected-shift error %.2e mm/m (%d cells at the grid edge)\n', 1e3*R.inj_err, R.n_edge);
  if R.inj_err > 1e-4, fprintf('  *** SELF-TEST FAIL ***\n'); end

  % seaward frame from the flexure fit, if that has been run
  x_sea = xb - xb(1);
  ff = fullfile(opts.out_dir, 'flexure_fit_cats.mat');
  if exist(ff, 'file')
    F = load(ff); F = F.(char(fieldnames(F)));
    pf = pn;
    ii = find(strcmp({F.lines.name}, pf), 1);
    if ~isempty(ii), x_sea = F.lines(ii).x_flip_sign * (xb - F.lines(ii).x_flip_ref); end
  end

  % the unwrapped chain's own slope at 100 m, same pairs, for comparison
  c_unw = nan(nb, 1);
  if exist(opts.net_dir, 'dir')
    c_unw = unwrapped_slope(pn, opts.net_dir, pass_ok, tide, nb, o.max_baseline);
  end

  [~, i100] = min(abs(zsel - 100));
  fprintf('  at %.0f m, all mm per m of tide; corr(dt, dtide) over the pairs %.2f, corr(dq, dtide) %.2f:\n', ...
    zsel(i100), R.r_tt, R.r_qt);
  fprintf('  %6s | %8s %6s %6s %5s | %8s %6s %7s | %8s %6s %8s %6s | %8s %7s %6s | %8s\n', 'x km', 'a_stack', 'sd_jk', 'sd_bt', 'F', ...
    'a_trend', 'sd', 'mm/day', 'a_q', 'sd', 'b_q', 'sd', 'a_ctrl', 'gain', 'r(d,T)', 'unwrap');
  [~, ord] = sort(x_sea);
  for b = ord(:).'
    fprintf('  %6.2f | %+8.2f %6.2f %6.2f %5.2f | %+8.2f %6.2f %+7.2f | %+8.2f %6.2f %+8.2f %6.2f | %+8.2f %+7.2f %6.2f | %+8.2f\n', x_sea(b)/1e3, ...
      1e3*R.a(i100,b), 1e3*R.a_sd(i100,b), 1e3*R.a_sd_boot(i100,b), R.F(i100,b), ...
      1e3*R.a_t(i100,b), 1e3*R.a_t_sd(i100,b), 1e3*R.v_t(i100,b), ...
      1e3*R.a_q(i100,b), 1e3*R.a_q_sd(i100,b), 1e3*R.b_q(i100,b), 1e3*R.b_q_sd(i100,b), ...
      1e3*R.a2(i100,b), 1e3*R.g2(i100,b), R.rcol(b), 1e3*c_unw(b));
  end
  OUT(end+1) = struct('name', pn, 'x_sea', x_sea, 'along', xb, 'lat', latb, 'lon', lonb, 'zsel', zsel, ...
    'base', base_b(:), 'base_col', base_col(:), ...
    'a', R.a, 'a_sd', R.a_sd, 'a_sd_boot', R.a_sd_boot, 'F', R.F, ...
    'a_t', R.a_t, 'a_t_sd', R.a_t_sd, 'v_t', R.v_t, 'F_t', R.F_t, 'r_tt', R.r_tt, ...
    'a_q', R.a_q, 'b_q', R.b_q, 'v_q', R.v_q, 'F_q', R.F_q, 'a_q_sd', R.a_q_sd, 'b_q_sd', R.b_q_sd, 'r_qt', R.r_qt, ...
    'a2', R.a2, 'g2', R.g2, 'F2', R.F2, 'rcol', R.rcol, 'c_unw', c_unw, 'npair', npair, ...
    'pairs', pairs, 'dt_days', DTD, 'kz', kz, 'inj_err', R.inj_err, 'main_idx', main_idx); %#ok<AGROW>
end
fprintf(['\na_stack: coherent 1-D scan on the tide alone. sd_jk: delete-one-PASS jackknife, the error bar\n' ...
         'to quote (a pass''s own error enters every pair it is in); sd_bt: the pair bootstrap, which\n' ...
         'treats those pairs as independent and understates it. F: peak normalised coherent sum (1 =\n' ...
         'every pair phase-locked). a_trend/mm/day: 2-D scan with the pair interval as second regressor\n' ...
         '- the secular rate co-estimated so it cannot leak into the tide through corr(dt, dtide).\n' ...
         'a_q/b_q: in-phase and quadrature-tide response with the trend (b_q = 0: elastic, strain rate\n' ...
         'a quarter period ahead of the tide).\n' ...
         'a_ctrl/gain: 2-D scan with the measured residual misalignment (compensated builds only; gain\n' ...
         'in mm/m per ns); r(d,T) its collinearity with the tide - above ~0.9 a_ctrl is a ridge.\n' ...
         'unwrap: the unwrapped chain, same pairs, slope through the origin.\n']);
if ~exist(opts.out_dir, 'dir'), mkdir(opts.out_dir); end
fn_out = fullfile(opts.out_dir, sprintf('tidal_stack%s.mat', opts.tag));
save(fn_out, 'OUT');
fprintf('saved %s\n', fn_out);
end

%% ========================================================================
function drop = resolve_drop(spec, pn, out_dir, qfloor)
% Turn opts.drop_passes into indices for this product.
drop = [];
if isempty(spec), return; end
if ischar(spec) || isstring(spec)
  if ~strcmpi(char(spec), 'auto')
    error('tidal_stack:drop', 'drop_passes must be [], ''auto'' or a struct');
  end
  f = fullfile(out_dir, 'pass_quality.mat');
  if ~exist(f, 'file')
    warning('tidal_stack:drop', ['drop_passes = auto needs %s; run ' ...
      'scripts/diagnostics/pass_quality.m first. Dropping nothing.'], f);
    return;
  end
  S = load(f); Q = S.OUT;
  k = find(strcmp({Q.name}, pn), 1);
  if isempty(k), return; end
  bad = (Q(k).q_med < qfloor) | (Q(k).frac_applied < 0.5);
  drop = find(bad);
  return;
end
k = find(strcmp({spec.name}, pn), 1);
if ~isempty(k), drop = spec(k).idx(:).'; end
end

%% ========================================================================
function c = unwrapped_slope(pn, net_dir, pass_ok, tide, nb, max_bl)
c = nan(nb, 1);
f = dir(fullfile(net_dir, [pn '_vvel_*.mat']));
keep = ~cellfun('isempty', regexp({f.name}, ['^' regexptranslate('escape', pn) '_vvel_\d+_\d+\.mat$'], 'once'));
f = f(keep); DH = nan(nb, 0); DT = [];
for q = 1:numel(f)
  ws = warning('off', 'MATLAB:load:variableNotFound');   % `aligned` is absent from older outputs
  o = load(fullfile(net_dir, f(q).name), 'dh_blk', 'depth_blk', 'pass_idx_ref', 'pass_idx_sec', ...
    'coalign_applied', 'aligned', 'baseline_y');
  warning(ws);
  if ~vdef.pairAligned(o), continue; end
  if max(abs(o.baseline_y)) > max_bl, continue; end
  if ~pass_ok(o.pass_idx_ref) || ~pass_ok(o.pass_idx_sec) || size(o.dh_blk, 2) ~= nb, continue; end
  v = nan(nb, 1);
  for b = 1:nb
    d = o.depth_blk(:,b); y = o.dh_blk(:,b); g = isfinite(d) & isfinite(y);
    if nnz(g) < 10 || max(d(g)) < 100, continue; end
    v(b) = interp1(d(g), y(g), 100, 'linear', NaN);
  end
  DH(:,end+1) = v; DT(end+1) = tide(o.pass_idx_sec) - tide(o.pass_idx_ref); %#ok<AGROW>
end
for b = 1:nb
  y = DH(b,:).'; g = isfinite(y);
  if nnz(g) >= 8, c(b) = DT(g).' \ y(g); end
end
end
