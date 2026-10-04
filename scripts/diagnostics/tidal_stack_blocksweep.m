function OUT = tidal_stack_blocksweep(opts)
%TIDAL_STACK_BLOCKSWEEP Is the tidal response an artefact of the block size?
%
%   THE QUESTION. Every estimate of the tidal column response in this
%   repository - the unwrapped chain, the HMC, the coherent stack - reads
%   its observations off along-track BLOCK AVERAGES. The block length is a
%   free choice (the vvel products use 200 columns, 500 m), and it is not
%   obviously innocent: the signal at 100 m is a phase of about 0.1 rad
%   per metre of tide, so anything that shifts a block's phase by that
%   much changes the answer, and a sign flip costs only 0.1 rad. This
%   function measures the sensitivity directly, by running the SAME
%   estimator on the SAME pairs at a range of block lengths.
%
%   WHY BLOCK LENGTH IS NOT INNOCENT, mechanically:
%     * vdef.complexBlocks forms an AMPLITUDE- and coherence-weighted
%       complex mean. It is a circular mean, not a linear one: across a
%       block carrying an along-track phase gradient the vector sum both
%       shrinks in modulus and pulls its phase toward the brightest
%       columns. Longer blocks span more gradient, so the bias grows with
%       block length while the noise falls. The sweep separates that
%       weighting from the averaging itself by running the whole
%       estimator a second time on UNIT phasors, where every column
%       counts alike, and reporting both answers (a_med and au_med). It
%       also reports the concentration - the modulus of the unit-phasor
%       block mean, 1 for a block whose columns all share a phase and 0
%       for a block that cancels itself - which is how much phase
%       dispersion the block is swallowing.
%     * The mis-registration artefact is separated from the signal only by
%       the part of the residual misalignment d(x) that varies ALONG
%       TRACK; in the tide direction the two are one regressor. Averaging
%       over a longer block removes exactly that along-track variation, so
%       d becomes more nearly a per-pair constant and hence MORE collinear
%       with the tide. Longer blocks should therefore make the artefact
%       LESS separable, which the returned rcol measures.
%   The two effects pull in opposite directions, so the honest statement
%   is a curve, not a number.
%
%   WHAT IS HELD FIXED. The pairs, the coalignment, the multilook, the
%   coherence threshold, the depth selection and the scan grid are all
%   identical across block lengths - the interferograms are formed ONCE
%   per pair and blocked at every length before being discarded, so the
%   only thing that varies is the averaging window. That is what makes
%   this a controlled sensitivity test rather than five separate runs.
%
%   NOISE OR BIAS? Short blocks carry less data and are noisier, so some
%   disagreement between block lengths is expected from noise alone and
%   means nothing. The sweep therefore also returns the bootstrap sd at
%   every block length, and sd_ratio - the change in the line-median
%   response between one block length and the reference, divided by the
%   bootstrap sd there. A sign flip worth reporting is one that moves the
%   estimate by more than its own error bar; a sign flip of a quantity
%   whose error bar already straddles zero is just noise being renamed.
%
%   HOW TO READ THE OUTPUT. a_med is the line median response at 100 m;
%   a_grid is the response interpolated onto a fixed seaward grid so
%   different block counts can be compared point by point; sign_agree is
%   the fraction of that grid whose sign matches the reference block
%   length. A leg whose sign_agree falls well below 1 is a leg whose
%   reported sign is a property of the processing choice, not of the ice.
%
%   opts, all optional:
%     .products   default the four calibrated legs and EAGER_2022
%     .blocks     columns per block to test (default [25 50 100 200 400 800],
%                 i.e. about 60 m to 2 km at the 2.5 m column spacing)
%     .ref_block  block length the sign is compared against (default 200)
%     .zsel, .mp_dir, .out_dir  as scripts/diagnostics/tidal_stack.m
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "addpath('<code>/scripts/diagnostics'); tidal_stack_blocksweep"
%
%   See also vdef.tidalStack, vdef.complexBlocks,
%   scripts/diagnostics/tidal_stack.m.

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
addpath(root); addpath(here); addpath(fullfile(root, 'opr_vvel'));
def = struct( ...
  'products', {vdef.surveyLines()}, ...
  'mp_dir', '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass', ...
  'blocks', [25 50 100 200 400 800], 'ref_block', 200, ...
  'zsel', 20:10:250, 'xgrid', 0:0.5:4, 'out_dir', vdef.figureDir());
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end, end
C = vdef.constants();
P = vdef.firnColumn(vdef.defaultParams());
nk = numel(opts.blocks);

OUT = struct('name', {}, 'blocks', {}, 'zsel', {}, 'a_med', {}, 'a_iqr', {}, ...
  'a_grid', {}, 'au_grid', {}, 'xgrid', {}, 'sign_agree', {}, 'rcol_med', {}, ...
  'F_med', {}, 'conc_med', {}, 'npair', {}, 'a2_med', {}, 'au_med', {}, ...
  'sd_med', {}, 'sd_grid', {}, 'sd_ratio', {}, 'sig_flip', {}, 'nb', {});
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
  kz = o.phase_sign * (4*pi*fc*nloc(zi)/C.c);
  ref_bin = round(interp1(Time, 1:Nt, Surface + o.ref_twtt_offset, 'linear', NaN));
  along = mp_.along_track(:);
  [tide, pass_ok] = pass_tide(pn, opts.mp_dir);
  Np = numel(D.pass);

  % block layout for every tested length, fixed before the pair loop
  st = cell(1,nk); xb = cell(1,nk); nbk = zeros(1,nk);
  for k = 1:nk
    st{k} = 1:opts.blocks(k):Nx; nbk(k) = numel(st{k});
    xb{k} = arrayfun(@(s) mean(along(s:min(s+opts.blocks(k)-1, Nx))), st{k}).';
  end
  fprintf('  %d passes (%d kept), %d columns; block lengths %s columns -> %s blocks\n', ...
    Np, nnz(pass_ok), Nx, mat2str(opts.blocks), mat2str(nbk));

  I = cell(1,nk); W = cell(1,nk); IU = cell(1,nk); DRES = cell(1,nk);
  for k = 1:nk, I{k} = []; W{k} = []; IU{k} = []; DRES{k} = []; end
  DT = []; npair = 0; t1 = tic;
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
      if zc   % residual misalignment, for the 2-D control (compensated build only)
        rz_s = D.pass(j).ref_z(:).'; rz_r = D.pass(i).ref_z(:).';
        nn = min([numel(rz_s) numel(rz_r) Nx]);
        pred = -(rz_s(1:nn) - rz_r(1:nn)) / (C.c/2);
        if isfield(ca, 'dtau_profile') && numel(ca.dtau_profile) >= 2
          pf = ca.dtau_profile(:).'; rem = interp1(linspace(1, nn, numel(pf)), pf, 1:nn, 'linear', 'extrap');
        else, rem = repmat(ca.dtau_bulk, 1, nn); end
        res = (pred - rem) * 1e9;
      end
      % one interferogram, every block length, then thrown away. The
      % unit-phasor copy strips the amplitude weighting without touching
      % anything else, so the two differ only in that weighting.
      igu = ig ./ max(abs(ig), realmin);
      for k = 1:nk
        [Ipb, Wpb] = vdef.complexBlocks(ig,  coh, ref_bin, st{k}, opts.blocks(k), o.coherence_threshold, zi);
        Iub        = vdef.complexBlocks(igu, coh, ref_bin, st{k}, opts.blocks(k), o.coherence_threshold, zi);
        I{k}(:,:,end+1) = Ipb; W{k}(:,:,end+1) = Wpb; IU{k}(:,:,end+1) = Iub; %#ok<AGROW>
        if zc
          DRES{k}(:,end+1) = arrayfun(@(b) mean(res(st{k}(b):min(st{k}(b)+opts.blocks(k)-1, nn)), 'omitnan'), 1:nbk(k)).'; %#ok<AGROW>
        end
      end
      clear ig coh igu
      DT(end+1) = tide(j) - tide(i); npair = npair + 1; %#ok<AGROW>
    end
  end
  clear D
  fprintf('  %d pairs in %.0f s\n', npair, toc(t1));
  if npair < 8, fprintf('  too few pairs\n'); continue; end

  % seaward frame, shared by every block length
  ff = fullfile(opts.out_dir, 'flexure_fit_cats.mat'); flip_s = 1; flip_r = along(1);
  if exist(ff, 'file')
    F = load(ff); F = F.(char(fieldnames(F)));
    pf = pn;
    ii = find(strcmp({F.lines.name}, pf), 1);
    if ~isempty(ii), flip_s = F.lines(ii).x_flip_sign; flip_r = F.lines(ii).x_flip_ref; end
  end

  [~, i100] = min(abs(zsel - 100));
  a_med = nan(1,nk); a_iqr = nan(1,nk); a2_med = nan(1,nk); rcol_med = nan(1,nk);
  F_med = nan(1,nk); conc_med = nan(1,nk); au_med = nan(1,nk); sd_med = nan(1,nk);
  a_grid = nan(numel(opts.xgrid), nk); au_grid = nan(numel(opts.xgrid), nk);
  sd_grid = nan(numel(opts.xgrid), nk);
  fprintf('  at %.0f m: %6s %5s | %8s %8s %8s %8s %8s | %7s %6s %6s\n', zsel(i100), ...
    'block', 'nb', 'a_med', 'sd_med', 'a_iqr', 'au_med', 'a2_med', '|r|med', 'F', 'conc');
  for k = 1:nk
    I{k} = I{k}(:,:,2:end); W{k} = W{k}(:,:,2:end); IU{k} = IU{k}(:,:,2:end);
    R = vdef.tidalStack(I{k}, W{k}, DT, kz, struct('dres', DRES{k}, 'self_test', 0));
    if ~isfield(R, 'a2'), R.a2 = nan(size(R.a)); R.rcol = nan(size(R.a, 2), 1); end   % surface-coupled: no 2-D control
    RU = vdef.tidalStack(IU{k}, W{k}, DT, kz, struct('self_test', 0));
    a = 1e3*R.a(i100,:); a2 = 1e3*R.a2(i100,:); au = 1e3*RU.a(i100,:);
    sd = 1e3*R.a_sd(i100,:);
    xs = flip_s * (xb{k} - flip_r) / 1e3;
    [xs_s, ord] = sort(xs); a_s = a(ord);
    g = isfinite(a_s);
    au_s = au(ord); sd_s = sd(ord);
    if nnz(g) >= 2
      a_grid(:,k)  = interp1(xs_s(g), a_s(g),  opts.xgrid(:), 'linear', NaN);
      sd_grid(:,k) = interp1(xs_s(g), sd_s(g), opts.xgrid(:), 'linear', NaN);
      gu = isfinite(au_s);
      if nnz(gu) >= 2
        au_grid(:,k) = interp1(xs_s(gu), au_s(gu), opts.xgrid(:), 'linear', NaN);
      end
    end
    a_med(k) = median(a, 'omitnan'); a2_med(k) = median(a2, 'omitnan');
    au_med(k) = median(au, 'omitnan');
    a_iqr(k) = iqr_omit(a); rcol_med(k) = median(abs(R.rcol), 'omitnan');
    F_med(k) = median(R.F(i100,:), 'omitnan');
    conc_med(k) = median(abs(reshape(IU{k}(i100,:,:), 1, [])), 'omitnan');
    sd_med(k) = median(sd, 'omitnan');
    fprintf('  %6d %5d | %+8.2f %8.2f %8.2f %+8.2f %+8.2f | %7.2f %6.2f %6.3f\n', ...
      opts.blocks(k), nbk(k), a_med(k), sd_med(k), a_iqr(k), au_med(k), a2_med(k), ...
      rcol_med(k), F_med(k), conc_med(k));
    I{k} = []; W{k} = []; IU{k} = [];
  end
  kref = find(opts.blocks == opts.ref_block, 1); if isempty(kref), kref = 1; end
  sign_agree = nan(1,nk); sd_ratio = nan(1,nk); sig_flip = nan(1,nk);
  for k = 1:nk
    g = isfinite(a_grid(:,k)) & isfinite(a_grid(:,kref));
    if ~any(g), continue; end
    sign_agree(k) = mean(sign(a_grid(g,k)) == sign(a_grid(g,kref)));
    % how far the estimate moved, in units of its own bootstrap sd
    den = sqrt(sd_grid(g,k).^2 + sd_grid(g,kref).^2);
    dz = abs(a_grid(g,k) - a_grid(g,kref)) ./ max(den, eps);
    sd_ratio(k) = median(dz, 'omitnan');
    % sign flips that are NOT explained by noise: opposite sign and both
    % estimates at least one sd clear of zero
    flip = sign(a_grid(g,k)) ~= sign(a_grid(g,kref));
    clear_k = abs(a_grid(g,k)) > sd_grid(g,k); clear_r = abs(a_grid(g,kref)) > sd_grid(g,kref);
    sig_flip(k) = mean(flip & clear_k & clear_r, 'omitnan');
  end
  fprintf('  vs the %d-column blocks on the %.1f-%.1f km grid:\n', opts.ref_block, opts.xgrid(1), opts.xgrid(end));
  fprintf('    sign agreement      %s\n', sprintf('%6.2f ', sign_agree));
  fprintf('    shift in sd units   %s\n', sprintf('%6.2f ', sd_ratio));
  fprintf('    flips beyond noise  %s\n', sprintf('%6.2f ', sig_flip));
  OUT(end+1) = struct('name', pn, 'blocks', opts.blocks, 'zsel', zsel, 'a_med', a_med, ...
    'a_iqr', a_iqr, 'a_grid', a_grid, 'au_grid', au_grid, 'xgrid', opts.xgrid, ...
    'sign_agree', sign_agree, 'rcol_med', rcol_med, 'F_med', F_med, ...
    'conc_med', conc_med, 'npair', npair, 'a2_med', a2_med, 'au_med', au_med, ...
    'sd_med', sd_med, 'sd_grid', sd_grid, 'sd_ratio', sd_ratio, 'sig_flip', sig_flip, ...
    'nb', nbk); %#ok<AGROW>
  fprintf('  %s done in %.0f s\n', nm, toc(t0));
end
if ~exist(opts.out_dir, 'dir'), mkdir(opts.out_dir); end
save(fullfile(opts.out_dir, 'tidal_stack_blocksweep.mat'), 'OUT');
fprintf('\nsaved %s\n', fullfile(opts.out_dir, 'tidal_stack_blocksweep.mat'));
end

%% ========================================================================
function v = iqr_omit(x)
x = x(isfinite(x));
if numel(x) < 4, v = NaN; return; end
q = quantile(x, [0.25 0.75]); v = q(2) - q(1);
end
