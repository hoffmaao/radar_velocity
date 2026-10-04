%RUN_MULTIPASS_SCRATCH Rerun multipass comp_mode 3 on the EAGER GL products.
%   Reruns the CReSIS +multipass processor (differential InSAR mode) over
%   the EAGER_2022_GL1..GL4 repeat-pass products, writing the
%   <name>_multipass03.mat outputs to the user's scratch season tree
%   instead of the shared dataproducts tree. multipass.m writes every
%   output next to param.multipass.fn, so the redirect works by
%   symlinking the archived combine_passes outputs
%   (mem1:/cresis/dataproducts/.../CSARP_multipass/EAGER_2022_GL*.mat)
%   into the scratch CSARP_multipass directory and pointing fn there. The
%   symlinks are created here if missing.
%
%   multipass.m is a SCRIPT, not a function - it consumes param /
%   param_override from the calling workspace and leaves its own locals
%   behind. Run ONE product per MATLAB session so no state leaks between
%   products (GL4, for instance, must NOT inherit GL3's
%   coregistration_time_shift, whose length would not even match):
%
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "product='EAGER_2022_GL3'; run('<code>/opr_vvel/server/run_multipass_scratch.m')"
%
%   The per-product settings (main pass, per-pass coregistration time
%   shifts and equalization) are copied verbatim from the comp_mode 3
%   branches of scripts/multipass_eastwind/run_multipass_EAGER.m, which
%   produced the archived Sep-Oct 2025 products. GL4 carries no
%   coregistration_time_shift there (it is commented out), and that is
%   reproduced here, not repaired.
%
%   Rerun-safe: exits early if the output multipass03.mat already exists
%   in scratch. Define force_rerun=true before running to overwrite.
%
%   BUILD WITHOUT Z-MOTION COMPENSATION (zmotion_off = true). The radar
%   rides on the ice, so the antenna-to-surface range is fixed and ref_z is
%   the ice's own motion - on this floating shelf, the tide (pass.surface is
%   constant with tide while ref_z swings by ~1 m). The standard build
%   compensates ref_z as if it were a range change, which misaligns every
%   pass in proportion to its tide; coalignment then has to undo that
%   downstream and whatever it misses stays tide-proportional. This build
%   skips the compensation at the source (multipass param
%   zmotion_comp_en = false) and writes to CSARP_multipass_nozc, leaving the
%   standard products untouched.
%
%   The frozen coregistration_time_shift and equalization vectors below
%   were estimated WITH the compensation, so they are not reused: the
%   calibration is re-derived in three stages, one fresh session each,
%   with the vectors passed between stages through
%   <product>_calib.mat rather than pasted by hand:
%     stage='coreg'     shift estimate (multipass coregistration
%                       estimation, coherent cross-correlation over the
%                       surface and ice column)
%     stage='equalize'  comp_mode 1, per-pass complex equalization: amplitude
%                       from ice-column POWER (not the toolbox's mean
%                       interferogram, which folds in coherence), phase
%                       from the toolbox's mean interferometric phase
%     stage='product'   comp_mode 3, the product the vvel chain reads
%   e.g.
%     matlab -batch "product='EAGER_2022_GL3'; zmotion_off=true; stage='coreg'; run('...')"
%   run_multipass_nozc.sh runs all three for every product.

global gRadar;

if ~exist('product','var')
  error(['Define product before running, e.g. ' ...
    'matlab -batch "product=''EAGER_2022_GL1''; run(''...run_multipass_scratch.m'')"']);
end
product = char(product);  % -batch double quotes build strings, not chars
if ~exist('force_rerun','var'), force_rerun = false; end

% master_override: rerun with a DIFFERENT main (master) pass, writing under
% a suffixed product name so the standard build is untouched:
%   matlab -batch "product='EAGER_2022_GL3'; master_override=6; run('...')"
% The frozen coregistration_time_shift and equalization vectors are kept.
% That is deliberate and first-order correct: changing the master shifts
% ref_z by ONE CONSTANT across all passes (z_masterA - z_masterB), and a
% fast-time or phase shift common to every pass cancels in any
% interferometric pair. What does NOT cancel - and what a master-change
% test therefore measures - is everything downstream that keys off the
% master itself: the resampling grid, the surface reference, and any
% master-specific residual.
if ~exist('master_override','var'), master_override = []; end
if ~exist('zmotion_off','var') || isempty(zmotion_off), zmotion_off = false; end
if ~exist('stage','var') || isempty(stage), stage = 'product'; end
stage = char(stage);
assert(ismember(stage, {'coreg','equalize','product'}), 'Unknown stage: %s', stage);
assert(zmotion_off || strcmp(stage, 'product'), ...
  'Calibration stages exist only for the zmotion_off build; the standard build uses the frozen vectors.');

scratch_dir = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_multipass';
archive_dir = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
if zmotion_off
  scratch_dir = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_multipass_nozc';
end

if ~exist(scratch_dir,'dir')
  mkdir(scratch_dir);
end

% Output name: suffixed when the master is overridden, so builds coexist
if isempty(master_override)
  product_out = product;
else
  product_out = sprintf('%s_m%02d', product, master_override);
end

% Input: symlink the archived combine_passes output into scratch so the
% processor reads the same input but writes everything here.
in_fn = fullfile(scratch_dir, [product_out '.mat']);
if ~exist(in_fn,'file')
  src_fn = fullfile(archive_dir, [product '.mat']);
  if ~exist(src_fn,'file')
    error('No combine_passes product: %s', src_fn);
  end
  [status,msg] = system(sprintf('ln -s ''%s'' ''%s''', src_fn, in_fn));
  if status ~= 0
    error('Symlink failed: %s', msg);
  end
end

out_fn = fullfile(scratch_dir, [product_out '_multipass03.mat']);
calib_fn = fullfile(scratch_dir, [product_out '_calib.mat']);
if zmotion_off && ~strcmp(stage, 'product')
  done_field = struct('coreg', 'coregistration_time_shift', 'equalize', 'equalization');
  if exist(calib_fn,'file') && ~force_rerun && ...
      ~isempty(whos('-file', calib_fn, done_field.(stage)))
    fprintf('[SKIP] %s already holds %s; define force_rerun=true to redo\n', ...
      calib_fn, done_field.(stage));
    return;
  end
elseif exist(out_fn,'file') && ~force_rerun
  fprintf('[SKIP] %s exists; define force_rerun=true to overwrite\n', out_fn);
  return;
end

%% Product settings (verbatim from run_multipass_EAGER.m, comp_mode 3)
param_override = [];
param = [];

param.multipass.fn = fullfile(scratch_dir, product_out);
param.multipass.rbins = [];
param.multipass.layer = struct('name',{'surface','bottom'},'source','layerdata','existence_check',false);
param.multipass.comp_mode = 3;
param.multipass.slope_correction_en = false;
param.multipass.output_fn_midfix = '';
param.multipass.time_gate = [];
% Stock debug plots, as in run_multipass_EAGER.m. Headless caveat: under
% matlab -batch, docking figure windows errors, and multipass.m docks
% figures both in its debug blocks AND in the unconditional plot loop
% that computes new_equalization - so debug_plots cannot dodge it. The
% four set(...,'WindowStyle','docked') calls in the server's OPR
% checkout are instead guarded with usejava('desktop') (local edit to
% mem1:~/scripts/opr, 4 Aug 2026); figures still render and save
% headless, they just do not dock.
param.multipass.debug_plots = {'debug','coherent'};

switch product
  case 'EAGER_2022'
    % NOT REBUILDABLE. The archived multipass03 product holds 13 passes,
    % but the archived combine_passes input EAGER_2022.mat now holds only
    % 10 - the input was replaced after the product was built. Any rerun
    % from the current archive produces a DIFFERENT product under the same
    % name, which is worse than no rerun. (Also the reason the 13-wide
    % pass_en_mask crashed multipass on the 10-pass input.)
    error(['EAGER_2022 cannot be rebuilt: archived combine_passes input ' ...
      'has 10 passes, archived multipass03 has 13. The input predating ' ...
      'the swap would be needed.']);
  case 'EAGER_2022_GL1'
    param.multipass.baseline_master_idx = 14;
    param.multipass.master_idx = 14;
    param.multipass.pass_en_mask = true(1,14);
    param.multipass.coregistration_time_shift = [0.7 -1.05 1 0.75 0 -1.3 -1.6 0.6 0.35 -0.3 -1.45 -0.3 0.2 0];
    param.multipass.equalization = 10.^([-5.1 -1.3 -1.6 3.7 -1.3 2.7 3.7 -1.0 4.0 0.5 3.6 -12.3 -0.2 4.6]/20) ...
      .* exp(1i*([136.2 17.4 98.8 36.3 80.9 -151.8 140.3 120.8 20.9 27.2 29.5 5.9 -151.3 0.0]/180*pi));
  case 'EAGER_2022_GL2'
    param.multipass.baseline_master_idx = 5;
    param.multipass.master_idx = 5;
    param.multipass.pass_en_mask = true(1,15);
    param.multipass.coregistration_time_shift = [-2 0.05 -1.4 -1.7 0 -0.75 -1.35 -2 0 -0.3 -1.3 -1.75 -0.7 -0.4 -0.4];
    param.multipass.equalization = 10.^([-40.6 -0.3 3.3 4.9 9.4 4.1 4.4 3.7 -5.8 5.3 2.4 3.1 1.2 3.6 1.3]/20) ...
      .* exp(1i*([-81.8 -98.3 -94.6 -137.9 0.0 45.9 -70.8 175.5 1.2 -78.3 147.4 -105.3 -73.1 -19.8 -45.2]/180*pi));
  case 'EAGER_2022_GL3'
    param.multipass.baseline_master_idx = 11;
    param.multipass.master_idx = 11;
    param.multipass.pass_en_mask = true(1,13);
    param.multipass.coregistration_time_shift = [-1 1 0.35 -0.05 -1.1 -1.45 0.5 0.2 -0.5 -1.3 0 0.3 0.3];
    param.multipass.equalization = 10.^([-1.2 2.1 -1.0 3.2 0.2 -0.1 0.4 -11.1 0.2 1.1 4.7 0.4 1.0]/20) ...
      .* exp(1i*([-89.1 65.5 -1.8 75.7 -55.1 -11.1 -178.7 -132.2 -129.9 -115.7 0.0 32.5 45.8]/180*pi));
  case 'EAGER_2022_GL4'
    param.multipass.baseline_master_idx = 9;
    param.multipass.master_idx = 9;
    param.multipass.pass_en_mask = true(1,14);
    % No coregistration_time_shift: commented out in run_multipass_EAGER.m
    param.multipass.equalization = 10.^([-2.8 4.1 3.5 8.5 3.8 -4.5 -6.6 6.5 10.0 2.1 -11.1 2.9 -7.3 -9.1]/20) ...
      .* exp(1i*([65.7 52.2 111.8 -20.5 -25.0 -126.3 -34.5 162.0 0.0 -164.0 44.9 -115.2 -49.5 -132.5]/180*pi));
  otherwise
    error('Unknown product: %s', product);
end

if ~isempty(master_override)
  fprintf('MASTER OVERRIDE: %d (product default %d)\n', ...
    master_override, param.multipass.baseline_master_idx);
  param.multipass.baseline_master_idx = master_override;
  param.multipass.master_idx = master_override;
end

%% Build without z-motion compensation: re-derive the calibration
if zmotion_off
  param.multipass.zmotion_comp_en = false;
  Np = numel(param.multipass.pass_en_mask);
  CAL = struct();
  if exist(calib_fn,'file'), CAL = load(calib_fn); end
  switch stage
    case 'coreg'
      % Search window: from just above the surface down through the
      % coherent ice column (the column ends at ~3.5 us, the base of the
      % floating shelf; 3.0 us stays inside it). Located on the archived
      % product's main pass, whose time axis the coregistered data share.
      A = load(fullfile(archive_dir, [product '_multipass03.mat']), 'pass');
      mp_ = A.pass(param.multipass.baseline_master_idx);
      dt_ = mp_.time(2) - mp_.time(1);
      sb_ = round(interp1(mp_.time, 1:numel(mp_.time), median(mp_.surface, 'omitnan')));
      param.multipass.rbins = max(1, sb_-20) : min(numel(mp_.time), sb_ + round(3.0e-6/dt_));
      clear A mp_
      param.multipass.coregistration_time_shift = zeros(1, Np);
      param.multipass.equalization = ones(1, Np);
      param.multipass.coregistration_estimation_enable = true;
      % wider than the -2:0.05:2 the compensated build needed: with no
      % compensation the expected shifts are small, but a shift at the
      % edge of the search is a failed estimate, so leave margin to see it
      param.multipass.coregistration_estimation_range = -3:0.05:3;
    case 'equalize'
      assert(isfield(CAL, 'coregistration_time_shift'), 'Run stage=''coreg'' first (%s)', calib_fn);
      param.multipass.comp_mode = 1;
      param.multipass.coregistration_time_shift = CAL.coregistration_time_shift;
      param.multipass.equalization = ones(1, Np);
    case 'product'
      assert(isfield(CAL, 'coregistration_time_shift') && isfield(CAL, 'equalization'), ...
        'Run stages coreg and equalize first (%s)', calib_fn);
      param.multipass.coregistration_time_shift = CAL.coregistration_time_shift;
      param.multipass.equalization = CAL.equalization;
  end
  fprintf('Z-MOTION COMPENSATION OFF, stage %s -> %s\n', stage, scratch_dir);
end

%% Run (mirrors the automated section of run_multipass_EAGER.m)
[~, param.multipass.pass_name] = fileparts(param.multipass.fn);
param_override = gRadar;

multipass.multipass

%% Keep what a calibration stage estimated (multipass is a script, so its
%% results are still in this workspace)
if zmotion_off && strcmp(stage, 'coreg')
  rng_ = param.multipass.coregistration_estimation_range;
  en_ = find(param.multipass.pass_en_mask);
  shift = zeros(1, Np); shift(en_) = coregistration_time_shift_est;
  at_edge = en_(abs(coregistration_time_shift_est) >= max(abs(rng_)) - 1e-9);
  if ~isempty(at_edge)
    error('Shift estimate at the search edge for pass(es) %s; widen coregistration_estimation_range', ...
      mat2str(at_edge));
  end
  CAL.coregistration_time_shift = shift;
  CAL.coreg_rbins = param.multipass.rbins;
  CAL.coreg_range = rng_;
  CAL.coreg_xcorr_sum = xcorr_sum;
  save(calib_fn, '-struct', 'CAL');
  fprintf('coregistration_time_shift (bins): %s\nsaved %s\n', mat2str(shift, 3), calib_fn);
elseif zmotion_off && strcmp(stage, 'equalize')
  % EQUALIZATION = EQUAL POWER PER PASS, PLUS A CONSTANT PHASE.
  % The toolbox's own estimate, new_equalization = mean(s_k .* conj(s_main))
  % over every bin, has amplitude gamma_k * sqrt(P_k * P_main): it is the
  % pass's gain TIMES its coherence with the main pass. Dividing by it
  % therefore boosts a decorrelated pass by 1/gamma_k - GL1's 20221209_01
  % came out at -11.2 dB, i.e. 3.6x in amplitude - and every downstream
  % stage that weights by interferogram amplitude (vdef.tidalStack weights
  % each pair by |I|) then gives that pass's pairs MORE weight, the
  % opposite of what its coherence warrants. So the amplitude here is the
  % pass's mean POWER over the ice column (10 bins below the surface down
  % to the end of the coreg window, ~3 us, where the interferograms are
  % formed), relative to the mean over enabled passes in dB, and the phase
  % is the toolbox's mean interferometric phase, which the bright surface
  % dominates. `data` is this mode-1 run's coregistered, UNequalized data,
  % still in the workspace because multipass is a script.
  en_ = find(param.multipass.pass_en_mask);
  rows_ = CAL.coreg_rbins; sb_ = rows_(1) + 20;           % the coreg window starts 20 bins above the surface
  col_rows_ = (sb_ + 10) : rows_(end);
  srf_rows_ = (sb_ - 3) : (sb_ + 3);
  Pcol_ = nan(1, Np); Psrf_ = nan(1, Np);
  for q_ = 1:numel(en_)
    s_ = data(col_rows_, :, q_); Pcol_(en_(q_)) = mean(abs(s_(isfinite(s_))).^2);
    s_ = data(srf_rows_, :, q_); Psrf_(en_(q_)) = mean(abs(s_(isfinite(s_))).^2);
  end
  rel_ = @(v) v - mean(v(en_));
  gain_dB_ = rel_(10*log10(Pcol_));
  tool_dB_ = rel_(db(new_equalization, 'voltage'));
  CAL.equalization = 10.^(gain_dB_/20) .* exp(1i*angle(new_equalization));
  CAL.equalization_gain_dB = gain_dB_;                     % ice-column power, applied
  CAL.equalization_surface_dB = rel_(10*log10(Psrf_));     % surface power, for the record
  CAL.equalization_toolbox_dB = tool_dB_;                  % gain x coherence, for the record
  CAL.equalization_phase_deg = angle(new_equalization)*180/pi;
  CAL.equalization_rows = col_rows_;
  save(calib_fn, '-struct', 'CAL');
  fprintf(['equalization gain dB (ice column, applied): %s\n' ...
           'equalization surface dB (record):           %s\n' ...
           'equalization toolbox dB (gain x coherence): %s\n' ...
           'equalization deg:                           %s\nsaved %s\n'], ...
    mat2str(gain_dB_, 3), mat2str(CAL.equalization_surface_dB, 3), mat2str(tool_dB_, 3), ...
    mat2str(CAL.equalization_phase_deg, 4), calib_fn);
end
