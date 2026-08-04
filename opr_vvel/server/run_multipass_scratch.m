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

global gRadar;

if ~exist('product','var')
  error(['Define product before running, e.g. ' ...
    'matlab -batch "product=''EAGER_2022_GL1''; run(''...run_multipass_scratch.m'')"']);
end
product = char(product);  % -batch double quotes build strings, not chars
if ~exist('force_rerun','var'), force_rerun = false; end

scratch_dir = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_multipass';
archive_dir = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';

if ~exist(scratch_dir,'dir')
  mkdir(scratch_dir);
end

% Input: symlink the archived combine_passes output into scratch so the
% processor reads the same input but writes everything here.
in_fn = fullfile(scratch_dir, [product '.mat']);
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

out_fn = fullfile(scratch_dir, [product '_multipass03.mat']);
if exist(out_fn,'file') && ~force_rerun
  fprintf('[SKIP] %s exists; define force_rerun=true to overwrite\n', out_fn);
  return;
end

%% Product settings (verbatim from run_multipass_EAGER.m, comp_mode 3)
param_override = [];
param = [];

param.multipass.fn = fullfile(scratch_dir, product);
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

%% Run (mirrors the automated section of run_multipass_EAGER.m)
[~, param.multipass.pass_name] = fileparts(param.multipass.fn);
param_override = gRadar;

multipass.multipass
