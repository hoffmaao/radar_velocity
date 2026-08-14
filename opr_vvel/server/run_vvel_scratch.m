%RUN_VVEL_SCRATCH Batch vertical-velocity inversion on CReSIS-server products.
%   Runs vvel.m over the EAGER 2022 repeat-pass multipass products, writing
%   CSARP_vvel outputs to the user's scratch instead of the shared season
%   tree. The input products are read in place through
%   param.vvel.in_fn (an absolute path), so unlike run_fabric_scratch.m
%   this needs no input symlinks - only the output root is redirected, by
%   the opr_* path stubs in ../test/stubs.
%
%   Rerun-safe: pairs whose output .mat already exists are skipped.
%
%   Launch on mem1 with:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('<code>/opr_vvel/server/run_vvel_scratch.m')"
%
%   THE PRODUCTS. All five were written by multipass.m comp_mode 3 in
%   Sep-Oct 2025 and hold 13-14 passes each over a ~4.8 km line crossing
%   the grounding zone at Windless Bight, spanning 2022-12-09 to
%   2022-12-12. accum3 at 750 MHz, 300 MHz sampling, 2.5 m along-track,
%   4250 fast-time bins reaching ~1020 m below the surface.
%
%   GL3 and GL4 are the clean ones: every pass sits within a few metres
%   cross-track of the main pass. EAGER_2022_GL2 contains passes with
%   cross-track baselines of tens to >100 m (20221210_06/07/08 especially),
%   which is why max_baseline is left at its default - those pairs will
%   warn, and their strain rates should not be believed without checking
%   the topographic-phase contribution first.

scratch = '/kucresis/scratch/hoffmana_sta/vvel';
mp_dir = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
season = '2022_Antarctica_Ground';

% One row per repeat-pass product: {pass_name, pairing, out_path}
product_tbl = { ...
  'EAGER_2022_GL3', 'main',     'vvel'; ...
  'EAGER_2022_GL4', 'main',     'vvel'; ...
  'EAGER_2022',     'main',     'vvel'; ...
  'EAGER_2022_GL1', 'main',     'vvel'; ...
  'EAGER_2022_GL2', 'main',     'vvel'; ...
  'EAGER_2022_GL3', 'sequential', 'vvel_seq'; ...
  'EAGER_2022_GL4', 'sequential', 'vvel_seq'};

% Restrict the run to some of the products by defining only_pass_names
% before this script runs, e.g.
%   matlab -batch "only_pass_names={'EAGER_2022_GL3'}; run('.../run_vvel_scratch.m')"
if exist('only_pass_names','var') && ~isempty(only_pass_names)
  % Compare as string, so a cell of char, a cell of string, a string array
  % and a bare char name all behave the same - a -batch launch line usually
  % ends up with double quotes, which build strings rather than char
  keep = ismember(string(product_tbl(:,1)), string(only_pass_names));
  fprintf('Restricting to %s: %d of %d table rows\n', ...
    strjoin(cellstr(string(only_pass_names)), ', '), nnz(keep), size(product_tbl,1));
  product_tbl = product_tbl(keep,:);
end

% Override the along-track block size and redirect the outputs, for a finer
% run alongside the standing one:
%   matlab -batch "block_size_override=50; out_suffix='_fine'; run('.../run_vvel_scratch.m')"
% A smaller block is a shorter along-track window - more measurements along
% the line, fewer looks in each, so the per-block strain is noisier.
if ~exist('block_size_override','var'), block_size_override = []; end
if ~exist('out_suffix','var'), out_suffix = ''; end

% pairing_override: process a different pairing than the table specifies.
% 'all' gives every unordered pair, which for 13-15 passes is 78-105 pairs
% against the 12-14 of the 'main' pairing - the redundancy a network
% inversion needs to average down error and reject inconsistent pairs.
%   matlab -batch "only_pass_names={'EAGER_2022_GL3'}; pairing_override='all'; out_suffix='_net'; run('.../run_vvel_scratch.m')"
if ~exist('pairing_override','var'), pairing_override = ''; end

this_dir  = fileparts(mfilename('fullpath'));
proj_root = fileparts(fileparts(this_dir));
addpath(proj_root);                                     % +vdef
addpath(fullfile(proj_root,'opr_vvel'));                % vvel, vvel_task, ...
addpath(fullfile(proj_root,'opr_vvel','test','stubs')); % shadow opr_* helpers

season_root = fullfile(scratch, season);
if ~exist(season_root,'dir')
  mkdir(season_root);
end

t0 = tic;
n_ok = 0; n_fail = 0;

for si = 1:size(product_tbl,1)
  [pass_name, pairing, out_path] = product_tbl{si,:};
  in_fn = fullfile(mp_dir, sprintf('%s_multipass03.mat', pass_name));
  fprintf('\n#####################################################################\n');
  fprintf('%s [%s] -> CSARP_%s%s\n', pass_name, pairing, out_path, out_suffix);
  fprintf('#####################################################################\n');

  if ~exist(in_fn,'file')
    fprintf('[FAIL] no such product: %s\n', in_fn);
    n_fail = n_fail + 1;
    continue;
  end

  param = [];
  param.radar_name    = 'accum3';
  param.season_name   = season;
  param.day_seg       = '20221211_09';  % only used to resolve season paths
  param.opr_file_lock = false;
  param.stub_out_root = season_root;

  pv = [];
  pv.pass_name   = pass_name;
  pv.in_fn       = in_fn;      % read the shared product in place
  pv.out_path    = [out_path out_suffix];
  pv.out_file_exts = {'.png'};
  if ~isempty(pairing_override), pairing = pairing_override; end
  pv.pairs       = pairing;
  pv.rerun_only  = true;

  % Tuned from the first GL3 run (3 Aug 2026). The coherence image shows the
  % coherent ice column ending at ~3.5 us below the surface (~295 m at
  % n_ice) - the ice base of the floating shelf, with the seabed return in a
  % separate band at 5.5-7 us. Fitting to 400 m reached past the base into
  % noise and produced a spurious sign reversal in dh below ~280 m, so the
  % fit stops at 250 m. The far-end blocks also ran the fast-time unwrap
  % away where coherence thins, hence the stricter coherence threshold and
  % gap limit. min_coverage stays at 0.3: coverage is legitimately ~0.4
  % below 150 m, and 0.5 truncated most blocks to the top 50-140 m.
  pv.mlook_window        = [5 15];
  pv.coherence_threshold = 0.5;
  pv.max_gap_bins        = 10;
  pv.min_coverage        = 0.3;
  pv.ref_twtt_offset     = 50e-9;
  pv.block_size          = 200;      % 500 m at 2.5 m along-track sampling
  if ~isempty(block_size_override), pv.block_size = block_size_override; end

  pv.order         = 2;
  pv.fit_top_depth = 20;
  pv.fit_bot_depth = 250;
  pv.norm_depth    = 250;
  pv.min_samples   = 50;
  pv.min_fit_range = 100;

  pv.vdef = struct('rho_sfc', 0.35, 'rho_bco', 0.81, 'bco_depth', 60);
  pv.densification_rate = 0;
  param.vvel = pv;

  t1 = tic;
  try
    ok = vvel(param);
    if ok
      n_ok = n_ok + 1;
      fprintf('[OK]   %s (%.1f min)\n', pass_name, toc(t1)/60);
    else
      n_fail = n_fail + 1;
      fprintf('[PART] %s: one or more pairs did not complete (%.1f min)\n', ...
        pass_name, toc(t1)/60);
    end
  catch ME
    n_fail = n_fail + 1;
    fprintf('[FAIL] %s: %s\n', pass_name, ME.message);
  end
  close all;
end

fprintf('\nBatch done: %d products clean, %d with failures (%.1f min)\n', ...
  n_ok, n_fail, toc(t0)/60);
fprintf('Outputs under %s\n', season_root);
