%TEST_VVEL_TASK End-to-end test of vvel/vvel_task on a synthetic product.
%   Builds a synthetic CSARP_multipass comp_mode 3 product from a known
%   vertical strain rate using the +vdef forward model - four passes at
%   realistic McMurdo repeat intervals, ONE OF THEM DISABLED so that the
%   pass index and the data-slice index genuinely differ - then runs the
%   real vvel.m (with the OPR support functions stubbed) and checks that
%   each pair recovers the truth.
%
%   What this covers that scripts/synthetic_vertical_velocity.m does not:
%     - pass index vs data slice index under a pass_en_mask with a gap
%     - the per-column repeat interval, recovered from each pass's own
%       gps_time/along_track (multipass does not resample gps_time), and
%       the per-block rescaling of velocity that follows from it
%     - the spatial baseline diagnostics
%     - the product/figure/output conventions of the OPR adapter
%
%   The truth strain rate is set well above firn-compaction rates so that
%   both the 5 h and the 3.2 day baseline carry an unambiguous signal in a
%   test small enough to run in a minute. Nothing in the chain depends on
%   the amplitude, so this does not change what is being tested.
%
%   Runs in MATLAB or Octave:
%     docker run --rm --platform linux/amd64 -v "$PWD":/work \
%       -w /work/opr_vvel/test gnuoctave/octave:latest \
%       octave --no-gui test_vvel_task.m

clear;
rand('seed', 11); randn('seed', 11);   %#ok<RAND> % Octave-compatible seeding
t0 = tic;

thisDir  = fileparts(mfilename('fullpath'));
projRoot = fullfile(thisDir, '..', '..');
addpath(projRoot);                    % +vdef
addpath(fullfile(thisDir, '..'));     % vvel.m, vvel_task.m, ...
addpath(fullfile(thisDir, 'stubs'));  % OPR support stubs

outRoot = fullfile(thisDir, 'output');
if exist(outRoot, 'dir'), rmdir(outRoot, 's'); end
mkdir(outRoot);

C = vdef.constants();

%% Geometry and truth
% =====================================================================
fc      = 750e6;         % accum3 centre frequency [Hz]
fs      = 300e6;         % sampling rate [Hz]
Nt      = 1500;          % 5 us record, ~420 m at firn/ice indices
Nx      = 2400;          % range lines over the 0.7 km Windless Bight line
Surface = 0.20e-6;       % surface twtt [s]

Time = (0:Nt-1).' / fs;

par = vdef.defaultParams();
par.bco_depth = 60;

H_norm  = 300;
S1_true = -5.0e-2;       % depth-averaged vertical strain rate [1/yr]
S2_true =  2.5e-2;       % top-to-bottom gradient [1/yr]
eps_zz_true = @(d) S1_true + S2_true*(2*d/H_norm - 1);

% Forward model at a repeat interval of exactly one year; dtau is linear in
% the interval, so each pass scales this by its own elapsed time
F = vdef.forwardDisplacement(eps_zz_true, Time, Surface, par, 1);
dtau_per_year = F.dtau(:,1);
dtau_per_year(~isfinite(dtau_per_year)) = 0;
depth_col = F.depth(:,1);

fprintf('Truth: S1 = %+.3e /yr, S2 = %+.3e /yr\n', S1_true, S2_true);

%% Pass timing
% =====================================================================
% Four passes over the same line. Pass 3 is walked more slowly than the
% others, so the repeat interval of the pair (1,3) varies by 1200 s along
% a 5 h baseline - the case the per-column interval logic exists for.
gps_epoch  = 1.04e9;
pass_start = [0, 2.5*3600, 5.0*3600, 3.2*86400];   % [s] after the main pass
pass_dur   = [5400, 5400,  6600,     5400];        % [s] to walk the line
pass_en    = [true, false, true,     true];        % pass 2 has no image

along_track = linspace(0, 700, Nx);                % common axis [m]

Npass = numel(pass_start);
gps_time = cell(1,Npass);
for k = 1:Npass
  gps_time{k} = gps_epoch + pass_start(k) + linspace(0, pass_dur(k), Nx);
end

% Repeat interval per column, exactly as vvel_task reconstructs it
dt_x = cell(1,Npass);
for k = 1:Npass
  dt_x{k} = gps_time{k} - gps_time{1};
end

%% Synthetic SLC stack
% =====================================================================
% Pass 1 (the main pass) is a fully developed speckle field with a depth-decaying envelope.
% Every other pass is that field delayed by its own dtau, mixed with an
% independent field so the coherence falls off with depth.
amp = 10.^(-depth_col/1200);
amp(~isfinite(amp)) = 0;

% Decay length chosen so the pairs that do NOT include the main pass still
% clear the coherence threshold at depth: two passes each correlated with
% the main pass at sqrt(gamma) are correlated with each other at gamma.
gamma = 0.95 * exp(-max(depth_col,0)/1200);
gamma(~isfinite(gamma)) = 0;

s_main = bsxfun(@times, amp, (randn(Nt,Nx) + 1i*randn(Nt,Nx))/sqrt(2));

en_idxs = find(pass_en);
data = complex(zeros(Nt, Nx, numel(en_idxs), 'single'));
for kk = 1:numel(en_idxs)
  k = en_idxs(kk);
  if k == 1
    data(:,:,kk) = single(s_main);
    continue;
  end
  % Matched-filter convention: a delay dtau multiplies by exp(-1i*2*pi*fc*dtau)
  dtau_map = dtau_per_year * (dt_x{k} / C.sec_per_year);
  carrier  = exp(-1i*2*pi*fc*dtau_map);
  indep    = bsxfun(@times, amp, (randn(Nt,Nx) + 1i*randn(Nt,Nx))/sqrt(2));
  data(:,:,kk) = single(bsxfun(@times, sqrt(gamma), carrier) .* s_main ...
    + bsxfun(@times, sqrt(max(1-gamma,0)), indep));
  clear dtau_map carrier indep;
end
clear s_main;

%% Pass struct, as combine_passes + multipass leave it
% =====================================================================
baseline_y = [0, 0.4, -0.9, 1.3];    % cross-track baseline vs the main pass [m]
baseline_z = [0, 0.2,  0.5, -0.3];   % vertical baseline vs the main pass [m]

pass = [];
for k = 1:Npass
  p = [];
  p.time        = Time;
  p.gps_time    = gps_time{k};
  p.along_track = along_track;
  p.surface     = Surface*ones(1,Nx);
  p.lat         = linspace(-77.6862, -77.6898, Nx);
  p.lon         = linspace(168.1200, 168.0963, Nx);
  p.elev        = 60*ones(1,Nx);
  p.ref_y       = baseline_y(k)*ones(1,Nx);
  p.ref_z       = baseline_z(k)*ones(1,Nx);
  w = []; w.fc = fc; w.time = Time;
  p.wfs         = w;
  p.wf          = 1;
  lay = []; lay.twtt_ref = Surface*ones(Nx,1);
  p.layers      = lay;
  p.input_type  = 'sar';
  p.param_pass  = struct('day_seg', sprintf('2022121%d_0%d', k, k));
  if k == 1
    pass = p;
  else
    pass(k) = p; %#ok<SAGROW>
  end
end

pass_name = 'windless_bight_20221209';
param_multipass = [];
param_multipass.multipass.pass_name = pass_name;
param_multipass.multipass.comp_mode = 3;
param_multipass.multipass.baseline_master_idx = 1;
param_multipass.multipass.pass_en_mask = pass_en;
param_multipass.multipass.output_fn_midfix = '';
param_combine_passes = [];
param_combine_passes.combine_passes.pass_name = pass_name;

in_dir = fullfile(outRoot, 'CSARP_multipass');
mkdir(in_dir);
in_fn = fullfile(in_dir, sprintf('%s_multipass03.mat', pass_name));
save(in_fn, '-v7', 'data','pass','param_multipass','param_combine_passes');
fprintf('Wrote synthetic product %s (%d passes, %d with images)\n', ...
  in_fn, Npass, numel(en_idxs));
clear data;

%% Assemble the param struct (what run_vvel would provide)
% =====================================================================
param = [];
param.day_seg       = '20221209_02';
param.season_name   = '2022_Antarctica_Ground';
param.radar_name    = 'accum3';
param.opr_file_lock = false;
param.stub_out_root = outRoot;

pv = [];
pv.pass_name     = pass_name;
pv.in_path       = 'multipass';
pv.out_path      = 'vvel';
pv.out_file_exts = {'.png'};
pv.fc            = [];             % exercise auto-detection from pass.wfs
pv.pairs         = 'main';
pv.mlook_window  = [5 15];
pv.block_size    = 1200;           % 2 blocks over the line
pv.coherence_threshold = 0.3;
pv.min_coverage  = 0.3;
pv.fit_top_depth = 20;
pv.fit_bot_depth = 300;
pv.norm_depth    = H_norm;
pv.vdef          = struct('bco_depth', 60);
param.vvel = pv;

%% Run
% =====================================================================
success = vvel(param);
assert(success, 'vvel did not succeed');

out_dir = fullfile(outRoot, 'CSARP_vvel');
assert(~exist(fullfile(out_dir, sprintf('%s_vvel_01_02.mat', pass_name)), 'file'), ...
  'pass 2 is disabled and must not have produced a product');

%% Check each pair against the truth
% =====================================================================
tol_S1 = 0.10;
tol_S2 = 0.20;   % the gradient is the weaker-constrained coefficient
for sec_idx = [3 4]
  out_fn = fullfile(out_dir, sprintf('%s_vvel_01_%02d.mat', pass_name, sec_idx));
  assert(exist(out_fn,'file') == 2, 'missing output %s', out_fn);
  out = load(out_fn);

  fprintf('\n--- pair 1 -> %d ---\n', sec_idx);
  fprintf('phase_sign %+d, fc %.0f MHz, baseline %.2f days\n', ...
    out.phase_sign, out.fc/1e6, out.delta_t_sec/86400);
  fprintf('%-6s %12s %12s %12s %10s\n', 'block', 'S1 [1/yr]', 'S2 [1/yr]', 'dt [h]', 'p_quad');
  for b = 1:numel(out.S1)
    fprintf('%-6d %12.4e %12.4e %12.4f %10.3g\n', ...
      b, out.S1(b), out.S2(b), out.delta_t_blk(b)*C.sec_per_year/3600, out.p_quad(b));
  end

  assert(out.phase_sign == -1, 'phase_sign is %+d, not -1', out.phase_sign);
  assert(out.fc == fc, 'fc was not read from the pass wfs struct');
  assert(out.pass_idx_ref == 1 && out.pass_idx_sec == sec_idx, ...
    'the product records the wrong pass pair');

  % Repeat interval, straight from gps_time
  expect_dt = mean(dt_x{sec_idx});
  assert(abs(out.delta_t_sec - expect_dt) < 1, ...
    'delta_t is %.1f s, expected %.1f s', out.delta_t_sec, expect_dt);

  % Per-block intervals: the reason the module reconstructs the interval
  % per column rather than taking one number for the whole line
  expect_blk = zeros(1, numel(out.delta_t_blk));
  for b = 1:numel(expect_blk)
    cidx = out.block_starts(b) : min(out.block_starts(b) + 1200 - 1, Nx);
    expect_blk(b) = mean(dt_x{sec_idx}(cidx));
  end
  got_blk = out.delta_t_blk * C.sec_per_year;
  assert(max(abs(got_blk - expect_blk)) < 1, ...
    'per-block repeat intervals are wrong by up to %.1f s', max(abs(got_blk - expect_blk)));

  % Spatial baseline, sec relative to ref (an index swap flips the sign)
  assert(max(abs(out.baseline_y - (baseline_y(sec_idx) - baseline_y(1)))) < 1e-6, ...
    'cross-track baseline is wrong');
  assert(max(abs(out.baseline_z - (baseline_z(sec_idx) - baseline_z(1)))) < 1e-6, ...
    'vertical baseline is wrong');

  % Strain rate
  assert(all(isfinite(out.S1)), 'not every block was inverted');
  S1_err = mean(out.S1) - S1_true;
  S2_err = mean(out.S2) - S2_true;
  fprintf('S1 error %+.3e /yr (%.1f%%), S2 error %+.3e /yr (%.1f%%)\n', ...
    S1_err, 100*S1_err/abs(S1_true), S2_err, 100*S2_err/abs(S2_true));
  assert(abs(S1_err) < tol_S1*abs(S1_true), ...
    'pair 1->%d: S1 not recovered within %.0f%%', sec_idx, 100*tol_S1);
  assert(abs(S2_err) < tol_S2*abs(S2_true), ...
    'pair 1->%d: S2 not recovered within %.0f%%', sec_idx, 100*tol_S2);

  % Sign convention: compression shortens the column
  deep = out.depth_blk(:,1) > 200 & isfinite(out.dh_blk(:,1));
  assert(mean(out.dh_blk(deep,1)) < 0, ...
    'compressive truth must give dh < 0 at depth');

  % dh_std is the 1-sigma of the block MEAN, not the within-block scatter:
  % it must sit below the scatter it is deflated from, and the displacement
  % signal must stand clear of it - that is what the averaging buys
  ok_std = isfinite(out.dtau_std_blk) & isfinite(out.dtau_scatter_blk);
  assert(any(ok_std(:)), 'no block carried an uncertainty at all');
  assert(all(out.dtau_std_blk(ok_std) <= out.dtau_scatter_blk(ok_std)*(1+1e-9)), ...
    'dtau_std must not exceed the within-block scatter it is deflated from');
  assert(all(out.neff_blk(ok_std) >= 1), ...
    'the effective sample count behind a block mean fell below 1');
  assert(any(out.neff_blk(ok_std) > 1.5), ...
    'the effective sample count was never above 1, so no deflation was exercised');
  deep_std = deep & isfinite(out.dh_std_blk(:,1));
  assert(abs(mean(out.dh_blk(deep_std,1))) > 5*mean(out.dh_std_blk(deep_std,1)), ...
    'the deep displacement signal must stand clear of its 1-sigma error bar');

  for ext = {'_coh.png','_dh.png','_epszz.png'}
    fig_fn = fullfile(out_dir, sprintf('%s_vvel_01_%02d%s', pass_name, sec_idx, ext{1}));
    assert(exist(fig_fn,'file') == 2, 'missing figure %s', fig_fn);
  end
end

% The pair (1,3) is the one with a genuinely varying repeat interval; make
% sure the test actually exercised that rather than passing on a constant
out13 = load(fullfile(out_dir, sprintf('%s_vvel_01_03.mat', pass_name)));
assert(out13.delta_t_spread_sec > 1000, ...
  'the 1->3 pair was meant to have a varying repeat interval (spread %.0f s)', ...
  out13.delta_t_spread_sec);

%% A pair naming a disabled pass is dropped, not processed
% =====================================================================
param_bad = param;
param_bad.vvel.pairs = [1 2];
success_bad = vvel(param_bad);
assert(~success_bad, 'a pair naming a pass with no image must not report success');
assert(~exist(fullfile(out_dir, sprintf('%s_vvel_01_02.mat', pass_name)), 'file'), ...
  'a pair naming a pass with no image must not write a product');

%% Sequential pairing reaches the pairs main pairing does not
% =====================================================================
param_seq = param;
param_seq.vvel.pairs = 'sequential';
success_seq = vvel(param_seq);
assert(success_seq, 'sequential pairing did not succeed');
assert(exist(fullfile(out_dir, sprintf('%s_vvel_03_04.mat', pass_name)), 'file') == 2, ...
  'sequential pairing did not produce the 3->4 pair');

out34 = load(fullfile(out_dir, sprintf('%s_vvel_03_04.mat', pass_name)));
expect_dt34 = mean(gps_time{4} - gps_time{3});
assert(abs(out34.delta_t_sec - expect_dt34) < 1, ...
  'the 3->4 interval is %.1f s, expected %.1f s', out34.delta_t_sec, expect_dt34);
S1_err34 = mean(out34.S1) - S1_true;
fprintf('\n3 -> 4 pair (neither is the main pass): S1 error %+.3e /yr (%.1f%%)\n', ...
  S1_err34, 100*S1_err34/abs(S1_true));
assert(abs(S1_err34) < tol_S1*abs(S1_true), ...
  'pair 3->4: S1 not recovered within %.0f%%', 100*tol_S1);

fprintf('\nPASS (%.1f s)\n', toc(t0));
