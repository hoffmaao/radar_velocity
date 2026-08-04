% script run_vvel
%
% Script for running vvel.m.
%
% Requires: the +vdef package on the MATLAB path, and a CSARP_multipass
% comp_mode 3 product from multipass.m (which in turn needs
% combine_passes.m).
%
% Unlike the per-frame OPR steps this is driven entirely from here rather
% than from a param spreadsheet worksheet: a repeat-pass experiment is
% identified by its pass_name and spans several segments, so there is no
% one day_seg row that owns it. param.day_seg is still set because the
% season/radar path helpers want it, but nothing is read from that segment.
%
% See also: vvel.m, vvel_task.m, run_multipass.m, multipass.m, +vdef

%% User Setup
% =====================================================================
param_override = [];

param = [];
param.radar_name  = 'accum3';
param.season_name = '2022_Antarctica_Ground';
param.day_seg     = '20221211_09';   % only used to resolve season paths
param.opr_file_lock = false;

pv = [];

% The multipass product to read. pass_name + in_path must match what
% run_multipass.m wrote; comp_mode 3 is the differential-InSAR mode, the
% only one that saves coregistered complex images per pass.
pv.pass_name = 'EAGER_2022_GL3';
pv.in_path   = 'multipass';
pv.comp_mode = 3;
pv.out_path  = 'vvel';
pv.out_file_exts = {'.jpg'};

% Which pairs. 'main' pairs every pass against the baseline main (the
% pass everything was coregistered onto, so the least resampled pair);
% 'sequential' pairs consecutive passes for the shortest baselines and the
% best coherence; 'all' does every combination; or give [ref sec] rows.
pv.pairs = 'main';

% Interferogram and traveltime estimation. The thresholds are tuned from the
% first GL3 run: the coherent ice column ends at ~3.5 us below the surface
% (~295 m, the base of the floating shelf), and where coherence thins toward
% the far end of the line a looser unwrap runs away. min_coverage stays at
% 0.3 - coverage is legitimately ~0.4 below 150 m.
pv.mlook_window        = [5 15];   % [fast-time along-track] samples
pv.coherence_threshold = 0.5;
pv.max_gap_bins        = 10;
pv.min_coverage        = 0.3;
pv.ref_twtt_offset     = 50e-9;    % dtau referenced to zero here below the surface

% Along-track blocks. The EAGER_2022 products are sampled at 2.5 m along a
% ~4.8 km line (~1930 range lines), so 200 lines is a 500 m block and about
% ten blocks per line. The OPR per-frame default of 1000 would leave two.
pv.block_size = 200;

% Inversion. Stop well above the ~295 m ice base; fitting into it produced a
% spurious sign reversal in dh below ~280 m. Fix norm_depth rather than
% letting each block use its own fitted bottom, otherwise S1/S2 are not
% comparable between blocks or pairs.
pv.order         = 2;
pv.fit_top_depth = 20;
pv.fit_bot_depth = 250;
pv.norm_depth    = 250;
pv.reg           = 0;
pv.min_samples   = 50;
pv.min_fit_range = 100;

% Firn column (see vdef.defaultParams). Only the density profile enters,
% through the refractive index that turns a traveltime difference into a
% distance, so this mainly matters over the top ~100 m.
pv.vdef = struct('rho_sfc', 0.35, 'rho_bco', 0.81, 'bco_depth', 60);

% Densification correction: leave at 0 for the hours-to-days McMurdo
% baselines. It grows linearly with the repeat interval and is only worth
% enabling season to season.
pv.densification_rate = 0;

% Flag pairs whose cross-track baseline is large enough to put topographic
% phase into the interferogram; surface referencing does not remove it.
pv.max_baseline = 10;

pv.rerun_only = false;

param.vvel = pv;

%% Automated Section
% =====================================================================

global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

success = vvel(param,param_override);
if ~success
  warning('One or more pairs did not complete; see the messages above.');
end
