%RUN_LOCAL_PAIR Single-pair vvel run against a locally held multipass product.
%   Runs the real vvel_task on one pass pair of a CSARP_multipass comp_mode
%   3 product copied off the CReSIS servers, with no cluster, no param
%   spreadsheet and no gRadar - the opr_* path helpers are shadowed by the
%   stubs in ../test/stubs, which map CSARP_* products under data_root.
%
%   The products are 1-1.5 GB each, so fetch one first, e.g.
%     mkdir -p ~/data/opr/accum/2022_Antarctica_Ground/CSARP_multipass
%     scp <user>@mem1.cresis.ku.edu:/cresis/dataproducts/opr_data/accum/\
%2022_Antarctica_Ground/CSARP_multipass/EAGER_2022_GL3_multipass03.mat \
%       ~/data/opr/accum/2022_Antarctica_Ground/CSARP_multipass/
%
%   Loading the stack needs roughly 1 GB of RAM per 4250 x 1930 x 3 passes,
%   so a 13-pass product wants ~4 GB free. Unlike the fabric module's local
%   driver this runs in Octave as well as MATLAB - nothing in the chain
%   uses an arguments block.
%
%   Outputs land in <data_root>/CSARP_vvel.

%% Paths (container defaults; falls back to host-side locations)
data_root = '/home/matlab/data/accum/2022_Antarctica_Ground';
if ~exist(data_root,'dir')
  data_root = fullfile(getenv('HOME'),'data','opr','accum','2022_Antarctica_Ground');
end

this_dir  = fileparts(mfilename('fullpath'));
proj_root = fileparts(fileparts(this_dir));

addpath(proj_root);                                     % +vdef
addpath(fullfile(proj_root,'opr_vvel'));                % vvel, vvel_task, ...
addpath(fullfile(proj_root,'opr_vvel','test','stubs')); % shadow opr_* helpers

%% Which product and which pair
% EAGER_2022_GL3: 13 passes, 2022-12-09 to 2022-12-12, baseline main 11
% (20221211_09). Pass 3 is 20221210_05 and pass 10 is 20221211_07, so
% [11 3] is a ~1.5 day baseline and [11 10] is ~5 hours.
param = [];
param.radar_name    = 'accum3';
param.season_name   = '2022_Antarctica_Ground';
param.day_seg       = '20221211_09';   % only used to resolve season paths
param.opr_file_lock = false;
param.stub_out_root = data_root;
param.load.pair     = [11 10];         % [ref sec] indices into pass

pv = [];
pv.pass_name   = 'EAGER_2022_GL3';
pv.in_path     = 'multipass';
pv.out_path    = 'vvel';
pv.out_file_exts = {'.png'};

pv.mlook_window        = [5 15];
pv.coherence_threshold = 0.3;
pv.max_gap_bins        = 20;
pv.min_coverage        = 0.3;
pv.ref_twtt_offset     = 50e-9;
pv.block_size          = 200;          % 500 m at 2.5 m along-track sampling

pv.order         = 2;
pv.fit_top_depth = 20;
pv.fit_bot_depth = 400;
pv.norm_depth    = 400;
pv.min_samples   = 50;

pv.vdef = struct('rho_sfc', 0.35, 'rho_bco', 0.81, 'bco_depth', 60);
param.vvel = pv;

%% Run the one pair (vvel_task loads the product itself when mp is omitted)
fprintf('=== vvel_task: %s pair [%d %d] ===\n', ...
  pv.pass_name, param.load.pair(1), param.load.pair(2));
success = vvel_task(param);
assert(success,'vvel_task failed');

fprintf('\nOutputs:\n  %s\n', fullfile(data_root,'CSARP_vvel'));
