function [success] = vvel_task(param, mp)
% [success] = vvel_task(param, mp)
%
% Processes ONE repeat-pass pair out of a multipass comp_mode 3 product
% into a vertical strain-rate profile, using the +vdef package:
%   vdef.multilook            interferogram + coherence from the two
%                             coregistered SLC slices
%   vdef.differentialRange    interferogram phase -> dtau(twtt, x),
%                             referenced to zero just below the surface
%   vdef.blockAverage         coherence-weighted along-track averaging
%   vdef.verticalDisplacement dtau -> vertical displacement relative to the
%                             surface, then relative vertical velocity
%   vdef.invertStrainRate     Legendre inversion for eps_zz(d)
% All processing options come from the param.vvel struct (see
% vvel_defaults.m, where every field is documented), which the +vdef
% functions read directly.
%
% mp is the loaded multipass product (see vvel.m, which loads the file once
% and reuses it across pairs). When it is omitted the task loads the file
% itself, so a single pair can be run standalone.
%
% param.load.pair = [ref_idx sec_idx], indices into the product's pass
% struct array (NOT into the third dimension of data - disabled passes are
% missing from data but still present in pass).
%
% Sign convention: the interferogram is sec .* conj(ref) (vdef.multilook),
% so with the matched-filter convention phase_sign = -1 and dtau is
% t_sec - t_ref. dh < 0 means the column between the surface and the
% reflector SHORTENED between the two passes (vertical compression).
%
% See also: vvel.m, run_vvel.m, multipass.m, +vdef

% er_ice, c = speed of light
physical_constants;
[output_dir,radar_type,radar_name] = opr_output_dir(param.radar_name);

C = vdef.constants();

% Idempotent: vvel.m has already done this for a dispatched batch, but a
% task run standalone must end up with the same options
param = vvel_defaults(param);

pass_name = param.vvel.pass_name;
midfix    = param.vvel.output_fn_midfix;

if ~isfield(param,'load') || ~isfield(param.load,'pair') || numel(param.load.pair) ~= 2
  error('vvel:pair','param.load.pair must be [ref_idx sec_idx].');
end
ref_idx = param.load.pair(1);
sec_idx = param.load.pair(2);

pair_id = sprintf('%s%s_%02d_%02d', pass_name, midfix, ref_idx, sec_idx);
fprintf('vvel processing pair %s (%s)\n', pair_id, datestr(now));

out_dir = opr_filename_out(param,param.vvel.out_path,'',1);
if ~exist(out_dir,'dir')
  mkdir(out_dir);
end
fn_name = sprintf('%s%s_vvel_%02d_%02d', pass_name, midfix, ref_idx, sec_idx);
out_fn  = fullfile(out_dir,[fn_name '.mat']);

%% Load the multipass product
% =========================================================================
if nargin < 2 || isempty(mp)
  mp = vvel_load_multipass(param);
  if isempty(mp)
    success = false;
    return;
  end
end
pass = mp.pass;
data = mp.data;
baseline_main_idx = mp.baseline_main_idx;
pass_en_idxs = mp.pass_en_idxs;

Npass_total = numel(pass);
if ref_idx < 1 || ref_idx > Npass_total || sec_idx < 1 || sec_idx > Npass_total
  error('vvel:pair','Pair [%d %d] is out of range for the %d passes in %s.', ...
    ref_idx, sec_idx, Npass_total, mp.in_fn);
end

% data(:,:,k) holds pass(pass_en_idxs(k)); disabled passes have no slice
k_ref = find(pass_en_idxs == ref_idx, 1);
k_sec = find(pass_en_idxs == sec_idx, 1);
if isempty(k_ref) || isempty(k_sec)
  warning('Pass %d or %d is disabled in pass_en_mask and has no image in the multipass product. Skipping pair %s.', ...
    ref_idx, sec_idx, pair_id);
  success = false;
  return;
end

%% Common axes
% =========================================================================
% multipass resamples every pass onto the baseline main's fast-time axis
% AND its along-track axis, so the main pass's own vectors ARE the common
% axes. (ref in the product is a copy of pass(baseline_main_idx) taken
% before the pass loop and may still carry a full SLC image, so it is
% deliberately not loaded.)
main_pass = pass(baseline_main_idx);

Time = main_pass.time(:);
Nt   = numel(Time);
Nx   = size(data,2);

if size(data,1) ~= Nt
  error('vvel:axes','data has %d fast-time bins but the baseline main time axis has %d.', ...
    size(data,1), Nt);
end

if isfield(main_pass,'surface') && ~isempty(main_pass.surface)
  Surface = main_pass.surface(:).';
elseif isfield(main_pass,'layers') && ~isempty(main_pass.layers) && isfield(main_pass.layers,'twtt_ref')
  Surface = main_pass.layers(1).twtt_ref(:).';
else
  error('vvel:surface','The baseline main pass carries neither .surface nor .layers(1).twtt_ref, so the surface reference cannot be placed.');
end
Surface = vvel_fit_length(Surface, Nx, 'Surface');

along_track = vvel_main_along_track(main_pass, Nx);

% Centre frequency: the pass carries it, param.vvel.fc overrides
if ~isempty(param.vvel.fc)
  fc = param.vvel.fc;
else
  fc = vvel_pass_fc(main_pass);
end

%% Repeat interval, per along-track column
% =========================================================================
% multipass resamples the images onto the main pass along-track axis but
% leaves each pass's gps_time on its own axis, so the two passes must be
% mapped through their own along_track vectors before differencing. On the
% EAGER 2022 lines a ~4.8 km line is walked in ~22 min, in either direction
% depending on the product, so the offset is close to but not exactly
% constant - spreads of 20-180 s are typical. The spread is recorded.
t_ref_x = vvel_interp_to_main(pass(ref_idx), along_track, Nx);
t_sec_x = vvel_interp_to_main(pass(sec_idx), along_track, Nx);
dt_x    = t_sec_x - t_ref_x;                       % [s]

delta_t_sec = mean(dt_x(isfinite(dt_x)));
delta_t     = delta_t_sec / C.sec_per_year;        % [yr]
delta_t_spread_sec = max(dt_x(isfinite(dt_x))) - min(dt_x(isfinite(dt_x)));

if ~isfinite(delta_t) || delta_t == 0
  error('vvel:delta_t','The repeat interval between passes %d and %d is %g s; check gps_time in the combine_passes product.', ...
    ref_idx, sec_idx, delta_t_sec);
end
fprintf('Repeat interval: %.4f days (along-track spread %.1f s)\n', ...
  delta_t_sec/86400, delta_t_spread_sec);

%% Spatial baseline diagnostics
% =========================================================================
% ref_y (cross-track) and ref_z (vertical) are already on the main pass
% along-track axis. multipass motion-compensates ref_z, but a residual
% cross-track baseline puts topographic phase into the interferogram that
% surface referencing does NOT remove, so it is reported and thresholded
% rather than silently absorbed into the strain rate.
baseline_y = vvel_baseline(pass, sec_idx, ref_idx, 'ref_y', Nx);
baseline_z = vvel_baseline(pass, sec_idx, ref_idx, 'ref_z', Nx);
max_baseline_y = max(abs(baseline_y));
if isfinite(max_baseline_y) && max_baseline_y > param.vvel.max_baseline
  warning('Cross-track baseline reaches %.1f m (threshold %.1f m) for pair %s. Topographic phase from a baseline this large is not removed by surface referencing and will bias the inferred strain rate.', ...
    max_baseline_y, param.vvel.max_baseline, pair_id);
end

%% Interferogram
% =========================================================================
opts = param.vvel;   % the +vdef functions read their options straight off it
opts.delta_t = delta_t;

s_ref = data(:,:,k_ref);
s_sec = data(:,:,k_sec);

% Coalignment: remove the residual bulk fast-time shift between the two
% slices before the interferogram. multipass motion-compensates each pass
% by ref_z/(c/2); on a floating shelf ref_z is essentially the tide, not a
% platform-to-surface range change, so the compensation misaligns the pair
% in proportion to the tide and the misalignment leaks into the inferred
% strain (~57 mm of apparent column displacement per metre of tide when
% this was left uncorrected). See vdef.coalignPair for the mechanism.
coalign = struct('dtau_bulk', NaN, 'dtau_profile', [], 'dtau_win', [], ...
  'x_win', [], 'quality_win', [], 'quality', NaN, 'peak_ratio', NaN, ...
  'n_win', 0, 'n_win_ok', 0, 'applied', false);
% multipass ADVANCES each pass by ref_z/(c/2) (spectrum times
% exp(+1i*2*pi*f*ref_z/(c/2)), multipass.m:515-522), so if none of that
% shift were absorbed downstream the secondary would trail the reference
% by MINUS the baseline: predicted = -(ref_z_sec - ref_z_ref)/(c/2).
% Measured values on the EAGER products come out at 1.2-1.3 times this,
% so slightly more than the full erroneous compensation survives.
dtau_bulk_pred = -mean(baseline_z, 'omitnan') / (c/2);
if param.vvel.coalign_en
  [s_sec, coalign] = vdef.coalignPair(s_ref, s_sec, ...
    struct('Time', Time, 'Surface', Surface, 'fc', fc), param.vvel);
  if coalign.applied
    fprintf('Coalign: mean %.3f ns, along-track range %.3f ns over %d/%d windows (predicted mean from ref_z: %.3f ns; quality %.2f, sidelobe %.2f)\n', ...
      coalign.dtau_bulk*1e9, ...
      (max(coalign.dtau_profile)-min(coalign.dtau_profile))*1e9, ...
      coalign.n_win_ok, coalign.n_win, dtau_bulk_pred*1e9, ...
      coalign.quality, coalign.peak_ratio);
  else
    warning('Coalignment failed for pair %s; proceeding on the unaligned pair. The tide-proportional artefact is NOT corrected for this pair.', pair_id);
  end
end
dtau_bulk          = coalign.dtau_bulk;
dtau_bulk_profile  = coalign.dtau_profile;
dtau_bulk_win      = coalign.dtau_win;
coalign_x_win      = coalign.x_win;
coalign_quality_win = coalign.quality_win;
coalign_quality    = coalign.quality;
coalign_peak_ratio = coalign.peak_ratio;
coalign_n_win      = coalign.n_win;
coalign_n_win_ok   = coalign.n_win_ok;
coalign_applied    = coalign.applied;

[igram, coh] = vdef.multilook(s_ref, s_sec, opts.mlook_window);
clear s_ref s_sec;

map = [];
map.Time      = Time;
map.Surface   = Surface;
map.fc        = fc;
map.phase     = angle(igram);
map.coherence = coh;
% multipass writes wrapped phase only and flattens no surface phase, so
% vdef.differentialRange does all the referencing. Why upstream flattening
% cannot be assumed: see opr_vvel/README.md.
map.phase_is_unwrapped = false;

%% Chain: dtau -> blocks -> displacement -> strain rate
% =========================================================================
[dtau, info] = vdef.differentialRange(map, opts);
fprintf('Phase sign: %+d\n', info.phase_sign);

blk  = vdef.blockAverage(dtau, map, info, opts);
Nblk = numel(blk.starts);

% Firn column model
par = vdef.defaultParams();
for fld = fieldnames(param.vvel.vdef).'
  par.(fld{1}) = param.vvel.vdef.(fld{1});
end

V = vdef.verticalDisplacement(blk, map, par, opts);

% Per-block repeat intervals. verticalDisplacement divides by the scalar
% opts.delta_t; rescale to each block's own interval, which differs from
% the line mean when the two traverses were not walked at the same speed.
% (The optional densification correction inside verticalDisplacement uses
% the scalar mean, which is far finer than its own accuracy.)
delta_t_blk = cellfun(@(cidx) mean(dt_x(cidx),'omitnan'), blk.cols) / C.sec_per_year;
V.v         = bsxfun(@rdivide, V.dh,         delta_t_blk);
V.v_std     = bsxfun(@rdivide, V.dh_std,     abs(delta_t_blk));
V.v_scatter = bsxfun(@rdivide, V.dh_scatter, abs(delta_t_blk));

S = vdef.invertStrainRate(V, blk, opts);

% Blocks whose VALID depth span came out much shorter than the requested
% fit window are not a column strain rate - a fit over the top 50 m of a
% 250 m window sees firn densification only, yet lands in the product
% looking like every other block. They are dropped here rather than in
% +vdef because this is a product-quality rule, not numerics; top_depth /
% bot_depth survive so the reason stays visible.
block_short = isfinite(S.S1) & (S.bot_depth - S.top_depth) < param.vvel.min_fit_range;
if any(block_short)
  fprintf('Dropping %d of %d block(s) whose fitted range fell below %.0f m\n', ...
    nnz(block_short), Nblk, param.vvel.min_fit_range);
  for fld = {'S1','S2','epszz_mean','p_quad','rms'}
    S.(fld{1})(block_short) = NaN;
  end
  for fld = {'coef','coef_std','eps_zz','v_fit'}
    S.(fld{1})(:,block_short) = NaN;
  end
end

fprintf('Blocks: %d, inverted: %d\n', Nblk, nnz(isfinite(S.S1)));

%% Per-block geolocation
% =========================================================================
GPS_time  = cellfun(@(cidx) mean(t_ref_x(cidx),'omitnan'), blk.cols);
Latitude  = vvel_block_mean(main_pass, 'lat',  blk.cols, Nx);
Longitude = vvel_block_mean(main_pass, 'lon',  blk.cols, Nx);
Elevation = vvel_block_mean(main_pass, 'elev', blk.cols, Nx);
Along_track = cellfun(@(cidx) mean(along_track(cidx),'omitnan'), blk.cols);

%% Plot results
% =========================================================================
h_fig = get_figures(3,true);
clear h_axes;
if exist('parula','file') || exist('parula','builtin')
  cmap = parula(256);
else
  cmap = jet(256); % Octave fallback
end

% Docking is unavailable in headless batch sessions (-nodisplay)
dock_en = exist('usejava','builtin') && usejava('desktop');

depth_axis = V.depth(:,1);
ylim_depth = [0 max(depth_axis(isfinite(depth_axis)))];
if ~all(isfinite(ylim_depth)) || diff(ylim_depth) <= 0
  ylim_depth = [0 1];
end

% Coherence image (decimated along track: the figure cannot resolve more)
dec = max(1, ceil(Nx/2000));
fig_idx = 1; clf(h_fig(fig_idx));
if dock_en, set(h_fig(fig_idx),'WindowStyle','docked'); end
% Octave's imagesc ignores 'parent' and draws into the CURRENT figure, so
% every figure is made current before its axes is populated. (set(0,...)
% rather than figure() so nothing is raised in an interactive session.)
set(0,'CurrentFigure',h_fig(fig_idx));
h_axes(fig_idx) = axes('parent',h_fig(fig_idx));
imagesc(1:dec:Nx, (Time - mean(Surface,'omitnan'))*1e6, coh(:,1:dec:Nx), ...
  'parent', h_axes(fig_idx));
set(h_axes(fig_idx),'YDir','reverse');
colormap(h_axes(fig_idx),cmap);
caxis(h_axes(fig_idx),[0 1]);
h_colorbar = colorbar(h_axes(fig_idx));
set(get(h_colorbar,'ylabel'),'string','Coherence');
title(h_axes(fig_idx),sprintf('Coherence %s (%.2f d baseline)', ...
  pair_id, delta_t_sec/86400),'Interpreter','none');
xlabel(h_axes(fig_idx),'Range line');
ylabel(h_axes(fig_idx),'TWTT below surface (\mus)');

% Block-averaged dtau and displacement
fig_idx = 2; clf(h_fig(fig_idx));
if dock_en, set(h_fig(fig_idx),'WindowStyle','docked'); end
set(0,'CurrentFigure',h_fig(fig_idx));
% axes(...,'Position') rather than subplot: subplot's 'parent' option is
% not portable to Octave, and these figures are built headless
h_axes(fig_idx) = axes('parent',h_fig(fig_idx),'Position',[0.09 0.12 0.38 0.74]);
plot(h_axes(fig_idx), 1e12*blk.dtau, V.depth);
set(h_axes(fig_idx),'YDir','reverse'); ylim(h_axes(fig_idx),ylim_depth);
% Interpreter 'none' so pair_id's underscores stay literal; the Greek
% symbols live on the axis labels, which have no underscores to escape.
% The titles stay short because the axis labels already name the quantity -
% a "Differential traveltime <pair_id>" title overruns into the right panel,
% and Octave's gnuplot backend silently drops the second line of a two-line
% title rather than wrapping it.
title(h_axes(fig_idx),pair_id,'Interpreter','none');
xlabel(h_axes(fig_idx),'\Delta\tau (ps)');
ylabel(h_axes(fig_idx),'Depth (m)');
grid(h_axes(fig_idx),'on');

h_axes(4) = axes('parent',h_fig(fig_idx),'Position',[0.58 0.12 0.38 0.74]);
plot(h_axes(4), 1e3*V.dh, V.depth);
set(h_axes(4),'YDir','reverse'); ylim(h_axes(4),ylim_depth);
title(h_axes(4),sprintf('%.2f day baseline', delta_t_sec/86400));
xlabel(h_axes(4),'\Deltah (mm)   [<0 = column shortened]');
grid(h_axes(4),'on');

% Strain-rate profiles
fig_idx = 3; clf(h_fig(fig_idx));
if dock_en, set(h_fig(fig_idx),'WindowStyle','docked'); end
set(0,'CurrentFigure',h_fig(fig_idx));
h_axes(fig_idx) = axes('parent',h_fig(fig_idx));
plot(h_axes(fig_idx), S.eps_zz, S.depth_grid, 'LineWidth', 1.5);
set(h_axes(fig_idx),'YDir','reverse'); ylim(h_axes(fig_idx),ylim_depth);
title(h_axes(fig_idx),['Vertical strain rate  ' pair_id],'Interpreter','none');
xlabel(h_axes(fig_idx),'\epsilon_{zz} (1/yr)   [<0 = compression]');
ylabel(h_axes(fig_idx),'Depth (m)');
grid(h_axes(fig_idx),'on');

%% Save outputs
% =========================================================================
eps_zz     = S.eps_zz;
depth_grid = S.depth_grid;
v_fit      = S.v_fit;
S1         = S.S1;
S2         = S.S2;
epszz_mean = S.epszz_mean;
p_quad     = S.p_quad;
coef       = S.coef;
coef_std   = S.coef_std;
fit_rms    = S.rms;
n_used     = S.n_used;
n_eff      = S.n_eff;
norm_depth = S.H;
fit_top_depth = S.top_depth;
fit_bot_depth = S.bot_depth;

dtau_blk     = blk.dtau;
% dtau_std_blk is the 1-sigma of the block MEAN; the undeflated
% within-block scatter and the effective sample count behind it are kept
% alongside it so the error bar can be audited.
dtau_std_blk     = blk.dtau_std;
dtau_scatter_blk = blk.dtau_scatter;
neff_blk         = blk.n_eff;
coh_blk      = blk.coh;
coverage_blk = blk.coverage;
block_starts = blk.starts;
depth_blk    = V.depth;
dh_blk       = V.dh;
dh_std_blk   = V.dh_std;
v_blk        = V.v;
v_std_blk    = V.v_std;
n_local      = V.n_local;
densification_applied = V.densification_applied;

phase_sign    = info.phase_sign;
max_valid_bin = info.max_valid_bin;

pass_idx_ref = ref_idx;
pass_idx_sec = sec_idx;
% Recorded, not read: multipass.m hardcodes this false, so stamping it into
% the output keeps a future toolbox version that exposes it from silently
% changing what these products mean. Rationale in opr_vvel/README.md.
surf_flatten_en = false;

param_vvel = param;
if isfield(mp,'param_multipass')
  param_multipass = mp.param_multipass;
else
  param_multipass = [];
end
if isfield(mp,'param_combine_passes')
  param_combine_passes = mp.param_combine_passes;
else
  param_combine_passes = [];
end
if isfield(param,'opr_file_lock') && param.opr_file_lock
  file_version = '1L';
else
  file_version = '1';
end
file_type = 'vvel';

fprintf('Saving output file:\n  %s\n', out_fn);
opr_save(out_fn,'eps_zz','depth_grid','v_fit','S1','S2','epszz_mean','p_quad', ...
  'coef','coef_std','fit_rms','n_used','n_eff','norm_depth', ...
  'fit_top_depth','fit_bot_depth', ...
  'dtau_blk','dtau_std_blk','dtau_scatter_blk','neff_blk', ...
  'coh_blk','coverage_blk','block_starts', ...
  'depth_blk','dh_blk','dh_std_blk','v_blk','v_std_blk','n_local', ...
  'densification_applied','phase_sign','max_valid_bin','block_short', ...
  'dtau_bulk','dtau_bulk_pred','dtau_bulk_profile','dtau_bulk_win', ...
  'coalign_x_win','coalign_quality_win','coalign_n_win','coalign_n_win_ok', ...
  'coalign_quality','coalign_peak_ratio','coalign_applied', ...
  'delta_t','delta_t_sec','delta_t_blk','delta_t_spread_sec', ...
  'baseline_y','baseline_z','fc','Time','Surface','GPS_time', ...
  'Latitude','Longitude','Elevation','Along_track', ...
  'pass_idx_ref','pass_idx_sec','baseline_main_idx','surf_flatten_en', ...
  'param_vvel','param_multipass','param_combine_passes', ...
  'file_type','file_version');

for file_ext = param.vvel.out_file_exts
  file_ext = file_ext{1};
  fprintf('Saving output images of type %s\n', file_ext);
  opr_saveas(h_fig(1),fullfile(out_dir,sprintf('%s_coh%s',fn_name,file_ext)));
  opr_saveas(h_fig(2),fullfile(out_dir,sprintf('%s_dh%s',fn_name,file_ext)));
  opr_saveas(h_fig(3),fullfile(out_dir,sprintf('%s_epszz%s',fn_name,file_ext)));
end

%% Done
% =========================================================================
fprintf('%s done %s\n', mfilename, datestr(now));

success = true;

end

%% ========================================================================
function x = vvel_fit_length(x, Nx, name)
% Trim or pad a per-column vector to the common along-track length. The
% main pass's own vectors should already be Nx long; anything else is a
% product inconsistency worth naming.
if numel(x) == Nx
  return;
end
if numel(x) > Nx
  warning('%s has %d columns but the image has %d; truncating.', name, numel(x), Nx);
  x = x(1:Nx);
else
  error('vvel:length','%s has %d columns but the image has %d.', name, numel(x), Nx);
end
end

%% ========================================================================
function along_track = vvel_main_along_track(main_pass, Nx)
% Common along-track axis. Every pass is resampled onto the main pass's, so
% its own vector is it; fall back to sample index when the field
% is absent (older combine_passes products).
if isfield(main_pass,'along_track') && numel(main_pass.along_track) >= Nx
  along_track = main_pass.along_track(:).';
  along_track = along_track(1:Nx);
else
  warning('The baseline main pass has no along_track vector; using sample index. Repeat intervals will be assigned by column number rather than by position.');
  along_track = 1:Nx;
end
end

%% ========================================================================
function fc = vvel_pass_fc(p)
% Centre frequency of a pass, from whichever of the two wfs layouts
% combine_passes wrote (a per-waveform array, or a single struct).
fc = [];
if isfield(p,'wfs') && ~isempty(p.wfs)
  if isfield(p,'wf') && ~isempty(p.wf) && numel(p.wfs) >= p.wf ...
      && isfield(p.wfs(p.wf),'fc') && ~isempty(p.wfs(p.wf).fc)
    fc = p.wfs(p.wf).fc;
  elseif isfield(p.wfs(1),'fc') && ~isempty(p.wfs(1).fc)
    fc = p.wfs(1).fc;
  end
end
if isempty(fc)
  error('vvel:fc','The centre frequency is not in the pass wfs struct; set param.vvel.fc explicitly.');
end
end

%% ========================================================================
function t = vvel_interp_to_main(p, along_track, Nx)
% Map a pass's gps_time onto the main pass along-track axis. gps_time is one
% of the fields multipass does NOT resample, so this has to be done here.
if ~isfield(p,'gps_time') || isempty(p.gps_time)
  error('vvel:gps_time','A pass has no gps_time; the repeat interval cannot be determined.');
end
g = p.gps_time(:).';
if isfield(p,'along_track') && numel(p.along_track) == numel(g)
  x = p.along_track(:).';
  % interp1 needs strictly monotonic sample points; the along-track
  % estimate can repeat a value where two range lines map to one main-pass
  % index
  [x, keep] = unique(x);
  g = g(keep);
  if numel(x) < 2
    t = repmat(mean(g), 1, Nx);
    return;
  end
  t = interp1(x, g, along_track, 'linear', 'extrap');
else
  % No position vector: assume the pass spans the same columns as the image
  t = interp1(linspace(along_track(1), along_track(end), numel(g)), g, ...
    along_track, 'linear', 'extrap');
end
t = t(:).';
end

%% ========================================================================
function b = vvel_baseline(pass, sec_idx, ref_idx, fld, Nx)
% Spatial baseline of sec relative to ref, on the main pass along-track axis
% (ref_y/ref_z ARE resampled by multipass, unlike gps_time).
b = nan(1, Nx);
if ~isfield(pass,fld)
  return;
end
a = pass(sec_idx).(fld);
c = pass(ref_idx).(fld);
if isempty(a) || isempty(c)
  return;
end
a = a(:).'; c = c(:).';
n = min([numel(a) numel(c) Nx]);
b(1:n) = a(1:n) - c(1:n);
end

%% ========================================================================
function m = vvel_block_mean(p, fld, cols, Nx)
% Per-block mean of a main-pass-axis geolocation vector, NaN when absent.
if ~isfield(p,fld) || isempty(p.(fld))
  m = nan(1, numel(cols));
  return;
end
v = p.(fld)(:).';
if numel(v) < Nx
  v(end+1:Nx) = NaN;
end
m = cellfun(@(cidx) mean(v(cidx),'omitnan'), cols);
end
