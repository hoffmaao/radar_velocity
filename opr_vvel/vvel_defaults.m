function param = vvel_defaults(param)
% param = vvel_defaults(param)
%
% Fills in every param.vvel default. Called by vvel.m before the batch and
% again by vvel_task.m, so a task run standalone (a test, a one-off pair)
% gets exactly the same options as one dispatched by vvel.m. Idempotent.
%
% See also: vvel.m, vvel_task.m

% The +vdef functions read this struct directly, so the field names here
% are the option names in the package.

% pass_name: the combine_passes pass name; names the input and output files
if ~isfield(param.vvel,'pass_name') || isempty(param.vvel.pass_name)
  error('vvel:pass_name','param.vvel.pass_name is required (it names the multipass product).');
end

% in_path: CSARP_multipass product to read (combine_passes.out_path)
if ~isfield(param.vvel,'in_path') || isempty(param.vvel.in_path)
  param.vvel.in_path = 'multipass';
end

% in_fn: absolute path override for the input file (skips in_path)
if ~isfield(param.vvel,'in_fn')
  param.vvel.in_fn = '';
end

% out_path: CSARP_vvel by default. Like multipass, the products are named
% by pass_name and are NOT filed under a day_seg directory: a repeat-pass
% experiment spans segments, so no single one owns the result.
if ~isfield(param.vvel,'out_path') || isempty(param.vvel.out_path)
  param.vvel.out_path = 'vvel';
end

if ~isfield(param.vvel,'out_file_exts') || isempty(param.vvel.out_file_exts)
  param.vvel.out_file_exts = {'.jpg'};
end

% output_fn_midfix: must match param.multipass.output_fn_midfix used to
% write the input product; it is carried into the output names too
if ~isfield(param.vvel,'output_fn_midfix') || isempty(param.vvel.output_fn_midfix)
  param.vvel.output_fn_midfix = '';
end

% comp_mode: multipass comp_mode of the input file. 3 (differential InSAR)
% is the only mode that saves coregistered complex images per pass.
if ~isfield(param.vvel,'comp_mode') || isempty(param.vvel.comp_mode)
  param.vvel.comp_mode = 3;
end

% baseline_main_idx: override the main pass recorded in param_multipass.
% Leave empty; it is a diagnostic escape hatch, not a processing choice
% (the images were coregistered onto whichever pass multipass used).
if ~isfield(param.vvel,'baseline_main_idx')
  param.vvel.baseline_main_idx = [];
end

% pairs: which pass pairs to process, as [ref sec] indices into the pass
% struct array. 'main' (default) pairs every other enabled pass against
% the baseline main, which is the natural choice because that is the
% pass everything was coregistered onto; 'sequential' pairs consecutive
% enabled passes (shortest baselines, best coherence); 'all' does every
% unordered pair; or give an explicit N x 2 matrix.
if ~isfield(param.vvel,'pairs') || isempty(param.vvel.pairs)
  param.vvel.pairs = 'main';
end

% rerun_only: skip pairs whose output .mat already exists
if ~isfield(param.vvel,'rerun_only') || isempty(param.vvel.rerun_only)
  param.vvel.rerun_only = false;
end

% continue_on_error: carry on to the next pair when one fails
if ~isfield(param.vvel,'continue_on_error') || isempty(param.vvel.continue_on_error)
  param.vvel.continue_on_error = true;
end

% fc: centre frequency [Hz] converting interferogram phase to traveltime.
% Default: read from the main pass's wfs struct.
if ~isfield(param.vvel,'fc')
  param.vvel.fc = [];
end

% mlook_window: [fast-time along-track] boxcar for the interferogram and
% coherence estimate
if ~isfield(param.vvel,'mlook_window') || isempty(param.vvel.mlook_window)
  param.vvel.mlook_window = [5 15];
end

% phase_sign: sign s in dtau = s*Phi/(2*pi*fc). FIXED at -1, not estimated:
% the interferogram is sec.*conj(ref) and the matched-filter convention
% gives Phi = -2*pi*fc*dtau. Unlike the polarimetric product there is no
% coregistration row_offset field here to regress the sign against, so
% changing this is asserting a different convention, not fitting one.
if ~isfield(param.vvel,'phase_sign') || isempty(param.vvel.phase_sign)
  param.vvel.phase_sign = -1;
end

% ref_twtt_offset: two-way time below the surface return where dtau is
% referenced to zero (clear of the surface sidelobes)
if ~isfield(param.vvel,'ref_twtt_offset') || isempty(param.vvel.ref_twtt_offset)
  param.vvel.ref_twtt_offset = 50e-9;
end

% coalign_en: measure and remove the residual bulk fast-time shift between
% the two slices before the interferogram (vdef.coalignPair). multipass
% motion-compensates each pass by ref_z/(c/2); on a floating shelf ref_z
% is essentially the tide, not a platform-to-surface range change, so the
% compensation misaligns the pair in proportion to the tide, and the
% misalignment leaks into the inferred strain (measured at roughly 57 mm
% of apparent column displacement per metre of tide before this fix).
% Leave enabled except when deliberately reproducing the artefact.
if ~isfield(param.vvel,'coalign_en') || isempty(param.vvel.coalign_en)
  param.vvel.coalign_en = true;
end

% coalign_max_lag: credibility bound on the measured shift [bins]. This is
% NOT a formality - the trace-averaged surface profile carries a range
% sidelobe at about +/-5.4 bins at 0.37-0.55 of the main peak, and on one
% EAGER pair (GL1 / 20221211_07) that sidelobe outranked the true peak and
% produced a -19 ns estimate against a true 0.8 ns. The bound is what keeps
% the search inside it. 3 bins is above the tidal range over c/2 (~3 bins)
% and well above every true shift measured on these products (max 1.75
% bins). A peak found ON the bound is rejected, not clamped. Raise this
% only with evidence that the real shift is larger, and check where the
% sidelobe sits first. See vdef.coalignPair.
if ~isfield(param.vvel,'coalign_max_lag') || isempty(param.vvel.coalign_max_lag)
  param.vvel.coalign_max_lag = 3;
end

% coalign_half_win: half-width of the surface window [bins]. 60 is the
% width the estimator was validated at on the EAGER products (ratio
% 0.97-1.07 against the known residual on the uncoregistered products,
% within +/-0.9 ns of zero on the coregistered ones).
if ~isfield(param.vvel,'coalign_half_win') || isempty(param.vvel.coalign_half_win)
  param.vvel.coalign_half_win = 60;
end

% coalign_upsample: FFT upsampling factor for the cross-correlation peak.
% 32 puts the quantisation at dt/32 (0.10 ns at fs = 300 MHz), about 2% of
% a typical 5 ns residual. Upsampling rather than a parabolic fit because a
% parabola through integer lags of a narrow peak biases sub-bin shifts low.
if ~isfield(param.vvel,'coalign_upsample') || isempty(param.vvel.coalign_upsample)
  param.vvel.coalign_upsample = 32;
end

% coalign_min_quality: floor on the normalised surface-profile correlation
% (info.quality) below which the measured shift is rejected as noise and
% the pair is left unaligned.
%
% NOTE this measures something different from the phase-slope consistency
% the previous group-delay estimator reported, so the old calibration of
% this floor does not carry over. The estimator correlates TRACE-AVERAGED
% POWER, which is insensitive to interferometric coherence - a fully
% decorrelated pair over the same surface still has a well-defined
% envelope position, and measuring it is correct, not a failure. What does
% drive the quality down is a window with no coherent surface return in
% it at all, where both profiles are flat noise. Real EAGER pairs measure
% 0.778-1.000; two independent noise fields land far below. 0.5 sits
% between those populations. Without the floor, a failed measurement is
% noise over the search range and can land inside coalign_max_lag,
% applying a spurious bulk shift larger than the artefact being removed.
if ~isfield(param.vvel,'coalign_min_quality') || isempty(param.vvel.coalign_min_quality)
  param.vvel.coalign_min_quality = 0.5;
end

% coherence_threshold: samples below this do not constrain the fast-time
% unwrap and are excluded from the block averages
if ~isfield(param.vvel,'coherence_threshold') || isempty(param.vvel.coherence_threshold)
  param.vvel.coherence_threshold = 0.3;
end

% max_gap_bins: longest incoherent run in fast time that may be bridged
% during the unwrap; everything below a longer run is marked invalid
if ~isfield(param.vvel,'max_gap_bins') || isempty(param.vvel.max_gap_bins)
  param.vvel.max_gap_bins = 20;
end

% min_coverage: minimum fraction of coherent samples for a (bin, block) to
% enter the averaged dtau profile
if ~isfield(param.vvel,'min_coverage') || isempty(param.vvel.min_coverage)
  param.vvel.min_coverage = 0.3;
end

% block_size: range lines averaged into one inversion. This is where the
% precision comes from - the per-pixel signal over a few-day baseline is a
% few hundredths of a radian, well under the per-pixel phase noise.
if ~isfield(param.vvel,'block_size') || isempty(param.vvel.block_size)
  param.vvel.block_size = 1000;
end

% order: Legendre order of the velocity fit. 2 gives a strain rate linear
% in depth (S1 mean, S2 gradient) and matches the EGIG 2011-2012 products.
if ~isfield(param.vvel,'order') || isempty(param.vvel.order)
  param.vvel.order = 2;
end

% fit_top_depth / fit_bot_depth: depth range included in the inversion [m].
% The top excludes the surface-return neighbourhood; an empty bottom takes
% the deepest valid sample per block.
if ~isfield(param.vvel,'fit_top_depth') || isempty(param.vvel.fit_top_depth)
  param.vvel.fit_top_depth = 20;
end
if ~isfield(param.vvel,'fit_bot_depth')
  param.vvel.fit_bot_depth = [];
end

% norm_depth: H normalising depth for the Legendre basis [m]. Empty uses
% each block's own fitted bottom, which makes S1/S2 block-dependent in
% meaning; set it explicitly when comparing blocks or pairs.
if ~isfield(param.vvel,'norm_depth')
  param.vvel.norm_depth = [];
end

% reg: ridge weight on the Legendre coefficients (0 = plain weighted least
% squares)
if ~isfield(param.vvel,'reg') || isempty(param.vvel.reg)
  param.vvel.reg = 0;
end

% bins_per_look: fast-time correlation length in bins, used to deflate the
% degrees of freedom. Adjacent bins of a multilooked interferogram are not
% independent; the default is the fast-time multilook length.
if ~isfield(param.vvel,'bins_per_look') || isempty(param.vvel.bins_per_look)
  param.vvel.bins_per_look = param.vvel.mlook_window(1);
end

% cols_per_look: along-track correlation length in range lines, the
% counterpart of bins_per_look. It sets how many of the columns in a block
% count as independent when the within-block scatter is turned into the
% 1-sigma of the block mean (dtau_std_blk / dh_std_blk / v_std_blk); the
% default is the along-track multilook length.
if ~isfield(param.vvel,'cols_per_look') || isempty(param.vvel.cols_per_look)
  param.vvel.cols_per_look = param.vvel.mlook_window(2);
end

% min_fit_range: a block whose VALID depth span comes out shorter than this
% [m] is dropped after the inversion. Real coherence gaps routinely leave a
% block reaching only the top few tens of metres, where the fit sees firn
% densification rather than a column strain rate - without this it lands in
% the product looking like a full-depth result. top_depth / bot_depth and
% the block_short mask survive so the reason stays visible.
if ~isfield(param.vvel,'min_fit_range') || isempty(param.vvel.min_fit_range)
  param.vvel.min_fit_range = 100;
end

% min_samples: minimum raw samples for a block to be inverted at all
if ~isfield(param.vvel,'min_samples') || isempty(param.vvel.min_samples)
  param.vvel.min_samples = 50;
end

% densification_rate: relative density per year, driving the correction for
% the refractive index rising as the firn above a reflector compacts. Off
% (0) by default: negligible over hours-to-days baselines, NOT negligible
% season to season, where it grows linearly with the repeat interval.
if ~isfield(param.vvel,'densification_rate') || isempty(param.vvel.densification_rate)
  param.vvel.densification_rate = 0;
end

% max_baseline: cross-track baseline [m] above which a pair is flagged.
% Surface referencing removes what is common to the column, but topographic
% phase from a cross-track baseline is not common to the column and does
% not come out.
if ~isfield(param.vvel,'max_baseline') || isempty(param.vvel.max_baseline)
  param.vvel.max_baseline = 10;
end

% vdef: firn column overrides (see vdef.defaultParams): rho_sfc, rho_bco,
% bco_depth, mixing, max_depth, nz. Only the density profile enters, via
% the refractive index that turns a traveltime difference into a distance,
% so this mainly matters for the top ~100 m.
if ~isfield(param.vvel,'vdef') || isempty(param.vvel.vdef)
  param.vvel.vdef = struct();
end

end
