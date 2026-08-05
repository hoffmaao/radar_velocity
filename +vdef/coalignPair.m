function [s_sec, info] = coalignPair(s_ref, s_sec, map, opts)
%COALIGNPAIR Measure and remove the bulk fast-time shift between a pair.
%   [s_sec, info] = COALIGNPAIR(s_ref, s_sec, map, opts) measures the bulk
%   misalignment of the secondary slice relative to the reference slice in
%   a window around the surface return, to well under a bin, and applies
%   the inverse time shift (envelope and carrier phase) to the secondary
%   so the pair is internally aligned before the interferogram is formed.
%
%   WHY THIS EXISTS. multipass comp_mode 3 motion-compensates each pass's
%   FCS z-motion by a time shift of ref_z/(c/2) (multipass.m:512-522,
%   envelope and phase), and pass.surface is not updated for that shift
%   (layers.twtt_ref explicitly is, multipass.m:667-669). On grounded ice
%   ref_z is real platform motion and the compensation is correct. On a
%   FLOATING shelf the platform and the surface ride the tide together,
%   the antenna-to-surface range does not change, and ref_z is essentially
%   the tide - so the compensation displaces the returns by a
%   tide-proportional amount that never was a range change.
%
%   WHETHER ANY OF THAT SURVIVES INTO A GIVEN PRODUCT DEPENDS ON ITS
%   CALIBRATION. multipass's own param.multipass.coregistration_time_shift
%   applies a per-pass fast-time shift derived from the comp_mode 2
%   coregistration stage, and it removes this misalignment as a side
%   effect. Measured on the EAGER 2022 products: EAGER_2022 and GL4 carry
%   coregistration_time_shift of ALL ZEROS and show the full residual
%   (0.97-1.07 times -(ref_z_sec - ref_z_ref)/(c/2)); GL1, GL2 and GL3
%   carry nonzero shifts of up to 2 bins and are already aligned to better
%   than half a nanosecond. So this correction is a REPAIR for
%   uncoregistered products, and on a properly coregistered one it must
%   measure ~0 and do nothing. An estimator that injects noise where the
%   truth is zero corrupts the good products, which is exactly what the
%   previous group-delay estimator did (see THE ESTIMATOR below).
%
%   The correction is EMPIRICAL, from the data itself, rather than the
%   deterministic inverse: how much survives depends on the product's
%   coregistration_time_shift, which is not knowable from ref_z alone.
%
%   KNOWN LIMITATION - this removes ONE SCALAR PER PAIR, and the real
%   misalignment varies ALONG TRACK. ref_z is a per-column vector, and on
%   the EAGER pairs its 5th-to-95th-percentile spread within a single pair
%   reaches 1.1 m, i.e. 7.4 ns or over two range bins - larger than the
%   mean offset the scalar removes. The two passes traverse the line at
%   different times and sometimes in opposite directions, so the tide
%   changes by different amounts at each end. A scalar can only take out
%   the mean, and what is left is a spatially varying, tide-proportional
%   residual. Measured on the leg-1 acceptance comparison, this correction
%   cuts the build-to-build artefact from 302 to 82 microstrain per metre
%   of tide; the remainder tracks that along-track structure (per-block
%   slopes running -19 to -134). Making the estimate per-column, or per
%   along-track block, is the known next step.
%
%   THE ESTIMATOR is the normalised cross-correlation of the two slices'
%   TRACE-AVERAGED POWER profiles over a window around the surface,
%   FFT-upsampled by opts.coalign_upsample so the peak is located to a
%   fraction of a bin without a parabolic fit.
%
%   Why the trace average is well conditioned: a single trace's envelope
%   is speckle and its correlation peak is about one bin wide, which is
%   what made a parabolic peak fit underestimate sub-bin shifts several
%   fold. Averaging |s|^2 over the ~2000 along-track columns first leaves
%   the smooth pulse shape of the surface return. The surface twtt does
%   vary along track, which smears that average - but both slices are
%   resampled onto the SAME along-track axis by multipass, so the smearing
%   is common to the two profiles and cancels in their cross-correlation.
%
%   Why not the cross-spectrum group delay used before: on the EAGER
%   products it ran 1.2-1.3 times the true residual on the uncoregistered
%   product and returned up to 8.5 ns of pure noise on the coregistered
%   ones, where the truth is ~0, all at quality 0.965-1.000. The envelope
%   correlation reproduces the truth at ratio 0.97-1.07 on the former and
%   stays inside +/-0.9 ns on the latter.
%
%   Why the envelope rather than the carrier: at fc = 750 MHz one carrier
%   cycle is 1.33 ns, so the 5 ns shifts seen here span about four cycles
%   and the carrier phase alone cannot resolve them. The misalignment is a
%   genuine full time shift, so the envelope determines it, and the
%   correction is then applied to envelope and carrier together.
%
%   THE LAG BOUND IS LOAD-BEARING, not a formality. The trace-averaged
%   surface profile carries a range sidelobe of the transmit pulse at
%   about +/-5.4 bins, at 0.37-0.55 of the main peak in every EAGER pair.
%   On one pair (GL1 / 20221211_07) that sidelobe outranked the true peak
%   and produced a -19 ns estimate against a true 0.8 ns. A peak-dominance
%   test does NOT separate these cases: that pair's peak-to-sidelobe ratio
%   was 1.8 while a correctly measured pair (GL1 / 20221211_01) sat at
%   1.03. What does separate them is physics - the residual cannot exceed
%   the tidal range over c/2 (about 3 bins here) plus whatever
%   coregistration_time_shift left, so opts.coalign_max_lag defaults to
%   3 bins, comfortably above every true shift measured (max 1.75 bins)
%   and safely inside the 5.4-bin sidelobe. A peak found AT the search
%   boundary is rejected rather than clamped, since the true peak may lie
%   outside the credible range.
%
%   map fields: .Time (Nt x 1) [s], .Surface (1 x Nx) [s], .fc [Hz]
%   opts fields (defaults in opr_vvel/vvel_defaults.m):
%     .coalign_max_lag      largest credible shift [bins]; a peak at or
%                           beyond this is rejected as a failed measurement
%     .coalign_half_win     half-width of the surface window [bins]
%     .coalign_upsample     cross-correlation upsampling factor
%     .coalign_min_quality  minimum info.quality for the estimate to be
%                           applied; below it the measurement is rejected
%                           as noise and the pair is left unaligned
%
%   Returns:
%     s_sec            the secondary, advanced by the measured delay
%     info.dtau_bulk   measured extra delay of sec relative to ref [s]
%     info.quality     normalised correlation at the peak, in [0..1]
%     info.peak_ratio  sidelobe-to-peak ratio, for QC only - see above for
%                      why it must not be used as an acceptance gate
%     info.applied     false when no reliable estimate was possible; the
%                      secondary is then returned unchanged and dtau_bulk
%                      is NaN. A quality-floor rejection still reports the
%                      measured quality so the product records why.
%
%   SIGN. Positive dtau_bulk means the secondary's returns arrive LATER
%   than the reference's. In the matched-filter convention a delay tau
%   multiplies the signal by exp(-1i*2*pi*fc*tau), so the correction
%   multiplies by exp(+1i*2*pi*fc*dtau_bulk) and advances the envelope by
%   dtau_bulk.
%
%   See also vdef.multilook, vdef.differentialRange.

info = struct('dtau_bulk', NaN, 'quality', NaN, 'peak_ratio', NaN, 'applied', false);

Time = map.Time(:);
Nt = size(s_ref, 1);
dt = Time(2) - Time(1);
maxlag  = max(1, opts.coalign_max_lag);
halfwin = max(4, round(opts.coalign_half_win));
minq    = opts.coalign_min_quality;
if isfield(opts,'coalign_upsample') && ~isempty(opts.coalign_upsample)
  UP = max(1, round(opts.coalign_upsample));
else
  UP = 32;
end

sfc = mean(map.Surface, 'omitnan');
if ~isfinite(sfc)
  warning('coalignPair: no finite surface twtt; pair left unaligned.');
  return;
end
b0 = round(interp1(Time, 1:Nt, sfc, 'linear', NaN));
if ~isfinite(b0)
  warning('coalignPair: surface twtt outside the fast-time axis; pair left unaligned.');
  return;
end
win = max(1, b0-halfwin) : min(Nt, b0+halfwin);
W = numel(win);
if W < 16
  warning('coalignPair: surface window of %d bins is too short; pair left unaligned.', W);
  return;
end

% Trace-averaged power profiles, mean-removed so the correlation is not
% dominated by the DC pedestal, and tapered so the window edges do not ring
taper = 0.5 - 0.5*cos(2*pi*(0:W-1).'/(W-1));   % hann, no toolbox dependency
p_ref = trace_mean_power(s_ref(win,:), taper);
p_sec = trace_mean_power(s_sec(win,:), taper);

den = sqrt(sum(p_ref.^2) * sum(p_sec.^2));
if ~(den > 0) || ~all(isfinite([p_ref; p_sec]))
  warning('coalignPair: empty surface window; pair left unaligned.');
  return;
end

% Cross-correlation, upsampled by zero-padding the cross-spectrum
Nfft = 2^nextpow2(2*W);
X  = fft(p_sec, Nfft) .* conj(fft(p_ref, Nfft));
Xu = zeros(Nfft*UP, 1);
h  = Nfft/2;
Xu(1:h)         = X(1:h);
Xu(end-h+1:end) = X(h+1:end);
xc = real(ifft(Xu)) * UP / den;

lags = ifftshift((-Nfft*UP/2 : Nfft*UP/2-1).') / UP;
[lags, ord] = sort(lags);
xc = xc(ord);

sel = abs(lags) <= maxlag;
if ~any(sel)
  warning('coalignPair: no lags within the %.2f bin credibility bound; pair left unaligned.', maxlag);
  return;
end
lag_in = lags(sel);
xc_in  = xc(sel);
[peak, im] = max(xc_in);
dtau_bins = lag_in(im);

if ~isfinite(dtau_bins) || ~isfinite(peak)
  warning('coalignPair: cross-correlation peak is not finite; pair left unaligned.');
  return;
end

if peak < minq
  warning('coalignPair: surface correlation %.2f is below the %.2f floor; pair left unaligned.', ...
    peak, minq);
  info.quality = peak;
  return;
end

% A peak sitting on the edge of the search means the true shift may lie
% outside the credible range; that is a failed measurement, not a clamp
if abs(dtau_bins) >= maxlag - 1/UP
  warning('coalignPair: peak at %.2f bins sits on the %.2f bin credibility bound, so the true shift may be outside it; pair left unaligned.', ...
    dtau_bins, maxlag);
  return;
end

% Highest competing local maximum, reported for QC only. Measured over a
% window WIDER than the search bound: the sidelobe worth knowing about
% sits at about +/-5.4 bins, outside the 3-bin bound that deliberately
% excludes it, so a ratio computed only inside the bound would almost
% always come back NaN and tell nobody anything.
qc_sel = abs(lags) <= max(2*maxlag, 8);
info.peak_ratio = sidelobe_ratio(xc(qc_sel), lags(qc_sel), dtau_bins, peak);

dtau = dtau_bins * dt;

% Advance the secondary by the measured delay: envelope through the
% baseband spectrum, carrier phase at fc
df = 1/(Nt*dt);
f_bb = df * ifftshift(-floor(Nt/2):floor((Nt-1)/2)).';
s_sec = ifft(bsxfun(@times, fft(double(s_sec)), exp(1i*2*pi*f_bb*dtau))) ...
  * exp(1i*2*pi*map.fc*dtau);

info.dtau_bulk = dtau;
info.quality   = peak;
info.applied   = true;

end

%% ========================================================================
function p = trace_mean_power(s, taper)
% Mean power over the along-track columns, mean-removed and tapered.
% Non-finite samples are ignored rather than zeroed, so a column with a
% dropout does not pull the profile toward zero at that bin.
p = mean(abs(double(s)).^2, 2, 'omitnan');
p(~isfinite(p)) = 0;
p = (p - mean(p)) .* taper;
end

%% ========================================================================
function r = sidelobe_ratio(xc, lags, peak_lag, peak)
% Height of the tallest local maximum that is not the main peak, relative
% to the peak. Local maxima within half a bin of the peak lag are part of
% the main lobe, not competitors, so they are excluded by position rather
% than by index.
r = NaN;
if ~(peak > 0), return; end
isloc = false(size(xc));
isloc(2:end-1) = xc(2:end-1) > xc(1:end-2) & xc(2:end-1) > xc(3:end);
isloc = isloc & abs(lags - peak_lag) > 0.5;
if ~any(isloc), return; end
r = max(xc(isloc)) / peak;
end
